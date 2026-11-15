"""
The MCP server: tools to irradiate a material with NIS, get the results, and plot them.

Over stdio (a local client) there is no authentication. Over HTTP each request needs an API key (see auth.py) or
an OAuth token from the GitHub login (see oauth.py), and each user sees only their own runs.
"""

import functools
import logging
from contextlib import asynccontextmanager
from typing import Annotated, Any

from pydantic import Field

from mcp.server.mcpserver import Context, MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.server.mcpserver.utilities.types import Image

from . import model as m
from .auth import key_from_headers
from .materials import atoms_from_mass, check_atoms
from .plot import MAX_SERIES, plot_png
from .units import format_time, parse_time, sig

log = logging.getLogger("nismcp.server")

LOCAL_USER = "local"

INSTRUCTIONS = """\
NIS (Nuclear Inventory Simulator) calculates the nuclides made when a material is irradiated with neutrons, and
how they decay: activity (Bq), decay heat (W), atoms and mass of each nuclide against time.

How to use the tools:
1. irradiate: give the material (mass fractions of nuclides or natural elements, and the mass), the neutron
   spectrum ('14MeV' for D-T fusion, 'thermal', or a custom group spectrum), the flux, the irradiation time and
   the cooling times. It gives a run_id and a summary.
2. get_results: the nuclides at a time of the run, ranked by a quantity.
3. plot_nuclides: a PNG plot of a quantity of nuclides against time (for example the top 3 by activity).

Always tell the user the assumptions of the run (the 'assumptions' in the result), in particular the flux if
the user did not give it. The model is 0-D: the flux is the same in all the material (no self-shielding, no
transport), so for a large block the results are an upper limit near the surface."""

TimeValue = str | float


def _errors_to_user(fn):
    """An error in the input (or a missing run, or data that is not ready) goes back to the model as a tool error
    with its message, so the model can correct the call. Other errors are bugs: the SDK logs them."""
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except (ValueError, KeyError, TimeoutError, RuntimeError, PermissionError) as e:
            raise ToolError(str(e)) from e
    return wrapper


def _user(ctx: Context, store, secured: bool) -> str:
    """The user of the request: the subject of the OAuth token (github:<login>) or of the API key (its name).
    Without authentication (stdio), the one local user."""
    request = getattr(ctx.request_context, "request", None) if ctx is not None else None
    # with OAuth, the SDK checks the token and puts the user in the request
    authenticated = getattr(request, "scope", {}).get("user") if request is not None else None
    access = getattr(authenticated, "access_token", None)
    if access is not None and access.subject:
        return access.subject
    if store is not None:
        headers = ctx.headers if ctx is not None else None
        key = key_from_headers({k.lower(): v for k, v in (headers or {}).items()})
        user = store.user(key) if key else None
        if user is not None:
            return user
    if not secured:
        return LOCAL_USER
    # the middleware lets no request through without a key or a token, so this is not expected
    raise PermissionError("no valid API key or token")


def _time_rows(result: m.RunResult, tables, top_n: int, quantity: str) -> list:
    """The totals and the top nuclides at the end of each step."""
    rows = []
    end = result.irradiation_end
    for p in result.points:
        if not p.end_of_step:
            continue
        row = {"time_since_start": format_time(p.time)}
        if p.time >= end and end > 0:
            row["cooling_time"] = format_time(p.time - end)
        row.update(m.totals(tables, p.atoms))
        if top_n > 0 and p.time > 0:
            names = m.top_nuclides(tables, p.atoms, quantity, top_n)
            row[f"top_by_{quantity}"] = m.nuclide_table(tables, p.atoms, quantity, names)
        rows.append(row)
    return rows


def _assumptions(result: m.RunResult, data, flux_given: bool) -> list:
    fluxes = sorted({s.flux for s in result.steps if s.flux > 0})
    lines = [f"spectrum: {result.spectrum}"]
    if fluxes:
        text = ", ".join(f"{f:.3g}" for f in fluxes)
        lines.append(f"flux: {text} n/cm2/s" + ("" if flux_given else " (the default, the user gave no flux)"))
    lines.append("0-D model: the flux is the same in all the material; no self-shielding, no neutron transport, "
                 "no burn-up of the flux")
    if data.has_fission_yields:
        lines.append("fission makes fission products (independent yields, at the energy of the spectrum)")
    else:
        lines.append("no fission yields are loaded: a fission removes the nuclide but makes no fission products")
    lines.append("decay heat from the mean decay energies of the decay data (beta, gamma, alpha)")
    lines.append(f"nuclear data: {', '.join(f'{k}={v}' for k, v in data.describe().items() if v)}")
    return lines


def _point(result: m.RunResult, time_since_start=None, cooling_time=None, what="cooling_time") -> m.Point:
    """The point of the run at a time since the start, or at a cooling time (after the end of the irradiation)."""
    if time_since_start is not None and cooling_time is not None:
        raise ValueError("give time_since_start or cooling_time, not both")
    if cooling_time is not None:
        return result.point_at(parse_time(cooling_time, what), after=result.irradiation_end)
    if time_since_start is not None:
        return result.point_at(parse_time(time_since_start, "time_since_start"))
    return result.point_at(None)


def build_server(loader, runner: m.Runner, key_store=None, oauth=None) -> MCPServer:
    """The MCP server. key_store: the API keys (HTTP). oauth: (GitHubOAuthProvider, AuthSettings), for OAuth with
    the GitHub login (HTTP); the provider also takes the API keys."""
    secured = key_store is not None or oauth is not None
    @asynccontextmanager
    async def lifespan(_server):
        # over stdio, start to read the data when the transport has the standard streams, so a message of the
        # C++ code on stdout cannot break the protocol (start does nothing if the data is already loading)
        loader.start()
        yield None

    provider, settings = oauth if oauth is not None else (None, None)
    server = MCPServer("nis", title="NIS nuclear inventory", instructions=INSTRUCTIONS,
                       description="Neutron activation, fission products and decay with NIS", lifespan=lifespan,
                       auth_server_provider=provider, auth=settings)
    if provider is not None:
        from .oauth import CALLBACK_PATH, CONSENT_PATH, callback_handler, consent_handler
        server.custom_route(CONSENT_PATH, methods=["GET", "POST"], include_in_schema=False)(consent_handler(provider))
        server.custom_route(CALLBACK_PATH, methods=["GET"], include_in_schema=False)(callback_handler(provider))

    @server.tool(title="NIS status and nuclear data")
    @_errors_to_user
    def nuclear_data_info() -> dict[str, Any]:
        """The status of the server, the nuclear data, the spectra that have a name, and the defaults.
        Use it if a run fails, or to tell the user which data the results come from."""
        info = {"status": loader.status(), "files": loader.files.describe(),
                "spectra": {name: text for name, (_, _, text) in m.SPECTRA.items()},
                "defaults": {"flux_n_cm2_s": m.DEFAULT_FLUX, "spectrum": "14MeV"},
                "limits": {"max_steps": m.MAX_STEPS, "max_total_time": "10000 y",
                           "runs_kept_per_user": runner.max_runs, "max_series_in_a_plot": MAX_SERIES}}
        if loader.ready:
            data, tables = loader.get()
            info["nuclides_in_decay_data"] = len(tables.names)
            info["has_reaction_data"] = data.has_reactions
            info["has_fission_yields"] = data.has_fission_yields
        return info

    @server.tool(title="Irradiate a material")
    @_errors_to_user
    def irradiate(
        material: Annotated[dict[str, float] | None, Field(
            description="Mass fractions by nuclide or natural element, for example {'U235': 1.0} or "
                        "{'Fe': 0.70, 'Cr': 0.18, 'Ni': 0.12}. Fractions or percent; they are normalised. "
                        "Give mass_kg with it.")] = None,
        mass_kg: Annotated[float | None, Field(description="The mass of the material in kg.")] = None,
        atoms: Annotated[dict[str, float] | None, Field(
            description="In place of material and mass_kg: the number of atoms of each nuclide.")] = None,
        irradiation_time: Annotated[TimeValue | None, Field(
            description="How long the irradiation is, for example '1 d', '5 min', '2 y', or seconds.")] = None,
        cooling_times: Annotated[list[TimeValue] | None, Field(
            description="Times after the end of the irradiation to get results at, for example "
                        "['1 h', '1 d', '1 y'].")] = None,
        flux_n_cm2_s: Annotated[float | None, Field(
            description=f"The total neutron flux in n/cm2/s. If not given: {m.DEFAULT_FLUX:g}.")] = None,
        spectrum: Annotated[str | dict[str, list[float]] | None, Field(
            description="'14MeV' (D-T fusion, the default), 'thermal', or a group spectrum "
                        "{'bounds_ev': [n+1 ascending group bounds in eV], 'flux': [n values]}; only the shape "
                        "is used, the flux sets the size.")] = None,
        schedule: Annotated[list[dict[str, Any]] | None, Field(
            description="In place of irradiation_time and cooling_times, for pulses or changes of flux: steps "
                        "in order, [{'duration': '1 h', 'flux': 1e14}, {'duration': '1 d', 'flux': 0}]. A flux "
                        "of 0 is cooling. A step with no flux uses flux_n_cm2_s.")] = None,
        top_n: Annotated[int, Field(ge=0, le=50,
                                    description="How many nuclides to list at each step end.")] = 5,
        rank_by: Annotated[str, Field(description="Rank the nuclides by 'activity', 'decay_heat', 'atoms' or "
                                                  "'mass'.")] = "activity",
        ctx: Context = None,
    ) -> dict[str, Any]:
        """Irradiate a material with neutrons, then let it cool, and calculate the nuclides with NIS.

        It gives the total activity and decay heat, and the top nuclides, at the end of each step. Use the
        run_id with get_results and plot_nuclides. Example: 1 kg of U235 in a 14 MeV flux for 1 day, results
        1 h and 1 d after: material={'U235': 1}, mass_kg=1, irradiation_time='1 d', cooling_times=['1 h', '1 d']."""
        user = _user(ctx, key_store, secured)
        if rank_by not in m.QUANTITIES:
            raise ValueError(f"rank_by must be one of {sorted(m.QUANTITIES)}")
        data, tables = loader.get()
        notes = []
        if atoms is not None:
            if material is not None or mass_kg is not None:
                raise ValueError("give material and mass_kg, or atoms, not both")
            initial = check_atoms(atoms, tables.index)
            described = {"atoms": {k: sig(v) for k, v in initial.items()}}
        else:
            if material is None or mass_kg is None:
                raise ValueError("give material and mass_kg (or atoms)")
            initial, notes = atoms_from_mass(float(mass_kg), material, tables.index)
            described = {"material": material, "mass_kg": mass_kg}
        steps = m.steps_from_input(irradiation_time, cooling_times, schedule, flux_n_cm2_s)
        if any(s.flux > 0 for s in steps) and not data.has_reactions:
            raise ValueError("the server has no reaction data, so it can do decay only: give a flux of 0")
        spec = m.spectrum_from_input(spectrum)
        result = runner.run(user, initial, steps, spec, described, notes)
        log.info("user %s: run %s, %d steps, %.2f s", user, result.run_id, len(steps), result.solve_seconds)
        flux_given = flux_n_cm2_s is not None or bool(schedule and all("flux" in s for s in schedule))
        return {
            "run_id": result.run_id,
            "initial_atoms": {k: sig(v) for k, v in sorted(initial.items(), key=lambda kv: -kv[1])[:20]},
            "steps": [{"kind": s.kind, "duration": format_time(s.duration), "flux_n_cm2_s": sig(s.flux)}
                      for s in steps],
            "results_at_step_ends": _time_rows(result, tables, top_n, rank_by),
            "assumptions": _assumptions(result, data, flux_given),
            "notes": notes,
            "solve_time_s": sig(result.solve_seconds, 3),
        }

    @server.tool(title="Results of a run")
    @_errors_to_user
    def get_results(
        run_id: Annotated[str, Field(description="The run_id from irradiate.")],
        cooling_time: Annotated[TimeValue | None, Field(
            description="Time after the end of the irradiation, for example '1 d'. Must be one of the times "
                        "of the run (a step end).")] = None,
        time_since_start: Annotated[TimeValue | None, Field(
            description="In place of cooling_time: the time since the start of the run.")] = None,
        quantity: Annotated[str, Field(description="'activity' (Bq), 'decay_heat' (W), 'atoms' or 'mass' "
                                                   "(g).")] = "activity",
        top_n: Annotated[int, Field(ge=1, le=200, description="How many nuclides, largest first.")] = 10,
        nuclides: Annotated[list[str] | None, Field(
            description="In place of top_n: only these nuclides, for example ['Cs137', 'Sr90'].")] = None,
        ctx: Context = None,
    ) -> dict[str, Any]:
        """The nuclides of a run at one time (the last time if no time is given), with the quantity, the percent
        of the total and the half-life, and the totals."""
        user = _user(ctx, key_store, secured)
        if quantity not in m.QUANTITIES:
            raise ValueError(f"quantity must be one of {sorted(m.QUANTITIES)}")
        _, tables = loader.get()
        result = runner.get(user, run_id)
        point = _point(result, time_since_start, cooling_time)
        if nuclides:
            from .materials import canonical_nuclide
            names = [canonical_nuclide(n, tables.index) for n in nuclides]
        else:
            names = m.top_nuclides(tables, point.atoms, quantity, top_n)
        end = result.irradiation_end
        output = {"run_id": run_id, "time_since_start": format_time(point.time)}
        if end > 0 and point.time >= end:
            output["cooling_time"] = format_time(point.time - end)
        output["totals"] = m.totals(tables, point.atoms)
        output["nuclides"] = m.nuclide_table(tables, point.atoms, quantity, names)
        if quantity == "mass":
            output["note"] = "the mass of a nuclide is from its mass number (A g/mol), so it is approximate"
        return output

    @server.tool(title="Plot nuclides against time", structured_output=False)
    @_errors_to_user
    def plot_nuclides(
        run_id: Annotated[str, Field(description="The run_id from irradiate.")],
        nuclides: Annotated[list[str] | None, Field(
            description=f"The nuclides to plot, at most {MAX_SERIES}. If not given, the top_n.")] = None,
        top_n: Annotated[int, Field(ge=1, le=MAX_SERIES, description="If no nuclides are given: plot the top_n "
                                                                     "nuclides by the quantity.")] = 3,
        quantity: Annotated[str, Field(description="'activity' (Bq), 'decay_heat' (W), 'atoms' or 'mass' "
                                                   "(g).")] = "activity",
        rank_at_cooling_time: Annotated[TimeValue | None, Field(
            description="Rank the nuclides for top_n at this time after the end of the irradiation (a time of "
                        "the run). If not given: at the last time.")] = None,
        include_initial: Annotated[bool, Field(
            description="For top_n: also rank the nuclides of the material at the start (for example U235). "
                        "If false, only the products are ranked.")] = False,
        log_time: Annotated[bool | None, Field(description="A log time axis. If not given: log if the times "
                                                           "span more than a factor of 100.")] = None,
        log_value: Annotated[bool, Field(description="A log axis for the quantity.")] = True,
        ctx: Context = None,
    ) -> list:
        """A PNG plot of a quantity (activity by default) of some nuclides against the time since the start.
        The irradiation is a shaded band. Example: the top 3 products by activity: top_n=3."""
        user = _user(ctx, key_store, secured)
        _, tables = loader.get()
        result = runner.get(user, run_id)
        if nuclides:
            from .materials import canonical_nuclide
            names = [canonical_nuclide(n, tables.index) for n in nuclides]
            ranked = ""
        else:
            point = _point(result, cooling_time=rank_at_cooling_time, what="rank_at_cooling_time")
            exclude = () if include_initial else set(result.initial_atoms)
            names = m.top_nuclides(tables, point.atoms, quantity, top_n, exclude=exclude)
            ranked = f", ranked by {quantity} at {format_time(point.time)} since the start"
            if not names:
                raise ValueError(f"no nuclides with {quantity} at that time")
        png = plot_png(result, tables, names, quantity, log_time=log_time, log_value=log_value)
        return [Image(data=png, format="png"),
                f"Plot of the {quantity} of {', '.join(names)} for run {run_id}{ranked}."]

    @server.tool(title="List my runs")
    @_errors_to_user
    def list_runs(ctx: Context = None) -> dict[str, Any]:
        """The runs of this user that the server keeps in memory (the last few)."""
        user = _user(ctx, key_store, secured)
        runs = [{"run_id": r.run_id, "input": r.material, "spectrum": r.spectrum.split(":")[0],
                 "steps": [f"{s.kind} {format_time(s.duration)}" for s in r.steps]} for r in runner.list(user)]
        return {"runs": runs}

    @server.custom_route("/health", methods=["GET"])
    async def health(request):
        from starlette.responses import JSONResponse
        return JSONResponse({"status": loader.status()})

    return server
