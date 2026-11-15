"""
One run: a material in a neutron spectrum, with irradiation and cooling steps, solved with CRAM48.

Each step is cut into a few sub-steps (log spaced in time), so the results have enough points for a plot of the
activity against time. The result keeps the atoms of all the nuclides at each point.
"""

import math
import threading
import time
import uuid
from collections import OrderedDict
from dataclasses import dataclass, field

import numpy as np

from .units import format_time, parse_time, sig

# the spectra the user can give by name: (group bounds in eV, flux per group). Only the shape is used, the
# flux of each step sets the size. A narrow group at the energy, as in the reaction tests.
FUSION_ENERGY_EV = 14.1e6
THERMAL_ENERGY_EV = 0.0253
SPECTRA = {
    "14MeV": ([FUSION_ENERGY_EV * 0.9999, FUSION_ENERGY_EV * 1.0001], [1.0],
              "D-T fusion neutrons, all at 14.1 MeV (one narrow group)"),
    "thermal": ([THERMAL_ENERGY_EV * 0.9999, THERMAL_ENERGY_EV * 1.0001], [1.0],
                "thermal neutrons, all at 0.0253 eV (2200 m/s, one narrow group)"),
}
SPECTRUM_ALIASES = {"14mev": "14MeV", "14.1mev": "14MeV", "fusion": "14MeV", "dt": "14MeV", "d-t": "14MeV",
                    "thermal": "thermal", "0.0253ev": "thermal"}

DEFAULT_FLUX = 1e14            # n/cm2/s
SUB_STEPS = 12                 # points in each step, for the plots
MAX_STEPS = 50
MAX_TOTAL_TIME = 1e4 * 3.15576e7  # 10000 years
# the groups of a custom spectrum: the time to collapse the cross sections grows with the number of groups, so
# there is a limit (the large group structures have 709 or 1102 groups)
MAX_SPECTRUM_GROUPS = 2000
# the limits of the solves: there is one solve at a time, so a user can wait for the solves of other users
MAX_RUNS_IN_PROGRESS = 8            # all the users: the run that solves, and the runs that wait for it
MAX_RUNS_IN_PROGRESS_PER_USER = 2
SOLVE_WAIT_SECONDS = 120            # the longest wait for the solves of others
MAX_SOLVE_SECONDS = 300             # a run that takes longer stops (it is checked between the sub-steps)


@dataclass
class Step:
    duration: float          # s
    flux: float              # n/cm2/s, 0 for cooling

    @property
    def kind(self) -> str:
        return "irradiation" if self.flux > 0 else "cooling"


def spectrum_from_input(spectrum) -> tuple:
    """(bounds, values, description) from a name ("14MeV", "thermal") or {"bounds_ev": [...], "flux": [...]}."""
    if spectrum is None:
        spectrum = "14MeV"
    if isinstance(spectrum, str):
        key = SPECTRUM_ALIASES.get(spectrum.replace(" ", "").lower())
        if key is None:
            raise ValueError(f"unknown spectrum {spectrum!r}, use one of {sorted(SPECTRA)}, "
                             f"or give {{'bounds_ev': [...], 'flux': [...]}}")
        bounds, values, text = SPECTRA[key]
        return list(bounds), list(values), f"{key}: {text}"
    if isinstance(spectrum, dict):
        raw_bounds, raw_values = spectrum.get("bounds_ev", []), spectrum.get("flux", [])
        if not isinstance(raw_bounds, list) or not isinstance(raw_values, list):
            raise ValueError("a spectrum needs bounds_ev (n+1 group bounds in eV) and flux (n values, one per group)")
        if len(raw_values) > MAX_SPECTRUM_GROUPS:
            raise ValueError(f"a spectrum can have at most {MAX_SPECTRUM_GROUPS} groups")
        try:
            bounds = [float(x) for x in raw_bounds]
            values = [float(x) for x in raw_values]
        except (TypeError, ValueError, OverflowError):
            raise ValueError("the group bounds and fluxes of the spectrum must be numbers") from None
        if len(bounds) < 2 or len(values) != len(bounds) - 1:
            raise ValueError("a spectrum needs bounds_ev (n+1 group bounds in eV) and flux (n values, one per group)")
        # NaN and inf first: a comparison with NaN is always false, so NaN gets past the tests below (and a NaN
        # bound stopped the C++ code)
        if not all(math.isfinite(x) for x in bounds + values):
            raise ValueError("the group bounds and fluxes of the spectrum must be finite numbers")
        if any(not b > a for a, b in zip(bounds, bounds[1:])) or not bounds[0] > 0:
            raise ValueError("the group bounds must be more than 0 and go up")
        if any(v < 0 for v in values) or not max(values) > 0:
            raise ValueError("the group fluxes must be 0 or more, and not all 0")
        # only the shape is used, so scale the fluxes to a largest value of 1: the sum then cannot overflow or be
        # so small that the flux of a step divided by it is inf
        largest = max(values)
        values = [v / largest for v in values]
        return bounds, values, f"custom: {len(values)} groups, {bounds[0]:g} eV to {bounds[-1]:g} eV"
    raise ValueError("the spectrum is a name ('14MeV', 'thermal') or {'bounds_ev': [...], 'flux': [...]}")


def steps_from_input(irradiation_time, cooling_times, schedule, flux) -> list:
    """The steps from a schedule [{"duration": "1 d", "flux": 1e14}, ...], or from an irradiation time and
    the cooling times to report (each is the time after the end of the irradiation)."""
    steps = []
    if schedule:
        if irradiation_time is not None or cooling_times:
            raise ValueError("give a schedule, or irradiation_time and cooling_times, not both")
        for i, item in enumerate(schedule):
            if not isinstance(item, dict) or "duration" not in item:
                raise ValueError(f"schedule step {i + 1}: give {{'duration': '1 d', 'flux': 1e14}}")
            step_flux = item.get("flux")
            if step_flux is None:
                step_flux = flux if flux is not None else DEFAULT_FLUX
            try:
                step_flux = float(step_flux)
            except (TypeError, ValueError):
                raise ValueError(f"schedule step {i + 1}: the flux must be a number in n/cm2/s") from None
            if step_flux < 0 or not math.isfinite(step_flux):
                raise ValueError(f"schedule step {i + 1}: the flux must be 0 or more")
            steps.append(Step(parse_time(item["duration"], f"schedule step {i + 1}"), step_flux))
    else:
        if irradiation_time is None:
            raise ValueError("give irradiation_time (for example '1 d'), or a schedule")
        flux = DEFAULT_FLUX if flux is None else float(flux)
        if flux < 0 or not math.isfinite(flux):
            raise ValueError("the flux must be 0 or more")
        steps.append(Step(parse_time(irradiation_time, "irradiation_time"), flux))
        times = sorted(parse_time(t, "cooling_times") for t in (cooling_times or []))
        previous = 0.0
        for t in times:
            if t > previous:
                steps.append(Step(t - previous, 0.0))
                previous = t
    steps = [s for s in steps if s.duration > 0]
    if not steps:
        raise ValueError("all the steps have a duration of 0")
    if len(steps) > MAX_STEPS:
        raise ValueError(f"give at most {MAX_STEPS} steps")
    if sum(s.duration for s in steps) > MAX_TOTAL_TIME:
        raise ValueError("the total time must be at most 10000 years")
    return steps


class Point:
    """The atoms of all the nuclides at one time. Most nuclides have no atoms, so only the others are kept: a run
    can have 600 points of 4000 nuclides."""

    __slots__ = ("time", "step", "end_of_step", "size", "_index", "_values")

    def __init__(self, time: float, step: int, end_of_step: bool, atoms: np.ndarray):
        self.time = time              # s since the start
        self.step = step              # index of the step
        self.end_of_step = end_of_step
        atoms = np.asarray(atoms, dtype=float)
        self.size = len(atoms)
        self._index = np.flatnonzero(atoms).astype(np.int32)
        self._values = atoms[self._index]

    @property
    def atoms(self) -> np.ndarray:
        """The atoms of each nuclide (no sink)."""
        atoms = np.zeros(self.size)
        atoms[self._index] = self._values
        return atoms

    @property
    def nbytes(self) -> int:
        return self._index.nbytes + self._values.nbytes + 64


@dataclass
class RunResult:
    run_id: str
    owner: str
    created: float
    material: dict           # the input, for the record
    initial_atoms: dict
    steps: list
    spectrum: str
    points: list
    notes: list = field(default_factory=list)
    solve_seconds: float = 0.0

    @property
    def nbytes(self) -> int:
        """About the memory of the result."""
        return sum(p.nbytes for p in self.points)

    @property
    def irradiation_end(self) -> float:
        """The time of the end of the last irradiation step (0 if no irradiation)."""
        t, end = 0.0, 0.0
        for s in self.steps:
            t += s.duration
            if s.flux > 0:
                end = t
        return end

    def point_at(self, time_s=None, after: float = 0.0) -> Point:
        """The point at time_s after the time `after` (both in s since the start; for a cooling time, `after` is
        the end of the irradiation), or the last point. The time must be a point of the run: the end of a step,
        or one of the sub-steps. The tolerance is relative to time_s, so a short cooling time after a long
        irradiation does not give a point a few seconds away."""
        if time_s is None:
            return self.points[-1]
        if time_s < 0:
            raise ValueError("the time must be 0 or more")
        best = min(self.points, key=lambda p: abs((p.time - after) - time_s))
        # the times are sums of the step durations, so they have a rounding error of the size of the total time
        tolerance = 1e-6 * max(time_s, 1e-3) + 1e-12 * max(after, best.time)
        if abs((best.time - after) - time_s) > tolerance:
            frame = "after the end of the irradiation" if after > 0 else "since the start"
            ends = [format_time(p.time - after) for p in self.points if p.end_of_step and p.time >= after]
            raise ValueError(f"the run has no result at {format_time(time_s)} {frame}; the step ends are at "
                             + ", ".join(ends) + ". Run the irradiation again with this time in cooling_times.")
        return best


def sub_step_times(duration: float, count: int = SUB_STEPS) -> list:
    """The times in a step (from its start) of the points of the results, log spaced, the last is the duration."""
    if count <= 1:
        return [duration]
    times = np.geomspace(duration * 1e-3, duration, count)
    times[-1] = duration
    return [float(t) for t in times]


# the memory for the results of all the users together; the oldest runs go first
DEFAULT_MAX_RESULTS_BYTES = 500 * 1024 * 1024


class Runner:
    """Runs the models and keeps the results: the last max_runs of each user, and at most max_results_bytes for
    all the users together (the oldest runs are removed first, the run just made is always kept)."""

    def __init__(self, loader, max_runs: int = 20, reaction_cache_size: int = 4,
                 max_results_bytes: int = DEFAULT_MAX_RESULTS_BYTES):
        self.loader = loader
        self.max_runs = max_runs
        self.max_results_bytes = max_results_bytes
        self._runs = OrderedDict()
        self._matrices = OrderedDict()
        self._reaction_cache_size = reaction_cache_size
        # one solve at a time: the runs are short, and the nuclear data is shared
        self._solve_lock = threading.Lock()
        self._runs_lock = threading.Lock()
        # the runs in progress of each user (solving, or waiting for the solve lock), so that one user cannot
        # take all the threads of the server with runs that wait
        self._in_progress = {}
        self.max_in_progress = MAX_RUNS_IN_PROGRESS
        self.max_in_progress_per_user = MAX_RUNS_IN_PROGRESS_PER_USER
        self.solve_wait_seconds = SOLVE_WAIT_SECONDS
        self.max_solve_seconds = MAX_SOLVE_SECONDS

    # ------------------------------------------------------------------ runs

    def _reaction_matrix(self, data, bounds, values):
        key = (tuple(bounds), tuple(values))
        if key in self._matrices:
            self._matrices.move_to_end(key)
            return self._matrices[key]
        mat = data.reaction_matrix(bounds, values)
        self._matrices[key] = mat
        while len(self._matrices) > self._reaction_cache_size:
            self._matrices.popitem(last=False)
        return mat

    def _enter(self, owner: str) -> None:
        with self._runs_lock:
            if self._in_progress.get(owner, 0) >= self.max_in_progress_per_user:
                raise RuntimeError(f"you have {self.max_in_progress_per_user} runs in progress already: wait for "
                                   "them to end, then try again")
            if sum(self._in_progress.values()) >= self.max_in_progress:
                raise RuntimeError("the server is busy with the runs of other users, try again in a minute")
            self._in_progress[owner] = self._in_progress.get(owner, 0) + 1

    def _leave(self, owner: str) -> None:
        with self._runs_lock:
            count = self._in_progress.get(owner, 0) - 1
            if count > 0:
                self._in_progress[owner] = count
            else:
                self._in_progress.pop(owner, None)

    def run(self, owner: str, initial_atoms: dict, steps: list, spectrum, material: dict, notes: list) -> RunResult:
        data, tables = self.loader.get()
        self._enter(owner)
        try:
            if not self._solve_lock.acquire(timeout=self.solve_wait_seconds):
                raise TimeoutError("the server is busy with the runs of other users, try again in a minute")
            try:
                return self._run(data, tables, owner, initial_atoms, steps, spectrum, material, notes)
            finally:
                self._solve_lock.release()
        finally:
            self._leave(owner)

    def _run(self, data, tables, owner, initial_atoms, steps, spectrum, material, notes) -> RunResult:
        """The solve, with the solve lock held."""
        nis = data.nis
        bounds, values, spectrum_text = spectrum
        start = time.time()
        reaction = None
        if any(s.flux > 0 for s in steps):
            reaction = self._reaction_matrix(data, bounds, values)
            total = float(sum(values))
        n = np.asarray(nis.core.vector_from_inventory(initial_atoms), dtype=float)
        size = len(tables.names)
        points = [Point(0.0, -1, True, n[:size].copy())]
        t = 0.0
        matrices = {}
        for index, step in enumerate(steps):
            if step.flux not in matrices:
                # the spectrum is in n/cm2/s, and its total is the flux of the step
                matrices[step.flux] = (data.decay_mat if step.flux == 0
                                       else data.decay_mat + reaction * (step.flux / total))
            mat = matrices[step.flux]
            previous = 0.0
            times = sub_step_times(step.duration)
            for k, sub in enumerate(times):
                if time.time() - start > self.max_solve_seconds:
                    raise RuntimeError(f"the run took more than {self.max_solve_seconds:g} s, so it was "
                                       "stopped: give fewer steps or a smaller material")
                n = np.asarray(nis.core.solve_reachable(mat, n, sub - previous, nis.core.solver_cram48))
                if not np.all(np.isfinite(n)):
                    raise RuntimeError("the solve gave numbers that are not finite; check the input")
                # a solve can give very small negative values, they are 0
                n = np.where(n < 0, 0.0, n)
                previous = sub
                points.append(Point(t + sub, index, k == len(times) - 1, n[:size].copy()))
            t += step.duration
        result = RunResult(run_id=uuid.uuid4().hex[:8], owner=owner, created=time.time(), material=material,
                           initial_atoms=initial_atoms, steps=steps, spectrum=spectrum_text, points=points,
                           notes=notes, solve_seconds=time.time() - start)
        with self._runs_lock:
            self._runs[result.run_id] = result
            mine = [k for k, r in self._runs.items() if r.owner == owner]
            for k in mine[:-self.max_runs]:
                del self._runs[k]
            total = sum(r.nbytes for r in self._runs.values())
            while total > self.max_results_bytes and len(self._runs) > 1:
                _, removed = self._runs.popitem(last=False)
                total -= removed.nbytes
        return result

    def get(self, owner: str, run_id: str) -> RunResult:
        with self._runs_lock:
            result = self._runs.get(run_id)
        if result is None or result.owner != owner:
            raise ValueError(f"no run {run_id!r}; the runs are kept in memory only, so run the irradiation again")
        return result

    def list(self, owner: str) -> list:
        with self._runs_lock:
            return [r for r in self._runs.values() if r.owner == owner]


# ---------------------------------------------------------------------- the quantities of a result

QUANTITIES = {
    "activity": "Bq",
    "decay_heat": "W",
    "atoms": "atoms",
    "mass": "g",
}


def nuclide_values(tables, atoms: np.ndarray, quantity: str) -> np.ndarray:
    """The quantity for each nuclide, from its atoms."""
    lam = tables.decay_constants[:len(atoms)]
    if quantity == "activity":
        return lam * atoms
    if quantity == "decay_heat":
        return lam * atoms * tables.energy["total"][:len(atoms)]
    if quantity == "atoms":
        return atoms
    if quantity == "mass":
        return atoms * tables.mass_numbers / 6.02214076e23
    raise ValueError(f"unknown quantity {quantity!r}, use one of {sorted(QUANTITIES)}")


def totals(tables, atoms: np.ndarray) -> dict:
    activities = tables.decay_constants[:len(atoms)] * atoms
    heat = {part: float(activities @ e[:len(atoms)]) for part, e in tables.energy.items()}
    return {"activity_Bq": sig(float(activities.sum())), "decay_heat_W": sig(heat["total"]),
            "beta_heat_W": sig(heat["light_particle"]), "gamma_heat_W": sig(heat["electromagnetic"]),
            "alpha_heat_W": sig(heat["heavy_particle"])}


def half_life_text(lam: float) -> str:
    if lam <= 0:
        return "stable"
    return format_time(math.log(2) / lam)


def top_nuclides(tables, atoms: np.ndarray, quantity: str, count: int, exclude=()) -> list:
    """The names of the count nuclides with the largest quantity, largest first."""
    values = nuclide_values(tables, atoms, quantity)
    order = np.argsort(values)[::-1]
    names = []
    for i in order:
        if values[i] <= 0 or len(names) >= count:
            break
        if tables.names[i] not in exclude:
            names.append(tables.names[i])
    return names


def nuclide_table(tables, atoms: np.ndarray, quantity: str, names: list) -> list:
    """Rows of the nuclides: the quantity, its fraction of the total, the half-life."""
    values = nuclide_values(tables, atoms, quantity)
    total = float(values.sum())
    unit = QUANTITIES[quantity]
    rows = []
    for name in names:
        i = tables.index[name]
        value = float(values[i])
        row = {"nuclide": name, f"{quantity}_{unit}": sig(value),
               "half_life": half_life_text(float(tables.decay_constants[i]))}
        if total > 0:
            row["percent_of_total"] = sig(100.0 * value / total, 3)
        rows.append(row)
    return rows
