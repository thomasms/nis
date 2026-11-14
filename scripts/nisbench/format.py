"""
Read, check and write the standard benchmark format (version 1). The format is in benchmarks/README.md.
"""

import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

FORMAT_VERSION = 1
BENCHMARK_FORMAT = "nis-benchmark"
CASE_FORMAT = "nis-benchmark-case"

BENCHMARKS_DIR = Path(__file__).resolve().parents[2] / "benchmarks"
GROUP_STRUCTURES_DIR = BENCHMARKS_DIR / "group_structures"

CATEGORIES = ("fusion", "fission", "other")

# quantity -> the part of the NIS result it is
QUANTITIES = {
    "decay_heat": "total decay heat",
    "beta_heat": "decay heat of the light particles (beta and electrons)",
    "gamma_heat": "decay heat of the electromagnetic radiation (gamma and X rays)",
    "alpha_heat": "decay heat of the heavy particles (alpha)",
    "activity": "total activity",
}

_MEV_TO_J = 1.602176634e-13

# unit -> (factor from the SI value: W for heat, Bq for activity, normalisation, multiply by the cooling time)
# normalisation: None (the whole sample), "mass" (divide by sample.mass_kg), "fission" (divide by "fissions")
# or "fission_rate" (divide by the fission rate of the last step with a fission source)
HEAT_UNITS = {
    "W": (1.0, None, False),
    "MeV/s": (1.0 / _MEV_TO_J, None, False),
    "W/kg": (1.0, "mass", False),
    "kW/kg": (1e-3, "mass", False),
    "mW/g": (1.0, "mass", False),
    "uW/g": (1e3, "mass", False),
    "W/fission": (1.0, "fission", False),
    "MeV/s/fission": (1.0 / _MEV_TO_J, "fission", False),
    # t * P(t), the usual way to give the decay heat after a fission pulse
    "MeV/fission": (1.0 / _MEV_TO_J, "fission", True),
    # the decay heat over the fission rate of the irradiation, the usual way after an irradiation of some length
    "MeV/(fission/s)": (1.0 / _MEV_TO_J, "fission_rate", False),
}
ACTIVITY_UNITS = {
    "Bq": (1.0, None, False),
    "Bq/kg": (1.0, "mass", False),
    "Bq/g": (1e-3, "mass", False),
    "Bq/fission": (1.0, "fission", False),
}


def units_for(quantity: str) -> dict:
    return ACTIVITY_UNITS if quantity == "activity" else HEAT_UNITS


_NUCLIDE = re.compile(r"^[A-Z][a-z]?\d{1,3}[mn]?$")
_NAME = re.compile(r"^[A-Za-z0-9_.+-]+(/[A-Za-z0-9_.+-]+)*$")


class FormatError(ValueError):
    pass


@dataclass
class Measurement:
    quantity: str
    unit: str
    # the step (0 is the first) at the end of which each value is measured
    steps: list
    values: list
    uncertainties: list
    # calculation name -> values, None where there is no value
    calculations: dict = field(default_factory=dict)


@dataclass
class Case:
    benchmark: str
    name: str
    title: str
    mass_kg: float
    atoms: dict
    # (dt in s, flux in n/m2/s), flux 0 for cooling
    steps: list
    measurements: list
    # (energy bounds in eV, ascending, values), or None for decay only
    spectrum: tuple = None
    group_structure: str = None
    fissions: float = None
    path: Path = None
    # a fission source in place of a neutron flux: {"nuclide", "neutron_energy_ev", "fissions_at_start"}.
    # Each fission adds the independent fission yields of the nuclide at that neutron energy.
    fission_source: dict = None
    # the fission rate (fissions/s) of each step, with a fission source
    fission_rates: list = None

    @property
    def deck_path(self) -> Path:
        return self.path.with_suffix(".inp")

    @property
    def has_deck(self) -> bool:
        """A NIS input deck cannot give a fission source, so such a case has no deck (pynis runs it)."""
        return self.fission_source is None

    def step_fission_rates(self) -> list:
        return self.fission_rates or [0.0] * len(self.steps)

    @property
    def last_fission_rate(self) -> float:
        rates = [r for r in self.step_fission_rates() if r > 0]
        return rates[-1] if rates else None

    def end_times(self) -> list:
        """The time at the end of each step, from the start, in s."""
        times, t = [], 0.0
        for dt, _ in self.steps:
            t += dt
            times.append(t)
        return times

    def cooling_times(self, steps: list) -> list:
        """The time of the end of each given step after the end of the last irradiation step, in s.
        With no irradiation, the time from the start."""
        end = self.end_times()
        irradiated = [i for i, ((_, flux), rate) in enumerate(zip(self.steps, self.step_fission_rates()))
                      if flux > 0 or rate > 0]
        shutdown = end[irradiated[-1]] if irradiated else 0.0
        return [end[s] - shutdown for s in steps]


@dataclass
class Benchmark:
    name: str
    title: str
    category: str
    description: str
    references: list
    source: str
    path: Path
    cases: list = field(default_factory=list)


def read_group_structure(name: str) -> list:
    """The energy bounds (eV, ascending) of a group structure in benchmarks/group_structures."""
    path = GROUP_STRUCTURES_DIR / f"{name}.txt"
    if not path.exists():
        known = sorted(p.stem for p in GROUP_STRUCTURES_DIR.glob("*.txt"))
        raise FormatError(f"unknown group structure {name!r}, known: {', '.join(known)}")
    values = []
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            values += [float(v) for v in line.split()]
    return values


def _require(condition: bool, where: str, message: str) -> None:
    if not condition:
        raise FormatError(f"{where}: {message}")


def _number(value, where: str, what: str, minimum: float = None, positive: bool = False) -> float:
    _require(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value),
             where, f"{what} must be a number")
    if positive:
        _require(value > 0, where, f"{what} must be more than 0")
    if minimum is not None:
        _require(value >= minimum, where, f"{what} must be {minimum} or more")
    return float(value)


def _header(data: dict, where: str, expected_format: str) -> None:
    _require(isinstance(data, dict), where, "must be a JSON object")
    _require(data.get("format") == expected_format, where, f'"format" must be "{expected_format}"')
    _require(data.get("version") == FORMAT_VERSION, where, f'"version" must be {FORMAT_VERSION}')


def _spectrum(data, where: str):
    if data is None:
        return None, None
    _require(isinstance(data, dict), where, '"spectrum" must be an object or null')
    values = [_number(v, where, "a spectrum value", minimum=0.0) for v in data.get("values", [])]
    name = data.get("group_structure")
    if name is not None:
        _require("energy_bounds_ev" not in data, where, 'give "group_structure" or "energy_bounds_ev", not both')
        bounds = read_group_structure(name)
    else:
        _require("energy_bounds_ev" in data, where, 'the spectrum needs "group_structure" or "energy_bounds_ev"')
        bounds = [_number(v, where, "an energy bound", minimum=0.0) for v in data["energy_bounds_ev"]]
    _require(len(bounds) == len(values) + 1, where,
             f"the spectrum has {len(values)} values, so it needs {len(values) + 1} energy bounds, not {len(bounds)}")
    _require(all(a < b for a, b in zip(bounds, bounds[1:])), where, "the energy bounds must go up")
    _require(sum(values) > 0, where, "the spectrum must not be all 0")
    return (bounds, values), name


def _measurement(data, where: str, nsteps: int) -> Measurement:
    _require(isinstance(data, dict), where, "must be an object")
    quantity = data.get("quantity")
    _require(quantity in QUANTITIES, where, f"unknown quantity {quantity!r}, known: {', '.join(QUANTITIES)}")
    unit = data.get("unit")
    units = units_for(quantity)
    _require(unit in units, where, f"unknown unit {unit!r} for {quantity}, known: {', '.join(units)}")

    steps = data.get("steps")
    _require(isinstance(steps, list) and steps, where, '"steps" must be a list that is not empty')
    _require(all(isinstance(s, int) and 0 <= s < nsteps for s in steps), where,
             f'each of "steps" must be a step index from 0 to {nsteps - 1}')
    _require(all(a < b for a, b in zip(steps, steps[1:])), where, '"steps" must go up')

    values = [_number(v, where, "a value") for v in data.get("values", [])]
    uncertainties = [_number(v, where, "an uncertainty", minimum=0.0) for v in data.get("uncertainties", [])]
    _require(len(values) == len(steps), where, f'"values" has {len(values)} items, "steps" has {len(steps)}')
    _require(len(uncertainties) == len(steps), where,
             f'"uncertainties" has {len(uncertainties)} items, "steps" has {len(steps)}')

    calculations = data.get("calculations", {})
    _require(isinstance(calculations, dict), where, '"calculations" must be an object')
    checked = {}
    for name, calc in calculations.items():
        _require(isinstance(calc, list) and len(calc) == len(steps), where,
                 f'calculation {name!r} must be a list with one value (or null) for each step')
        checked[name] = [None if v is None else _number(v, where, f"a value of {name!r}") for v in calc]
    return Measurement(quantity, unit, list(steps), values, uncertainties, checked)


def parse_case(data: dict, where: str = "case") -> Case:
    _header(data, where, CASE_FORMAT)
    for key in ("benchmark", "name"):
        _require(isinstance(data.get(key), str) and _NAME.match(data[key]), where,
                 f'"{key}" must be a name of letters, digits, "_.+-" (and "/" in "name")')

    sample = data.get("sample")
    _require(isinstance(sample, dict), where, '"sample" must be an object')
    mass = sample.get("mass_kg")
    mass = None if mass is None else _number(mass, where, "sample.mass_kg", positive=True)
    atoms = sample.get("atoms")
    source = data.get("fission_source")
    _require(isinstance(atoms, dict) and (atoms or source is not None), where,
             '"sample.atoms" must be an object that is not empty (it can be empty with a "fission_source")')
    for nuclide, n in atoms.items():
        _require(bool(_NUCLIDE.match(nuclide)), where, f"{nuclide!r} is not a nuclide name, such as Co59 or Ag110m")
        _number(n, where, f"the atoms of {nuclide}", minimum=0.0)

    raw_steps = data.get("steps")
    _require(isinstance(raw_steps, dict), where, '"steps" must be an object with "dt_s" and "flux_n_m2_s"')
    dts = raw_steps.get("dt_s")
    fluxes = raw_steps.get("flux_n_m2_s")
    _require(isinstance(dts, list) and dts, where, '"steps.dt_s" must be a list that is not empty')
    _require(isinstance(fluxes, list) and len(fluxes) == len(dts), where,
             '"steps.flux_n_m2_s" must be a list with one flux for each step')
    steps = [(_number(dt, where, "dt_s", positive=True), _number(flux, where, "flux_n_m2_s", minimum=0.0))
             for dt, flux in zip(dts, fluxes)]

    rates = raw_steps.get("fission_rate_per_s")
    if rates is not None:
        _require(source is not None, where, '"steps.fission_rate_per_s" needs a "fission_source"')
        _require(isinstance(rates, list) and len(rates) == len(dts), where,
                 '"steps.fission_rate_per_s" must be a list with one rate for each step')
        rates = [_number(r, where, "fission_rate_per_s", minimum=0.0) for r in rates]
    if source is not None:
        _require(isinstance(source, dict), where, '"fission_source" must be an object')
        _require(isinstance(source.get("nuclide"), str) and bool(_NUCLIDE.match(source["nuclide"])), where,
                 '"fission_source.nuclide" must be a nuclide name, such as U235')
        source = {"nuclide": source["nuclide"],
                  "neutron_energy_ev": _number(source.get("neutron_energy_ev"), where,
                                               "fission_source.neutron_energy_ev", minimum=0.0),
                  "fissions_at_start": _number(source.get("fissions_at_start", 0.0), where,
                                               "fission_source.fissions_at_start", minimum=0.0)}
        _require(source["fissions_at_start"] > 0 or any(r > 0 for r in rates or []), where,
                 'the fission source needs "fissions_at_start" or a step with a fission rate')

    spectrum, group_structure = _spectrum(data.get("spectrum"), where)
    _require(spectrum is not None or all(flux == 0 for _, flux in steps), where,
             "a step has a flux, so the case needs a spectrum")

    fissions = data.get("fissions")
    fissions = None if fissions is None else _number(fissions, where, "fissions", positive=True)
    if source is not None:
        made = source["fissions_at_start"] + sum(dt * r for (dt, _), r in zip(steps, rates or [0.0] * len(steps)))
        _require(fissions is None, where, 'with a "fission_source", "fissions" is the fissions of the source')
        fissions = made

    raw = data.get("measurements")
    _require(isinstance(raw, list) and raw, where, '"measurements" must be a list that is not empty')
    measurements = [_measurement(m, f"{where}, measurement {i}", len(steps)) for i, m in enumerate(raw)]
    for m in measurements:
        normalisation = units_for(m.quantity)[m.unit][1]
        _require(normalisation != "mass" or mass is not None, where, f"the unit {m.unit} needs sample.mass_kg")
        _require(normalisation != "fission" or fissions is not None, where, f'the unit {m.unit} needs "fissions"')
        _require(normalisation != "fission_rate" or any(r > 0 for r in rates or []), where,
                 f"the unit {m.unit} needs a step with a fission rate")

    return Case(benchmark=data["benchmark"], name=data["name"], title=data.get("title", data["name"]),
                mass_kg=mass, atoms={k: float(v) for k, v in atoms.items()}, steps=steps,
                measurements=measurements, spectrum=spectrum, group_structure=group_structure, fissions=fissions,
                fission_source=source, fission_rates=rates)


def load_case(path: Path) -> Case:
    path = Path(path)
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError as e:
        raise FormatError(f"{path}: not valid JSON: {e}") from e
    case = parse_case(data, str(path))
    case.path = path
    return case


def parse_benchmark(data: dict, where: str = "benchmark") -> dict:
    _header(data, where, BENCHMARK_FORMAT)
    _require(isinstance(data.get("name"), str) and _NAME.match(data["name"]) and "/" not in data["name"], where,
             '"name" must be a name of letters, digits and "_.+-"')
    _require(data.get("category") in CATEGORIES, where, f'"category" must be one of {", ".join(CATEGORIES)}')
    for key in ("title", "description", "source"):
        _require(isinstance(data.get(key), str) and data[key], where, f'"{key}" must be a text that is not empty')
    _require(isinstance(data.get("references", []), list), where, '"references" must be a list')
    return data


def load_benchmark(path: Path, cases: list = None) -> Benchmark:
    """path is the benchmark folder (with benchmark.json). cases: only these case names (or a prefix of
    them, such as "Co" for "Co/2000exp_5min"), all if empty."""
    path = Path(path)
    meta_path = path / "benchmark.json"
    if not meta_path.exists():
        raise FormatError(f"{path}: no benchmark.json")
    data = parse_benchmark(json.loads(meta_path.read_text()), str(meta_path))
    _require(data["name"] == path.name, str(meta_path), f'"name" must be the folder name, {path.name!r}')

    benchmark = Benchmark(name=data["name"], title=data["title"], category=data["category"],
                          description=data["description"], references=data.get("references", []),
                          source=data["source"], path=path)
    for case_path in sorted(path.rglob("*.json")):
        if case_path == meta_path:
            continue
        case = load_case(case_path)
        expected = case_path.relative_to(path).with_suffix("").as_posix()
        _require(case.name == expected, str(case_path), f'"name" must be the path in the benchmark, {expected!r}')
        _require(case.benchmark == benchmark.name, str(case_path), f'"benchmark" must be {benchmark.name!r}')
        if not cases or any(case.name == c or case.name.startswith(c.rstrip("/") + "/") for c in cases):
            benchmark.cases.append(case)
    return benchmark


def load_benchmarks(names: list = None, root: Path = BENCHMARKS_DIR, cases: list = None) -> list:
    """All the benchmarks in root (a folder with a benchmark.json), or only the given names."""
    available = sorted(p.parent.name for p in Path(root).glob("*/benchmark.json"))
    names = names or available
    unknown = [n for n in names if n not in available]
    if unknown:
        raise FormatError(f"unknown benchmark {', '.join(unknown)}, known: {', '.join(available)}")
    return [load_benchmark(Path(root) / n, cases) for n in names]


def _compact(value: float):
    # 10 significant digits, so the files do not change from one write to the next
    return float(f"{value:.10g}")


def case_to_json(case: Case) -> dict:
    """The JSON object of a case (the inverse of parse_case)."""
    data = {
        "format": CASE_FORMAT,
        "version": FORMAT_VERSION,
        "benchmark": case.benchmark,
        "name": case.name,
        "title": case.title,
        "sample": {"atoms": {k: _compact(v) for k, v in sorted(case.atoms.items())}},
        "steps": {"dt_s": [_compact(dt) for dt, _ in case.steps],
                  "flux_n_m2_s": [_compact(flux) for _, flux in case.steps]},
    }
    if case.fission_rates is not None:
        data["steps"]["fission_rate_per_s"] = [_compact(r) for r in case.fission_rates]
    if case.mass_kg is not None:
        data["sample"] = {"mass_kg": _compact(case.mass_kg), **data["sample"]}
    if case.spectrum is None:
        data["spectrum"] = None
    elif case.group_structure:
        data["spectrum"] = {"group_structure": case.group_structure, "values": [_compact(v) for v in case.spectrum[1]]}
    else:
        data["spectrum"] = {"energy_bounds_ev": [_compact(v) for v in case.spectrum[0]],
                            "values": [_compact(v) for v in case.spectrum[1]]}
    if case.fission_source is not None:
        data["fission_source"] = {k: _compact(v) if isinstance(v, float) else v for k, v in case.fission_source.items()}
    elif case.fissions is not None:
        data["fissions"] = _compact(case.fissions)
    data["measurements"] = [{
        "quantity": m.quantity,
        "unit": m.unit,
        "steps": m.steps,
        "values": [_compact(v) for v in m.values],
        "uncertainties": [_compact(v) for v in m.uncertainties],
        "calculations": {k: [None if v is None else _compact(v) for v in vs] for k, vs in m.calculations.items()},
    } for m in case.measurements]
    return data


def _dump(data: dict) -> str:
    # one line for each list of numbers, so the files are short and easy to diff
    text = json.dumps(data, indent=2)
    # 1.116e+14, not 111600000000000.0
    text = re.sub(r"(?<![\w.])(\d{7,}\.0)(?![\w.])", lambda m: f"{float(m.group(1)):.10g}", text)
    return re.sub(r"\[\s+([^\[\]{}]*?)\s+\]", lambda m: "[" + re.sub(r",\s+", ", ", m.group(1)) + "]", text) + "\n"


def write_case(case: Case, path: Path) -> None:
    data = case_to_json(case)
    # check what we write, so a converter can never make a file that does not load
    parse_case(json.loads(json.dumps(data)), str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_dump(data))
    case.path = path


def write_benchmark(data: dict, folder: Path) -> None:
    data = {"format": BENCHMARK_FORMAT, "version": FORMAT_VERSION, **data}
    parse_benchmark(data, str(folder))
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "benchmark.json").write_text(json.dumps(data, indent=2) + "\n")
