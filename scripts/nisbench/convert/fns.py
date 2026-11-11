"""
FNS decay heat benchmark: from the FISPACT-II validation inputs to the standard format.

The FNS (JAEA) experiments irradiate thin foils with 14 MeV neutrons, then measure the decay heat.
The original files are the FISPACT-II validation inputs for these experiments:

    <source>/<material>/
        TENDL-2017_<case>.i          FISPACT-II input (mass, flux amplitude, times)
        TENDL-2017_<case>.out        FISPACT-II output (the initial atoms)
        TENDL-2017_<case>.nuclides   FISPACT-II decay heat at each step (kW/kg)
        <case>_fluxes                709 group flux (high to low energy)
        <case>.exp                   measured decay heat: time, uW/g, uncertainty

Each case is the irradiation step, then one cooling step to each measured time. The FISPACT-II
heat is interpolated to the measured times (log of the heat linear in time).
"""

import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

from ..format import Case, Measurement

NAME = "fns"
METADATA = {
    "name": NAME,
    "title": "FNS decay heat",
    "category": "fusion",
    "description": "Decay heat of thin foils after a 5 minute or 7 hour irradiation with 14 MeV (D-T) neutrons "
                   "at the Fusion Neutronics Source (FNS) of JAEA, measured from about 1 minute to about 400 days "
                   "after shutdown. 132 cases of elements and alloys.",
    "references": [
        "F. Maekawa, M. Wada, Y. Ikeda, Decay heat experiment and validation of calculation code systems "
        "for fusion reactor, JAERI 99-055 (1999)",
        "J.-Ch. Sublet et al., FISPACT-II: an advanced simulation system for activation, transmutation and "
        "material modelling, Nuclear Data Sheets 139 (2017) 77-137",
    ],
    "source": "The FISPACT-II validation inputs for FNS (UKAEA), converted with scripts/bench.py convert fns. "
              "The spectrum is the 709 group flux of the FISPACT-II input. FISPACT-II (TENDL-2017) is the "
              "FISPACT-II output, interpolated to the measured times.",
}
SOURCE_HELP = "the folder of the FISPACT-II FNS validation inputs (<material>/TENDL-2017_<case>.i, ...)"

GROUP_STRUCTURE = "CCFE-709"
CALCULATION = "FISPACT-II (TENDL-2017)"

_TIME_UNITS = {"SECS": 1.0, "MINS": 60.0, "HOURS": 3600.0, "DAYS": 86400.0, "YEARS": 365.25 * 86400.0}
# 1 kW/kg = 1e6 uW/g
_KW_PER_KG_TO_UW_PER_G = 1e6
# e.g. "  Co 58m      1.40935E+10 ..." or "  Co 59    #> 1.02186E+22 ..."
_NUCLIDE_LINE = re.compile(r"^\s+([A-Z][a-z]?)\s*(\d+)([a-z]?)\s+[#>]*\s*([0-9.]+E[+-]\d+)")


@dataclass
class _Original:
    mass_kg: float = 0.0
    # n/cm2/s
    flux_amplitude: float = 0.0
    irradiation_time: float = 0.0
    cooling_steps: list = field(default_factory=list)
    initial_atoms: dict = field(default_factory=dict)
    # uW/g at the end of each cooling step
    fispact_heat: list = field(default_factory=list)
    exp_time: list = field(default_factory=list)
    exp_heat: list = field(default_factory=list)
    exp_uncertainty: list = field(default_factory=list)
    # low to high energy
    flux_shape: list = field(default_factory=list)


def _to_seconds(tokens: list) -> float:
    unit = tokens[2].upper() if len(tokens) > 2 else "SECS"
    return float(tokens[1]) * _TIME_UNITS.get(unit, 1.0)


def _parse_input(o: _Original, path: Path) -> None:
    """The simple form of all the FNS inputs: one irradiation, then cooling."""
    cooling = False
    for line in path.read_text().splitlines():
        tokens = line.split()
        if not tokens:
            continue
        key = tokens[0].upper()
        if key == "MASS":
            o.mass_kg = float(tokens[1])
        elif key == "FLUX" and float(tokens[1]) > 0:
            o.flux_amplitude = float(tokens[1])
        elif key == "ZERO":
            cooling = True
        elif key == "TIME":
            if cooling:
                o.cooling_steps.append(_to_seconds(tokens))
            else:
                o.irradiation_time += _to_seconds(tokens)


def _parse_initial_atoms(o: _Original, path: Path) -> None:
    """The first nuclide inventory in the FISPACT-II output."""
    current = None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("  NUCLIDE        ATOMS"):
            current = {}
        elif current is not None:
            if line.startswith("0  TOTAL NUMBER OF NUCLIDES"):
                o.initial_atoms = current
                return
            match = _NUCLIDE_LINE.match(line)
            if match:
                element, mass, state, atoms = match.groups()
                current[f"{element}{int(mass)}{state}"] = float(atoms)
    raise ValueError(f"no inventory in {path}")


def _parse_heat(o: _Original, path: Path) -> None:
    """The first table in the .nuclides file is the heat output, step 1 is shutdown."""
    total = None
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            if "Heat Output" in line:
                total = []
            elif total:
                break
            continue
        if total is not None and line.strip():
            total.append(float(line.split()[3]) * _KW_PER_KG_TO_UW_PER_G)
    o.fispact_heat = total[1:]


def _cumulative(steps: list) -> list:
    times, total = [], 0.0
    for step in steps:
        total += step
        times.append(total)
    return times


def _parse_experiment(o: _Original, path: Path) -> None:
    raw_times = []
    for line in path.read_text().splitlines():
        tokens = line.split()
        # some files end with rows of zeros as padding
        if len(tokens) >= 3 and float(tokens[0]) > 0 and float(tokens[1]) > 0:
            raw_times.append(float(tokens[0]))
            o.exp_heat.append(float(tokens[1]))
            o.exp_uncertainty.append(float(tokens[2]))

    # the measured times are in minutes, hours or days (it depends on the case):
    # take the unit for which they are closest to the FISPACT-II cooling steps
    step_times = _cumulative(o.cooling_steps)

    def mismatch(unit):
        return sum(min(abs(math.log(t * unit / s)) for s in step_times) for t in raw_times if t > 0)
    unit = min((_TIME_UNITS["MINS"], _TIME_UNITS["HOURS"], _TIME_UNITS["DAYS"]), key=mismatch)
    o.exp_time = [t * unit for t in raw_times]


def _parse_fluxes(o: _Original, path: Path) -> None:
    # FISPACT-II gives high to low energy, NIS needs low to high
    o.flux_shape = list(reversed([float(v) for v in path.read_text().split()[:709]]))


def interpolate_log(times: list, values: list, t: float) -> float:
    """Log of the value is linear in time, NaN outside the range."""
    for i in range(len(times) - 1):
        t0, t1 = times[i], times[i + 1]
        if t0 <= t <= t1 and values[i] > 0 and values[i + 1] > 0:
            f = (t - t0) / (t1 - t0)
            return math.exp((1 - f) * math.log(values[i]) + f * math.log(values[i + 1]))
    return math.nan


def _title(material: str, name: str) -> str:
    match = re.match(r"(\d{4})exp_(\d+)(min|hour)$", name)
    if not match:
        return f"{material} {name}"
    year, length, unit = match.groups()
    return f"{material}, {length} {'minute' if unit == 'min' else 'hour'} irradiation, {year} experiment"


def _case(material: str, name: str, o: _Original) -> Case:
    cooling = [t1 - t0 for t0, t1 in zip([0.0] + o.exp_time[:-1], o.exp_time)]
    # FISPACT-II gives the flux in n/cm2/s, the format in n/m2/s
    steps = [(o.irradiation_time, o.flux_amplitude * 1e4)] + [(dt, 0.0) for dt in cooling]
    step_times = _cumulative(o.cooling_steps)
    fispact = [interpolate_log(step_times, o.fispact_heat, t) for t in o.exp_time]
    measurement = Measurement(
        quantity="decay_heat", unit="uW/g", steps=list(range(1, len(steps))),
        values=o.exp_heat, uncertainties=o.exp_uncertainty,
        calculations={CALCULATION: [None if math.isnan(v) else v for v in fispact]})
    return Case(benchmark=NAME, name=f"{material}/{name}", title=_title(material, name), mass_kg=o.mass_kg,
                atoms=o.initial_atoms, steps=steps, measurements=[measurement],
                spectrum=(None, o.flux_shape), group_structure=GROUP_STRUCTURE)


def convert(source: Path, only: list = None) -> list:
    """only: material names, all if empty."""
    cases = []
    for input_file in sorted(Path(source).glob("*/TENDL-2017_*.i")):
        directory = input_file.parent
        material, name = directory.name, input_file.stem.replace("TENDL-2017_", "")
        if only and material not in only:
            continue
        o = _Original()
        try:
            _parse_input(o, input_file)
            _parse_initial_atoms(o, input_file.with_suffix(".out"))
            _parse_heat(o, input_file.with_suffix(".nuclides"))
            _parse_experiment(o, directory / f"{name}.exp")
            _parse_fluxes(o, directory / f"{name}_fluxes")
        except (OSError, IndexError, ValueError) as e:
            print(f"skip {material} {name}: {e}", file=sys.stderr)
            continue
        if not o.exp_time or sorted(o.exp_time) != o.exp_time or len(set(o.exp_time)) != len(o.exp_time):
            print(f"skip {material} {name}: no measurements in time order", file=sys.stderr)
            continue
        cases.append(_case(material, name, o))
    return cases
