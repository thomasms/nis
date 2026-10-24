"""
Run the FNS decay heat benchmarks with NIS and compare with the experiment and FISPACT-II.

The FNS (JAEA) experiments irradiate thin foils with 14 MeV neutrons, then measure the
decay heat. This script uses the FISPACT-II validation inputs for these experiments:

    <fns-dir>/<material>/
        TENDL-2017_<case>.i          FISPACT-II input (flux amplitude, times)
        TENDL-2017_<case>.out        FISPACT-II output (initial atoms, reference)
        TENDL-2017_<case>.nuclides   FISPACT-II decay heat per step (kW/kg)
        <case>_fluxes                709 group flux (high to low energy)
        <case>.exp                   measured decay heat: time, uW/g, uncertainty

The measured times are in minutes or days (it depends on the case), we take the unit
that matches the FISPACT-II cooling steps best. NIS gives the decay heat at the measured
times, FISPACT-II is interpolated to them from its cooling steps.

Usage:
    export PYTHONPATH=./build/bin/python3
    python scripts/fns_benchmark.py --fns-dir path/to/fns \\
        --decay-data data/jeff311.dec --reaction-data data/jeff4.reac \\
        --output fns_results.csv [--materials Co Fe ...] [--jobs 4]
"""

import argparse
import csv
import math
import multiprocessing
import os
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

# pynis is imported where it is used, so the parsers work without the python bindings

_BOUNDS_FILE = Path(__file__).parent / "data" / "ccfe709_bounds.txt"

_TIME_UNITS = {
    "SECS": 1.0,
    "MINS": 60.0,
    "HOURS": 3600.0,
    "DAYS": 86400.0,
    "YEARS": 365.25 * 86400.0,
}

# 1 W/kg = 1e3 uW/g
_W_PER_KG_TO_UW_PER_G = 1e3
# 1 kW/kg = 1e6 uW/g
_KW_PER_KG_TO_UW_PER_G = 1e6

# e.g. "  Co 58m      1.40935E+10 ..." or "  Co 59    #> 1.02186E+22 ..."
_NUCLIDE_LINE = re.compile(r"^\s+([A-Z][a-z]?)\s*(\d+)([a-z]?)\s+[#>]*\s*([0-9.]+E[+-]\d+)")


@dataclass
class Case:
    material: str
    name: str
    mass_kg: float = 0.0
    flux_amplitude: float = 0.0
    irradiation_time: float = 0.0
    # each step length in seconds
    cooling_steps: list = field(default_factory=list)
    initial_atoms: dict = field(default_factory=dict)
    # after irradiation, from FISPACT-II
    fispact_atoms: dict = field(default_factory=dict)
    # uW/g at each cooling step
    fispact_heat: list = field(default_factory=list)
    # seconds after shutdown
    exp_time: list = field(default_factory=list)
    exp_heat: list = field(default_factory=list)
    exp_uncertainty: list = field(default_factory=list)
    flux_shape: list = field(default_factory=list)


def read_bounds() -> list:
    values = []
    for line in _BOUNDS_FILE.read_text().splitlines():
        if line.startswith("#"):
            continue
        values.extend(float(v) for v in line.split())
    assert len(values) == 710
    return values


def _to_seconds(tokens: list) -> float:
    value = float(tokens[1])
    unit = tokens[2].upper() if len(tokens) > 2 else "SECS"
    return value * _TIME_UNITS.get(unit, 1.0)


def parse_input(case: Case, path: Path) -> None:
    """Supports the simple form used by all FNS inputs: one irradiation, then cooling."""
    cooling = False
    for line in path.read_text().splitlines():
        tokens = line.split()
        if not tokens:
            continue
        key = tokens[0].upper()
        if key == "MASS":
            case.mass_kg = float(tokens[1])
        elif key == "FLUX":
            flux = float(tokens[1])
            if flux > 0:
                case.flux_amplitude = flux
        elif key == "ZERO":
            cooling = True
        elif key == "TIME":
            if cooling:
                case.cooling_steps.append(_to_seconds(tokens))
            else:
                case.irradiation_time += _to_seconds(tokens)


def _read_inventory_blocks(path: Path) -> list:
    """All the nuclide inventories in a FISPACT-II output, in order."""
    blocks = []
    current = None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("  NUCLIDE        ATOMS"):
            current = {}
            continue
        if current is not None:
            if line.startswith("0  TOTAL NUMBER OF NUCLIDES"):
                blocks.append(current)
                current = None
                continue
            match = _NUCLIDE_LINE.match(line)
            if match:
                element, mass, state, atoms = match.groups()
                current[f"{element}{int(mass)}{state}"] = float(atoms)
    return blocks


def parse_output(case: Case, path: Path) -> None:
    blocks = _read_inventory_blocks(path)
    case.initial_atoms = blocks[0]
    case.fispact_atoms = blocks[1]


def parse_heat(case: Case, path: Path) -> None:
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
    # drop step 1 (shutdown), the rest are the cooling steps
    case.fispact_heat = total[1:]


def cumulative(steps: list) -> list:
    times = []
    total = 0.0
    for step in steps:
        total += step
        times.append(total)
    return times


def _time_unit(raw_times: list, step_times: list) -> float:
    """The unit (in s) for which the measured times are closest to the cooling steps."""
    def mismatch(unit):
        return sum(min(abs(math.log(t * unit / s)) for s in step_times) for t in raw_times if t > 0)
    return min((_TIME_UNITS["MINS"], _TIME_UNITS["HOURS"], _TIME_UNITS["DAYS"]), key=mismatch)


def parse_experiment(case: Case, path: Path) -> None:
    raw_times = []
    for line in path.read_text().splitlines():
        tokens = line.split()
        # some files end with rows of zeros as padding
        if len(tokens) >= 3 and float(tokens[0]) > 0 and float(tokens[1]) > 0:
            raw_times.append(float(tokens[0]))
            case.exp_heat.append(float(tokens[1]))
            case.exp_uncertainty.append(float(tokens[2]))
    unit = _time_unit(raw_times, cumulative(case.cooling_steps))
    case.exp_time = [t * unit for t in raw_times]


def interpolate_log(times: list, values: list, t: float) -> float:
    """Log of the value is linear in time, NaN outside the range."""
    for i in range(len(times) - 1):
        t0, t1 = times[i], times[i + 1]
        if t0 <= t <= t1 and values[i] > 0 and values[i + 1] > 0:
            f = (t - t0) / (t1 - t0)
            return math.exp((1 - f) * math.log(values[i]) + f * math.log(values[i + 1]))
    return math.nan


def parse_fluxes(case: Case, path: Path) -> None:
    values = [float(v) for v in path.read_text().split()[:709]]
    # FISPACT-II gives high to low energy, NIS needs low to high
    case.flux_shape = list(reversed(values))


def load_cases(fns_dir: Path, materials: list) -> list:
    cases = []
    for material_dir in sorted(fns_dir.iterdir()):
        if not material_dir.is_dir():
            continue
        if materials and material_dir.name not in materials:
            continue
        for input_file in sorted(material_dir.glob("TENDL-2017_*.i")):
            name = input_file.stem.replace("TENDL-2017_", "")
            case = Case(material=material_dir.name, name=name)
            try:
                parse_input(case, input_file)
                parse_output(case, input_file.with_suffix(".out"))
                parse_heat(case, input_file.with_suffix(".nuclides"))
                parse_experiment(case, material_dir / f"{name}.exp")
                parse_fluxes(case, material_dir / f"{name}_fluxes")
            except (OSError, IndexError, ValueError) as e:
                print(f"skip {material_dir.name} {name}: {e}", file=sys.stderr)
                continue
            if not case.exp_time or sorted(case.exp_time) != case.exp_time:
                print(f"skip {material_dir.name} {name}: no measurements in time order", file=sys.stderr)
                continue
            cases.append(case)
    return cases


def run_case(case: Case, bounds: list, reaction_data, decay_mat, energy_data, times: list = None) -> list:
    """Decay heat (uW/g) at each cooling time (s, ascending), by default the measured times."""
    import pynis as nis
    times = case.exp_time if times is None else times
    # scale the shape to the total flux in the input, as FISPACT-II does
    shape_total = sum(case.flux_shape)
    flux = [v * case.flux_amplitude / shape_total for v in case.flux_shape]
    spectrum = nis.core.Spectrum(bounds, flux)
    irradiation_mat = decay_mat + nis.core.reaction_mat(reaction_data, spectrum)

    model = nis.core.BasicModel()
    model.nuclides = case.initial_atoms
    model.dt = case.irradiation_time
    result = nis.core.compute_inventory(model, irradiation_mat, nis.core.solver_cram48)
    shutdown_atoms = result.nuclides

    heats = []
    nuclides = shutdown_atoms
    steps = [t1 - t0 for t0, t1 in zip([0.0] + times[:-1], times)]
    for step in steps:
        model = nis.core.BasicModel()
        model.nuclides = nuclides
        model.dt = step
        nuclides = nis.core.compute_inventory(model, decay_mat, nis.core.solver_cram48).nuclides
        inventory = nis.core.vector_from_inventory(nuclides)
        heat = nis.core.compute_decay_heat(energy_data, inventory)
        heats.append(heat.total / case.mass_kg * _W_PER_KG_TO_UW_PER_G)
    return heats, shutdown_atoms


# set before the worker processes are forked, so they share the data and do not read it again
_SHARED = {}


def _run_shared(index: int):
    d = _SHARED
    case = d["cases"][index]
    times = d["times"](case) if d["times"] else None
    return run_case(case, d["bounds"], d["reaction_data"], d["decay_mat"], d["energy_data"], times)


def run_cases(cases: list, bounds: list, reaction_data, decay_mat, energy_data,
              jobs: int = 1, times=None):
    """run_case for each case, with jobs processes. times(case) gives the cooling times, if set.
    Yields the results in the same order as the cases."""
    _SHARED.update(cases=cases, bounds=bounds, reaction_data=reaction_data,
                   decay_mat=decay_mat, energy_data=energy_data, times=times)
    if jobs <= 1:
        yield from (_run_shared(i) for i in range(len(cases)))
        return
    # fork, so the workers get the data that is already read (Linux and macOS only)
    with multiprocessing.get_context("fork").Pool(jobs) as pool:
        yield from pool.imap(_run_shared, range(len(cases)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fns-dir", type=Path, required=True)
    parser.add_argument("--decay-data", type=Path, required=True)
    parser.add_argument("--reaction-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("fns_results.csv"))
    parser.add_argument("--nuclides-output", type=Path, default=None,
                        help="CSV of atoms after irradiation, NIS against FISPACT-II")
    parser.add_argument("--materials", nargs="*", default=[])
    parser.add_argument("--jobs", type=int, default=os.cpu_count(), help="number of cases to run at the same time")
    args = parser.parse_args()
    import pynis as nis

    cases = load_cases(args.fns_dir, args.materials)
    print(f"{len(cases)} cases", file=sys.stderr)

    start = time.time()
    bounds = read_bounds()
    reaction_data = nis.core.reaction_data_from_file(str(args.reaction_data))
    decay_mat = nis.core.decay_mat_from_file(str(args.decay_data))
    energy_data = nis.core.decay_energy_data_from_file(str(args.decay_data))
    print(f"data read in {time.time() - start:.1f} s", file=sys.stderr)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["material", "case", "cooling_time_s",
                         "exp_uW_per_g", "exp_unc_uW_per_g", "fispact_uW_per_g", "nis_uW_per_g"])
        nuclide_rows = []
        results = run_cases(cases, bounds, reaction_data, decay_mat, energy_data, args.jobs)
        for case, (heats, shutdown_atoms) in zip(cases, results):
            step_times = cumulative(case.cooling_steps)
            fispact = [interpolate_log(step_times, case.fispact_heat, t) for t in case.exp_time]
            for t, e, u, fh, nh in zip(case.exp_time, case.exp_heat, case.exp_uncertainty, fispact, heats):
                writer.writerow([case.material, case.name, f"{t:.6g}", e, u, f"{fh:.6e}", f"{nh:.6e}"])
            for nuclide, atoms in case.fispact_atoms.items():
                nuclide_rows.append([case.material, case.name, nuclide, atoms,
                                     shutdown_atoms.get(nuclide, 0.0)])
            ce = [h / e for h, e in zip(heats, case.exp_heat) if e > 0]
            fce = [h / e for h, e in zip(fispact, case.exp_heat) if e > 0 and not math.isnan(h)]
            fce_text = f"{sum(fce)/len(fce):7.3f}" if fce else "    n/a"
            print(f"{case.material:7s} {case.name:14s} NIS C/E {sum(ce)/len(ce):7.3f}  "
                  f"FISPACT-II C/E {fce_text}", file=sys.stderr)

    if args.nuclides_output:
        with open(args.nuclides_output, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["material", "case", "nuclide", "fispact_atoms", "nis_atoms"])
            writer.writerows(nuclide_rows)

    print(f"done in {time.time() - start:.1f} s, results in {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
