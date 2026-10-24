"""
Make NIS input decks for the FNS decay heat benchmark.

Reads the FISPACT-II validation inputs (the same as fns_benchmark.py) and writes, for each case:

    <output-dir>/<material>/<case>.inp   NIS input: irradiation (dt, nuclides, spectrum)
                                         and a cooling step to each measured time
    <output-dir>/<material>/<case>.csv   measured decay heat and FISPACT-II (TENDL-2017)
                                         at the same times, in uW/g

The foil mass is in a comment in the deck ("# mass_kg ..."), fns_run_decks.py uses it
to give the heat per gram. The flux is scaled to the total flux of the FISPACT-II input.

Usage:
    python scripts/fns_make_decks.py --fns-dir path/to/fns --output-dir benchmarks/fns [--materials Co Fe ...]
"""

import argparse
import csv
import sys
from pathlib import Path

from fns_benchmark import cumulative, interpolate_log, load_cases, read_bounds


def _numbers(values: list) -> str:
    return " ".join(f"{v:.9e}" for v in values)


def write_deck(case, bounds: list, path: Path) -> None:
    shape_total = sum(case.flux_shape)
    flux = [v * case.flux_amplitude / shape_total for v in case.flux_shape]
    steps = [t1 - t0 for t0, t1 in zip([0.0] + case.exp_time[:-1], case.exp_time)]

    lines = [
        f"# FNS decay heat benchmark, {case.material} {case.name}",
        "# made by scripts/fns_make_decks.py from the FISPACT-II validation inputs",
        f"# mass_kg {case.mass_kg:.9e}",
        f"# total flux {case.flux_amplitude:.9e} n/cm2/s, 709 groups (CCFE-709)",
        "# cooling: one step to each measured time after shutdown",
        f"dt {case.irradiation_time:.9e}",
        "",
        "nuclides",
    ]
    lines += [f"    {name} {atoms:.9e}" for name, atoms in sorted(case.initial_atoms.items())]
    lines += ["", "spectrum", f"    {_numbers(bounds)}", f"    {_numbers(flux)}"]
    lines += ["", "cooling", f"    {_numbers(steps)}", ""]
    path.write_text("\n".join(lines))


def write_reference(case, path: Path) -> None:
    step_times = cumulative(case.cooling_steps)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cooling_time_s", "exp_uW_per_g", "exp_unc_uW_per_g", "fispact_uW_per_g"])
        for t, e, u in zip(case.exp_time, case.exp_heat, case.exp_uncertainty):
            writer.writerow([f"{t:.9e}", e, u, f"{interpolate_log(step_times, case.fispact_heat, t):.6e}"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fns-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("benchmarks/fns"))
    parser.add_argument("--materials", nargs="*", default=[])
    args = parser.parse_args()

    bounds = read_bounds()
    cases = load_cases(args.fns_dir, args.materials)
    for case in cases:
        directory = args.output_dir / case.material
        directory.mkdir(parents=True, exist_ok=True)
        write_deck(case, bounds, directory / f"{case.name}.inp")
        write_reference(case, directory / f"{case.name}.csv")
    print(f"{len(cases)} decks in {args.output_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
