"""
Run the FNS decks from fns_make_decks.py with the nis program, and compare with the
experiment and FISPACT-II.

Each deck is one nis run (the irradiation step, then the cooling steps). The decks run in
parallel (--jobs). Each run reads the nuclear data again, so this is slower than
fns_benchmark.py, which reads the data once.

Writes the same CSV as fns_benchmark.py, so fns_stats.py can read it.

Usage:
    python scripts/fns_run_decks.py --nis build/bin/nis --decks benchmarks/fns \\
        --decay-data data/jeff311.dec --reaction-data data/jeff4.reac \\
        --output fns_decks.csv [--materials Co Fe ...] [--jobs 4]
"""

import argparse
import csv
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# 1 W/kg = 1e3 uW/g
_W_PER_KG_TO_UW_PER_G = 1e3


def read_mass(deck: Path) -> float:
    for line in deck.read_text().splitlines():
        if line.startswith("# mass_kg"):
            return float(line.split()[2])
    raise ValueError(f"no mass in {deck}")


def decay_heats(output: str) -> list:
    """The total decay heat (W) after each step in the nis output: irradiation, then each cooling step."""
    heats = []
    lines = output.splitlines()
    for i, line in enumerate(lines):
        if line == "decay_heat":
            heats.append(float(lines[i + 1].split()[1]))
    return heats


def run_deck(nis: Path, deck: Path, decay_data: Path, reaction_data: Path) -> list:
    # a large threshold, so the nis output has no nuclide list
    command = [str(nis), "-d", str(decay_data), "-r", str(reaction_data), "-i", str(deck), "-t", "1e300"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{deck}: {result.stderr.strip()}")
    # drop the irradiation step, the measured times are the cooling steps
    return decay_heats(result.stdout)[1:]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--nis", type=Path, required=True, help="the nis program")
    parser.add_argument("--decks", type=Path, default=Path("benchmarks/fns"))
    parser.add_argument("--decay-data", type=Path, required=True)
    parser.add_argument("--reaction-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("fns_decks.csv"))
    parser.add_argument("--materials", nargs="*", default=[])
    parser.add_argument("--jobs", type=int, default=os.cpu_count(), help="number of decks to run at the same time")
    args = parser.parse_args()

    decks = sorted(d for d in args.decks.glob("*/*.inp") if not args.materials or d.parent.name in args.materials)
    print(f"{len(decks)} decks", file=sys.stderr)

    start = time.time()
    with ThreadPoolExecutor(args.jobs) as pool:
        results = pool.map(lambda d: run_deck(args.nis, d, args.decay_data, args.reaction_data), decks)

        with open(args.output, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["material", "case", "cooling_time_s",
                             "exp_uW_per_g", "exp_unc_uW_per_g", "fispact_uW_per_g", "nis_uW_per_g"])
            for deck, heats in zip(decks, results):
                mass = read_mass(deck)
                with open(deck.with_suffix(".csv")) as ref:
                    rows = list(csv.DictReader(ref))
                if len(rows) != len(heats):
                    raise RuntimeError(f"{deck}: {len(heats)} cooling steps, {len(rows)} measured times")
                ce = []
                for r, heat in zip(rows, heats):
                    nis_heat = heat / mass * _W_PER_KG_TO_UW_PER_G
                    exp = float(r["exp_uW_per_g"])
                    writer.writerow([deck.parent.name, deck.stem, f"{float(r['cooling_time_s']):.6g}", r["exp_uW_per_g"],
                                     r["exp_unc_uW_per_g"], r["fispact_uW_per_g"], f"{nis_heat:.6e}"])
                    if exp > 0:
                        ce.append(nis_heat / exp)
                fce = [float(r["fispact_uW_per_g"]) / float(r["exp_uW_per_g"]) for r in rows
                       if float(r["exp_uW_per_g"]) > 0 and not math.isnan(float(r["fispact_uW_per_g"]))]
                fce_text = f"{sum(fce)/len(fce):7.3f}" if fce else "    n/a"
                print(f"{deck.parent.name:7s} {deck.stem:14s} NIS C/E {sum(ce)/len(ce):7.3f}  "
                      f"FISPACT-II C/E {fce_text}", file=sys.stderr)

    print(f"done in {time.time() - start:.1f} s, results in {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
