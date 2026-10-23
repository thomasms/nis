"""
Find the nuclides that cause a difference in decay heat between NIS and FISPACT-II.

Uses the nuclides CSV from fns_benchmark.py (atoms at shutdown for both codes) and
the decay data, to give the heat of each nuclide at shutdown as a fraction of the
FISPACT-II total. Both codes use the same decay energies here, so a difference
comes from the inventory only.

Usage:
    export PYTHONPATH=./build/bin/python3
    python scripts/fns_heat_by_nuclide.py fns_nuclides.csv --decay-data data/jeff311.dec \\
        [--cases Ce:2000exp_5min Mn:2000exp_5min] [--top 5] [--tolerance 0.2]

With no --cases, it lists each case where the NIS heat at shutdown is not
within the tolerance of FISPACT-II.
"""

import argparse
import collections
import csv
import math
import sys
from pathlib import Path

import pynis as nis

_EV_TO_J = 1.602176634e-19


def read_atoms(path: Path) -> dict:
    """(material, case) -> {nuclide: (fispact atoms, nis atoms)}"""
    atoms = collections.defaultdict(dict)
    with open(path) as f:
        for r in csv.DictReader(f):
            atoms[(r["material"], r["case"])][r["nuclide"]] = (
                float(r["fispact_atoms"] or 0.0), float(r["nis_atoms"] or 0.0))
    return atoms


def heat_per_atom(decay_data: Path) -> dict:
    """ZAI -> decay heat per atom in W"""
    energy_data = nis.core.decay_energy_data_from_file(str(decay_data))
    energies = energy_data.light_particle + energy_data.electromagnetic + energy_data.heavy_particle
    heat = energy_data.decay_constants * energies * _EV_TO_J
    return {zai: heat[i] for i, zai in enumerate(nis.data.ALL_ZAIS())}


def nuclide_heats(case_atoms: dict, heats: dict) -> list:
    """[(nuclide, FISPACT-II heat, NIS heat, FISPACT-II atoms, NIS atoms)]"""
    out = []
    for nuclide, (fispact, nis_atoms) in case_atoms.items():
        zai = nis.zai.get_from_name(nuclide)
        if zai not in heats:
            continue
        out.append((nuclide, fispact * heats[zai], nis_atoms * heats[zai], fispact, nis_atoms))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("nuclides_csv", type=Path)
    parser.add_argument("--decay-data", type=Path, required=True)
    parser.add_argument("--cases", nargs="*", default=[], help="material:case, e.g. Ce:2000exp_5min")
    parser.add_argument("--top", type=int, default=5, help="number of nuclides to show for each case")
    parser.add_argument("--tolerance", type=float, default=0.2)
    args = parser.parse_args()

    atoms = read_atoms(args.nuclides_csv)
    heats = heat_per_atom(args.decay_data)

    keys = [tuple(c.split(":", 1)) for c in args.cases] if args.cases else sorted(atoms)
    for key in keys:
        if key not in atoms:
            print(f"{key[0]}:{key[1]} not found", file=sys.stderr)
            continue
        rows = nuclide_heats(atoms[key], heats)
        fispact_total = sum(r[1] for r in rows)
        nis_total = sum(r[2] for r in rows)
        if fispact_total <= 0:
            continue
        ratio = nis_total / fispact_total
        if not args.cases and ratio > 0 and abs(math.log(ratio)) <= math.log(1 + args.tolerance):
            continue

        print(f"{key[0]} {key[1]}: heat at shutdown NIS/FISPACT-II {ratio:.3f}")
        rows.sort(key=lambda r: -max(r[1], r[2]))
        for nuclide, fispact_heat, nis_heat, fispact_atoms, nis_atoms in rows[:args.top]:
            atoms_ratio = f"{nis_atoms / fispact_atoms:.3g}" if fispact_atoms > 0 else "inf"
            print(f"    {nuclide:8s} FISPACT-II {fispact_heat / fispact_total:6.1%}  "
                  f"NIS {nis_heat / fispact_total:6.1%}  atoms NIS/FISPACT-II {atoms_ratio}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
