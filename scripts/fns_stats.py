"""
Summary statistics of the FNS benchmark, from the CSV written by fns_benchmark.py.

For each case, the C/E is the geometric mean over the measured times.
Only times where both NIS and FISPACT-II have a value are used, so the two codes
are compared on the same points.

Usage:
    python scripts/fns_stats.py fns_results.csv [more.csv ...] [--worst 25]
"""

import argparse
import collections
import csv
import math
import sys
from pathlib import Path


def read_ratios(paths: list) -> dict:
    """(material, case) -> {"nis": [C/E, ...], "fispact": [C/E, ...]}"""
    # a later file replaces the same (material, case, time) from an earlier one
    rows = {}
    for path in paths:
        with open(path) as f:
            for r in csv.DictReader(f):
                rows[(r["material"], r["case"], r["cooling_time_s"])] = r

    ratios = collections.defaultdict(lambda: {"nis": [], "fispact": []})
    for (material, case, _), r in rows.items():
        exp = float(r["exp_uW_per_g"])
        nis = float(r["nis_uW_per_g"])
        fispact = float(r["fispact_uW_per_g"])
        if exp <= 0 or nis <= 0 or fispact <= 0 or math.isnan(nis) or math.isnan(fispact):
            continue
        ratios[(material, case)]["nis"].append(nis / exp)
        ratios[(material, case)]["fispact"].append(fispact / exp)
    return ratios


def geometric_mean(values: list) -> float:
    return math.exp(sum(math.log(v) for v in values) / len(values))


def fraction_within(values: list, tolerance: float) -> float:
    # symmetric in log, so 0.8 and 1.25 are both within 25%
    return sum(1 for v in values if abs(math.log(v)) <= math.log(1 + tolerance)) / len(values)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="+", type=Path)
    parser.add_argument("--worst", type=int, default=25, help="number of cases to list where NIS and FISPACT-II differ most")
    args = parser.parse_args()

    ratios = read_ratios(args.csv)
    cases = sorted(ratios)
    if not cases:
        print("no cases found", file=sys.stderr)
        return 1

    case_nis = [geometric_mean(ratios[c]["nis"]) for c in cases]
    case_fispact = [geometric_mean(ratios[c]["fispact"]) for c in cases]
    point_nis = [v for c in cases for v in ratios[c]["nis"]]
    point_fispact = [v for c in cases for v in ratios[c]["fispact"]]

    print(f"{len(cases)} cases, {len(point_nis)} points")
    print(f"{'C/E within':>12s} {'cases NIS':>10s} {'FISPACT-II':>10s} {'points NIS':>11s} {'FISPACT-II':>10s}")
    for tolerance in (0.1, 0.2, 0.5):
        print(f"{tolerance:>12.0%} {fraction_within(case_nis, tolerance):>10.0%} "
              f"{fraction_within(case_fispact, tolerance):>10.0%} "
              f"{fraction_within(point_nis, tolerance):>11.0%} {fraction_within(point_fispact, tolerance):>10.0%}")

    nis_over_fispact = [n / f for n, f in zip(case_nis, case_fispact)]
    print(f"NIS/FISPACT-II within 10%: {fraction_within(nis_over_fispact, 0.1):.0%}, "
          f"within 20%: {fraction_within(nis_over_fispact, 0.2):.0%}")

    print(f"\ncases where NIS and FISPACT-II differ most:")
    worst = sorted(zip(cases, case_nis, case_fispact, nis_over_fispact), key=lambda x: -abs(math.log(x[3])))
    for (material, case), n, f, ratio in worst[:args.worst]:
        print(f"  {material:7s} {case:15s} NIS C/E {n:7.3f}  FISPACT-II C/E {f:7.3f}  NIS/FISPACT-II {ratio:7.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
