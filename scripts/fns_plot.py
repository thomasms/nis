"""
Plot the FNS decay heat benchmark: NIS and FISPACT-II against the experiment.

For each case, the top panel shows the decay heat against cooling time:
the measured values with their uncertainty, NIS as a smooth curve (it is run at
log spaced times and at the measured times) and FISPACT-II at its cooling steps.
The bottom panel shows C/E at the measured times, with the experimental
uncertainty as a band around 1.

Uses the same inputs as fns_benchmark.py. Each case needs one CRAM solve per time, so the cases
run in parallel (--jobs), and --materials plots only some of them.

Usage:
    export PYTHONPATH=./build/bin/python3
    python scripts/fns_plot.py --fns-dir path/to/fns \\
        --decay-data data/jeff311.dec --reaction-data data/jeff4.reac \\
        --output-dir fns_plots [--materials Co Fe ...] [--points 60] [--format png] [--jobs 4]
"""

import argparse
import math
import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pynis as nis

from fns_benchmark import cumulative, interpolate_log, load_cases, reaction_cache, read_bounds, run_cases

# colours: the experiment is ink, the codes are categorical slots 1 and 2
_INK = "#0b0b0b"
_MUTED = "#52514e"
_GRID = "#e4e3df"
_BAND = "#ecebe7"
_NIS = "#2a78d6"
_FISPACT = "#eb6834"


def smooth_times(case, points: int) -> list:
    """Log spaced cooling times over the measured range, plus the measured times."""
    t0 = case.exp_time[0] * 0.9
    t1 = case.exp_time[-1] * 1.1
    grid = [t0 * (t1 / t0) ** (i / (points - 1)) for i in range(points)]
    return sorted(set(grid) | set(case.exp_time))


def time_unit(times: list) -> tuple:
    """(seconds per unit, label) to plot the times in"""
    if times[-1] <= 2 * 86400:
        return 60.0, "min"
    return 86400.0, "d"


def plot_case(case, times: list, nis_heat: list, output: Path) -> None:
    unit, unit_label = time_unit(case.exp_time)
    exp_t = [t / unit for t in case.exp_time]

    step_times = cumulative(case.cooling_steps)
    fispact = [(t, h) for t, h in zip(step_times, case.fispact_heat)
               if exp_t[0] * unit * 0.5 <= t <= exp_t[-1] * unit * 2 and h > 0]

    nis_at = dict(zip(times, nis_heat))
    nis_ce = [nis_at[t] / e for t, e in zip(case.exp_time, case.exp_heat)]
    fispact_ce = [interpolate_log(step_times, case.fispact_heat, t) / e
                  for t, e in zip(case.exp_time, case.exp_heat)]

    fig, (ax, ax_ce) = plt.subplots(2, 1, sharex=True, figsize=(7, 6),
                                    gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.08})

    ax.errorbar(exp_t, case.exp_heat, yerr=case.exp_uncertainty, fmt="o", ms=4,
                color=_INK, ecolor=_INK, elinewidth=1, capsize=2, label="Experiment", zorder=3)
    ax.plot([t / unit for t in times], nis_heat, color=_NIS, lw=2, label="NIS", zorder=2)
    if fispact:
        ax.plot([t / unit for t, _ in fispact], [h for _, h in fispact], color=_FISPACT, lw=2,
                ls="--", marker="s", ms=3, label="FISPACT-II (TENDL-2017)", zorder=2)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Decay heat (µW/g)", color=_INK)
    ax.set_title(f"FNS {case.material}, {case.name}", color=_INK, loc="left")
    ax.legend(frameon=False, labelcolor=_INK)

    # experimental uncertainty as a band around C/E = 1
    band = [u / e for u, e in zip(case.exp_uncertainty, case.exp_heat)]
    ax_ce.fill_between(exp_t, [1 - b for b in band], [1 + b for b in band], color=_BAND, lw=0, zorder=1)
    ax_ce.axhline(1.0, color=_MUTED, lw=1, zorder=1)
    ax_ce.plot(exp_t, nis_ce, "o", ms=4, color=_NIS, label="NIS", zorder=3)
    ax_ce.plot(exp_t, fispact_ce, "s", ms=4, color=_FISPACT, label="FISPACT-II", zorder=3)
    ax_ce.set_ylabel("C/E", color=_INK)
    ax_ce.set_xlabel(f"Cooling time ({unit_label})", color=_INK)

    finite = [v for v in nis_ce + fispact_ce if v > 0 and not math.isnan(v)]
    if finite:
        low, high = min(finite + [0.8]), max(finite + [1.2])
        if high / low > 5:
            ax_ce.set_yscale("log")
        ax_ce.set_ylim(low * 0.9, high * 1.1)

    for a in (ax, ax_ce):
        a.grid(True, which="major", color=_GRID, lw=0.8)
        a.set_axisbelow(True)
        a.tick_params(colors=_MUTED)
        for spine in a.spines.values():
            spine.set_color(_GRID)

    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fns-dir", type=Path, required=True)
    parser.add_argument("--decay-data", type=Path, required=True)
    parser.add_argument("--reaction-data", type=Path, required=True,
                        help="ENDF reaction data or a cache of it, the cache is made first if needed")
    parser.add_argument("--output-dir", type=Path, default=Path("fns_plots"))
    parser.add_argument("--materials", nargs="*", default=[])
    parser.add_argument("--points", type=int, default=60, help="number of log spaced times for the NIS curve")
    parser.add_argument("--format", default="png", help="file format for matplotlib, e.g. png, pdf, svg")
    parser.add_argument("--jobs", type=int, default=os.cpu_count(), help="number of cases to run at the same time")
    args = parser.parse_args()

    cases = load_cases(args.fns_dir, args.materials)
    print(f"{len(cases)} cases", file=sys.stderr)

    start = time.time()
    bounds = read_bounds()
    reaction_data = nis.core.reaction_data_from_file(str(reaction_cache(args.reaction_data)))
    decay_mat = nis.core.decay_mat_from_file(str(args.decay_data))
    energy_data = nis.core.decay_energy_data_from_file(str(args.decay_data))
    print(f"data read in {time.time() - start:.1f} s", file=sys.stderr)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    def times(case):
        return smooth_times(case, args.points)

    results = run_cases(cases, bounds, reaction_data, decay_mat, energy_data, args.jobs, times)
    for case, (heats, _) in zip(cases, results):
        output = args.output_dir / f"{case.material}_{case.name}.{args.format}"
        plot_case(case, times(case), heats, output)
        print(f"{output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
