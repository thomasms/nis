"""
One plot of the C/E of all the FNS cases, NIS against FISPACT-II.

For each case, the point is the geometric mean of C/E over the measured times (the same as fns_stats.py,
only the times where NIS and FISPACT-II both have a value). The gray band is C/E within 20 %. A value
outside the plot range is shown at the edge as an arrow, with its value.

Reads the CSV of fns_benchmark.py or fns_run_decks.py. Needs matplotlib.

Usage:
    python scripts/fns_summary_plot.py fns.csv --output fns_summary.png [--columns 2]
"""

import argparse
import math
import sys
from pathlib import Path

from fns_stats import fraction_within, geometric_mean, read_ratios

# the default categorical palette, slots 1 and 2 (validated for colour blindness)
NIS_COLOUR = "#2a78d6"
FISPACT_COLOUR = "#eb6834"
TEXT_COLOUR = "#3a3a38"
MUTED_COLOUR = "#8a8984"
GRID_COLOUR = "#e4e3df"
BAND_COLOUR = "#efeeea"
SURFACE_COLOUR = "#fcfcfb"

X_MIN, X_MAX = 0.1, 10.0
ROW_OFFSET = 0.17


def _plot_column(ax, rows: list, show_legend: bool) -> None:
    """rows: [(label, nis C/E, fispact C/E), ...], top to bottom."""
    y = list(range(len(rows)))[::-1]
    ax.axvspan(1 / 1.2, 1.2, color=BAND_COLOUR, zorder=0, lw=0)
    ax.axvline(1.0, color=MUTED_COLOUR, lw=1, zorder=1)

    # NIS a little above the row, FISPACT-II a little below, so a point never hides the other
    for (label, nis, fispact), yi in zip(rows, y):
        # a thin line between the two codes, so the difference is easy to see
        ax.plot([_clip(nis), _clip(fispact)], [yi + ROW_OFFSET, yi - ROW_OFFSET], color=GRID_COLOUR, lw=1.5, zorder=2)

    for values, colour, marker, name, offset in (([r[1] for r in rows], NIS_COLOUR, "o", "NIS", ROW_OFFSET),
                                                 ([r[2] for r in rows], FISPACT_COLOUR, "D", "FISPACT-II", -ROW_OFFSET)):
        inside = [(v, yi + offset) for v, yi in zip(values, y) if X_MIN <= v <= X_MAX]
        ax.scatter([v for v, _ in inside], [yi for _, yi in inside], s=26, marker=marker, color=colour,
                   edgecolors=SURFACE_COLOUR, linewidths=0.8, zorder=3, label=name)
        for v, yi in zip(values, y):
            if v < X_MIN or v > X_MAX:
                # an arrow at the edge, and the value inside the plot next to it
                edge = X_MIN if v < X_MIN else X_MAX
                ax.scatter([edge], [yi + offset], s=40, marker="<" if v < X_MIN else ">", color=colour,
                           edgecolors=SURFACE_COLOUR, linewidths=0.8, zorder=3, clip_on=False)
                # the NIS value a little up, FISPACT-II a little down, so both can be read
                ax.annotate(f"{v:.2g}", (edge, yi + offset), xytext=(7 if v < X_MIN else -7, 3.5 if offset > 0 else -3.5),
                            textcoords="offset points", ha="left" if v < X_MIN else "right", va="center",
                            fontsize=7, color=TEXT_COLOUR)

    ax.set_xscale("log")
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7, color=TEXT_COLOUR)
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5, 10])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5", "10"], fontsize=8, color=TEXT_COLOUR)
    ax.minorticks_off()
    ax.set_xlabel("C/E (geometric mean over the measured times)", fontsize=9, color=TEXT_COLOUR)
    ax.grid(axis="x", color=GRID_COLOUR, lw=0.8, zorder=0)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor(SURFACE_COLOUR)
    if show_legend:
        ax.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=2, frameon=False, fontsize=9,
                  handletextpad=0.3, columnspacing=1.2, borderaxespad=0.2)


def _clip(v: float) -> float:
    return min(max(v, X_MIN), X_MAX)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="+", type=Path, help="results of fns_benchmark.py or fns_run_decks.py")
    parser.add_argument("--output", type=Path, default=Path("fns_summary.png"))
    parser.add_argument("--columns", type=int, default=2)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ratios = read_ratios(args.csv)
    cases = [c for c in sorted(ratios) if ratios[c]["nis"]]
    if not cases:
        print("no cases found", file=sys.stderr)
        return 1
    rows = [(f"{material} {case}", geometric_mean(ratios[(material, case)]["nis"]),
             geometric_mean(ratios[(material, case)]["fispact"])) for material, case in cases]

    per_column = math.ceil(len(rows) / args.columns)
    columns = [rows[i:i + per_column] for i in range(0, len(rows), per_column)]
    fig, axes = plt.subplots(1, len(columns), figsize=(5.2 * len(columns), 0.16 * per_column + 1.6),
                             squeeze=False, facecolor=SURFACE_COLOUR)
    for i, (ax, column) in enumerate(zip(axes[0], columns)):
        _plot_column(ax, column, show_legend=(i == 0))

    nis_within = fraction_within([r[1] for r in rows], 0.2)
    fispact_within = fraction_within([r[2] for r in rows], 0.2)
    height = fig.get_figheight()
    fig.text(0.01, 1 - 0.25 / height, "FNS decay heat benchmark: C/E for each case",
             fontsize=12, color=TEXT_COLOUR, ha="left", va="top")
    fig.text(0.01, 1 - 0.5 / height,
             f"{len(rows)} cases. Within 20 % of the experiment (gray band, 0.83 to 1.2): "
             f"NIS {nis_within:.0%} of cases, FISPACT-II {fispact_within:.0%}.",
             fontsize=9, color=MUTED_COLOUR, ha="left", va="top")
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.75 / height))
    fig.savefig(args.output, dpi=150, facecolor=SURFACE_COLOUR)
    print(f"plot in {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
