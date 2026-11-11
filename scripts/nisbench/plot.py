"""
Plots of a results file. Needs matplotlib.

    summary   one plot of the C/E of all the cases of a benchmark and quantity: NIS and each calculation
    cases     one plot for each case: the values against cooling time, and C/E
"""

import collections
import math
from pathlib import Path

from .stats import case_ratios, fraction_within, geometric_mean

# the default categorical palette (validated for colour blindness): NIS, then the calculations
CODE_COLOURS = ["#2a78d6", "#eb6834", "#1f9e74", "#a35fc7"]
CODE_MARKERS = ["o", "D", "s", "^"]
INK = "#0b0b0b"
TEXT = "#3a3a38"
MUTED = "#8a8984"
GRID = "#e4e3df"
BAND = "#efeeea"
SURFACE = "#fcfcfb"

X_MIN, X_MAX = 0.1, 10.0

_QUANTITY_LABELS = {"decay_heat": "Decay heat", "beta_heat": "Beta decay heat", "gamma_heat": "Gamma decay heat",
                    "alpha_heat": "Alpha decay heat", "activity": "Activity"}


def _pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _label(code: str) -> str:
    return "NIS" if code == "nis" else code


def _clip(v: float) -> float:
    return min(max(v, X_MIN), X_MAX)


def _summary_column(ax, rows: list, codes: list, show_legend: bool) -> None:
    """rows: [(label, {code: C/E}), ...], top to bottom."""
    y = list(range(len(rows)))[::-1]
    ax.axvspan(1 / 1.2, 1.2, color=BAND, zorder=0, lw=0)
    ax.axvline(1.0, color=MUTED, lw=1, zorder=1)

    # the codes are spread over the row, so a point never hides another
    spread = 0.34
    offsets = [spread / 2 - spread * i / (len(codes) - 1) for i in range(len(codes))] if len(codes) > 1 else [0.0]
    if len(codes) > 1:
        for (_, values), yi in zip(rows, y):
            present = [(_clip(values[c]), yi + o) for c, o in zip(codes, offsets) if c in values]
            ax.plot([p[0] for p in present], [p[1] for p in present], color=GRID, lw=1.5, zorder=2)

    for i, (code, offset) in enumerate(zip(codes, offsets)):
        colour, marker = CODE_COLOURS[i % len(CODE_COLOURS)], CODE_MARKERS[i % len(CODE_MARKERS)]
        points = [(values[code], yi + offset) for (_, values), yi in zip(rows, y) if code in values]
        inside = [(v, yy) for v, yy in points if X_MIN <= v <= X_MAX]
        ax.scatter([v for v, _ in inside], [yy for _, yy in inside], s=26, marker=marker, color=colour,
                   edgecolors=SURFACE, linewidths=0.8, zorder=3, label=_label(code))
        for v, yy in points:
            if X_MIN <= v <= X_MAX:
                continue
            # an arrow at the edge, and the value inside the plot next to it
            low = v < X_MIN
            edge = X_MIN if low else X_MAX
            ax.scatter([edge], [yy], s=40, marker="<" if low else ">", color=colour,
                       edgecolors=SURFACE, linewidths=0.8, zorder=3, clip_on=False)
            ax.annotate(f"{v:.2g}", (edge, yy), xytext=(7 if low else -7, 3.5 if offset >= 0 else -3.5),
                        textcoords="offset points", ha="left" if low else "right", va="center",
                        fontsize=7, color=TEXT)

    ax.set_xscale("log")
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7, color=TEXT)
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5, 10])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5", "10"], fontsize=8, color=TEXT)
    ax.minorticks_off()
    ax.set_xlabel("C/E (geometric mean over the measured times)", fontsize=9, color=TEXT)
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor(SURFACE)
    if show_legend:
        ax.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=len(codes), frameon=False, fontsize=9,
                  handletextpad=0.3, columnspacing=1.2, borderaxespad=0.2)


def summary(rows: list, calculations: list, benchmark: str, quantity: str, output: Path,
            columns: int = 2, title: str = None) -> None:
    ratios = case_ratios(rows, calculations).get((benchmark, quantity))
    if not ratios:
        raise ValueError(f"no results for {benchmark}, {quantity}")
    codes = ["nis"] + [c for c in calculations if any(c in r for r in ratios.values())]
    table = [(case.replace("/", " "), {code: geometric_mean(v) for code, v in r.items()}) for case, r in ratios.items()]

    plt = _pyplot()
    columns = max(1, min(columns, len(table)))
    per_column = math.ceil(len(table) / columns)
    parts = [table[i:i + per_column] for i in range(0, len(table), per_column)]
    fig, axes = plt.subplots(1, len(parts), figsize=(5.2 * len(parts), 0.16 * per_column + 1.6),
                             squeeze=False, facecolor=SURFACE)
    for i, (ax, part) in enumerate(zip(axes[0], parts)):
        _summary_column(ax, part, codes, show_legend=(i == 0))

    within = ", ".join(f"{_label(c)} {fraction_within([r[1][c] for r in table if c in r[1]], 0.2):.0%}" for c in codes)
    height = fig.get_figheight()
    fig.text(0.01, 1 - 0.25 / height, title or f"{benchmark}: C/E of the {_QUANTITY_LABELS[quantity].lower()} for each case",
             fontsize=12, color=TEXT, ha="left", va="top")
    fig.text(0.01, 1 - 0.5 / height,
             f"{len(table)} cases. Cases within 20 % of the experiment (gray band, 0.83 to 1.2): {within}.",
             fontsize=9, color=MUTED, ha="left", va="top")
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.75 / height))
    fig.savefig(output, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def _unit_for(times: list) -> tuple:
    """(seconds per unit, label) to plot the times in"""
    if times[-1] <= 600:
        return 1.0, "s"
    if times[-1] <= 2 * 86400:
        return 60.0, "min"
    return 86400.0, "d"


def cases(rows: list, calculations: list, output_dir: Path, file_format: str = "png") -> list:
    """One plot for each (benchmark, case, quantity). Returns the files."""
    plt = _pyplot()
    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r["benchmark"], r["case"], r["quantity"])].append(r)

    files = []
    output_dir.mkdir(parents=True, exist_ok=True)
    for (benchmark, case, quantity), group in sorted(groups.items()):
        group.sort(key=lambda r: r["time_s"])
        times = [r["time_s"] for r in group]
        scale, unit_label = _unit_for(times)
        t = [x / scale for x in times]
        codes = ["nis"] + [c for c in calculations if any(not math.isnan(r[c]) for r in group)]

        fig, (ax, ax_ce) = plt.subplots(2, 1, sharex=True, figsize=(7, 6),
                                        gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.08})
        exp = [r["experiment"] for r in group]
        ax.errorbar(t, exp, yerr=[r["uncertainty"] for r in group], fmt="o", ms=4, color=INK, ecolor=INK,
                    elinewidth=1, capsize=2, label="Experiment", zorder=3)
        band = [r["uncertainty"] / r["experiment"] if r["experiment"] > 0 else 0 for r in group]
        ax_ce.fill_between(t, [1 - b for b in band], [1 + b for b in band], color=BAND, lw=0, zorder=1)
        ax_ce.axhline(1.0, color=MUTED, lw=1, zorder=1)

        all_ce = []
        for i, code in enumerate(codes):
            colour, marker = CODE_COLOURS[i % len(CODE_COLOURS)], CODE_MARKERS[i % len(CODE_MARKERS)]
            values = [r[code] for r in group]
            ax.plot(t, values, color=colour, lw=2, ls="-" if i == 0 else "--", marker=marker, ms=3,
                    label=_label(code), zorder=2)
            ce = [v / e if e > 0 else math.nan for v, e in zip(values, exp)]
            ax_ce.plot(t, ce, marker, ms=4, color=colour, zorder=3)
            all_ce += [v for v in ce if v > 0 and not math.isnan(v)]

        unit = group[0]["unit"].replace("uW", "µW")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel(f"{_QUANTITY_LABELS[quantity]} ({unit})", color=INK)
        ax.set_title(f"{benchmark} {case}", color=INK, loc="left")
        ax.legend(frameon=False, labelcolor=INK)
        ax_ce.set_ylabel("C/E", color=INK)
        ax_ce.set_xlabel(f"Cooling time ({unit_label})", color=INK)
        if all_ce:
            low, high = min(all_ce + [0.8]), max(all_ce + [1.2])
            if high / low > 5:
                ax_ce.set_yscale("log")
            ax_ce.set_ylim(low * 0.9, high * 1.1)
        for a in (ax, ax_ce):
            a.grid(True, which="major", color=GRID, lw=0.8)
            a.set_axisbelow(True)
            a.tick_params(colors=MUTED)
            for spine in a.spines.values():
                spine.set_color(GRID)

        suffix = "" if quantity == "decay_heat" else f"_{quantity}"
        path = output_dir / f"{benchmark}_{case.replace('/', '_')}{suffix}.{file_format}"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        files.append(path)
    return files
