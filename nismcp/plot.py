"""PNG plots of a run: a quantity of some nuclides against the time."""

import io

import numpy as np

from .model import QUANTITIES, nuclide_values
from .units import DAY, HOUR, MINUTE, YEAR

# the categorical colours, always in this order (light mode); a 9th series is not allowed
COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MAX_SERIES = len(COLOURS)
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e4e3df"
IRRADIATION = "#ebeae6"

_LABELS = {"activity": "Activity (Bq)", "decay_heat": "Decay heat (W)", "atoms": "Atoms", "mass": "Mass (g)"}


def _time_unit(seconds: float) -> tuple:
    for name, size in (("y", YEAR), ("d", DAY), ("h", HOUR), ("min", MINUTE)):
        if seconds >= 3 * size:
            return name, size
    return "s", 1.0


def plot_png(result, tables, names: list, quantity: str = "activity", log_time: bool = None,
             log_value: bool = True, title: str = None) -> bytes:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if quantity not in QUANTITIES:
        raise ValueError(f"unknown quantity {quantity!r}, use one of {sorted(QUANTITIES)}")
    if not names:
        raise ValueError("no nuclides to plot")
    if len(names) > MAX_SERIES:
        raise ValueError(f"plot at most {MAX_SERIES} nuclides, so each one has its own colour")

    points = [p for p in result.points if p.time > 0]
    times = np.array([p.time for p in points])
    unit, size = _time_unit(times[-1])
    x = times / size
    values = np.array([nuclide_values(tables, p.atoms, quantity) for p in points])
    if log_time is None:
        log_time = times[-1] / times[0] > 100

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": TEXT_SECONDARY, "axes.labelcolor": TEXT,
                         "xtick.color": TEXT_SECONDARY, "ytick.color": TEXT_SECONDARY, "text.color": TEXT})
    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # the irradiation steps, as shaded bands under the lines
    t0 = 0.0
    shaded = False
    for step in result.steps:
        if step.flux > 0:
            ax.axvspan(max(t0, times[0]) / size, (t0 + step.duration) / size, color=IRRADIATION, lw=0, zorder=0,
                       label=None if shaded else "irradiation")
            shaded = True
        t0 += step.duration

    plotted = []
    for colour, name in zip(COLOURS, names):
        y = values[:, tables.index[name]]
        if log_value:
            y = np.where(y > 0, y, np.nan)
        empty = np.all(np.isnan(y)) if log_value else not np.any(y)
        if empty:
            continue
        ax.plot(x, y, color=colour, lw=2, label=name, zorder=3, solid_capstyle="round")
        plotted.append((name, colour, y))
    if not plotted:
        plt.close(fig)
        raise ValueError(f"the nuclides {', '.join(names)} have no {quantity} in this run")

    if log_time:
        ax.set_xscale("log")
    if log_value:
        ax.set_yscale("log")
    ax.set_xlabel(f"Time since the start of the irradiation ({unit})")
    ax.set_ylabel(_LABELS[quantity])
    ax.grid(True, which="major", color=GRID, lw=0.8, zorder=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # a direct label at the end of each line, when there are few lines
    if len(plotted) <= 4:
        for name, colour, y in plotted:
            last = np.flatnonzero(~np.isnan(y)) if log_value else np.arange(len(y))
            if len(last):
                i = last[-1]
                ax.annotate(name, (x[i], y[i]), xytext=(6, 0), textcoords="offset points", va="center",
                            fontsize=9, color=TEXT)
        ax.margins(x=0.08)
    ax.legend(frameon=False, loc="best", fontsize=9)
    ax.set_title(title or f"{_LABELS[quantity].split(' (')[0]} of {', '.join(n for n, _, _ in plotted)}",
                 loc="left", fontsize=11, color=TEXT)
    fig.tight_layout()
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", facecolor=SURFACE)
    plt.close(fig)
    return buffer.getvalue()
