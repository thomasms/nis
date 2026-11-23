"""
C/E statistics of a results file, for each benchmark and quantity.

For each case, the C/E is the geometric mean over the measured times. Only the times where NIS and all
the calculations of the case have a value are used, so the codes are compared on the same points.

report() gives tables with borders. With color=True, the values have colours: green for a C/E within 10% of 1
(or two runs within 10%), yellow within 50%, red for more; the best value of a row is bold, and the cases where the
experiment has a large uncertainty are yellow.
"""

import collections
import math
import re
import statistics

# a case with a larger median uncertainty of the experiment than this is marked as uncertain
UNCERTAIN = 0.5


def geometric_mean(values: list) -> float:
    return math.exp(sum(math.log(v) for v in values) / len(values))


def fraction_within(values: list, tolerance: float) -> float:
    # symmetric in log, so 0.8 and 1.25 are both within 25%
    return sum(1 for v in values if abs(math.log(v)) <= math.log(1 + tolerance)) / len(values)


def _ok(v: float) -> bool:
    return v is not None and not math.isnan(v) and v > 0


def case_ratios(rows: list, calculations: list) -> dict:
    """(benchmark, quantity) -> {case: {code: [C/E, ...]}}, code is "nis" or a calculation name."""
    by_case = collections.defaultdict(list)
    for r in rows:
        by_case[(r["benchmark"], r["quantity"], r["case"])].append(r)

    out = collections.defaultdict(dict)
    for (benchmark, quantity, case), case_rows in sorted(by_case.items()):
        codes = ["nis"] + [c for c in calculations if any(not math.isnan(r[c]) for r in case_rows)]
        ratios = {code: [] for code in codes}
        for r in case_rows:
            if not _ok(r["experiment"]) or not all(_ok(r[c]) for c in codes):
                continue
            for code in codes:
                ratios[code].append(r[code] / r["experiment"])
        if ratios["nis"]:
            out[(benchmark, quantity)][case] = ratios
    return out


def case_uncertainties(rows: list) -> dict:
    """(benchmark, quantity, case) -> the median relative uncertainty of the experiment (nan if none)."""
    values = collections.defaultdict(list)
    for r in rows:
        u = r.get("uncertainty", math.nan)
        if _ok(r["experiment"]) and u is not None and not math.isnan(u):
            values[(r["benchmark"], r["quantity"], r["case"])].append(u / r["experiment"])
    return {k: statistics.median(v) for k, v in values.items()}


# ---------------------------------------------------------------- the text

_COLORS = {"green": "32", "yellow": "33", "red": "31", "bold": "1", "dim": "2", "cyan": "36"}
_ANSI = re.compile(r"\033\[[0-9;]*m")
_BOX = {"unicode": dict(h="─", v="│", tl="┌", tm="┬", tr="┐", ml="├", mm="┼", mr="┤", bl="└", bm="┴", br="┘"),
        "ascii": dict(h="-", v="|", tl="+", tm="+", tr="+", ml="+", mm="+", mr="+", bl="+", bm="+", br="+")}


class _Style:
    def __init__(self, color: bool, unicode: bool):
        self.color, self.box = color, _BOX["unicode" if unicode else "ascii"]

    def __call__(self, text: str, *styles: str) -> str:
        styles = [s for s in styles if s]
        if not self.color or not styles:
            return text
        return f"\033[{';'.join(_COLORS[s] for s in styles)}m{text}\033[0m"


def _width(text: str) -> int:
    return len(_ANSI.sub("", text))


def _pad(text: str, width: int, align: str) -> str:
    space = " " * (width - _width(text))
    return text + space if align == "<" else space + text if align == ">" else \
        space[:len(space) // 2] + text + space[len(space) // 2:]


def _table(style: _Style, headers: list, rows: list, align: list, title: str = None) -> list:
    """A table with borders. headers: rows of header cells; rows: rows of cells; align: '<', '>' or '^'."""
    b = style.box
    widths = [max(_width(r[i]) for r in headers + rows) for i in range(len(align))]

    def line(left, middle, right):
        return left + middle.join(b["h"] * (w + 2) for w in widths) + right

    def cells(row, aligns):
        return b["v"] + b["v"].join(f" {_pad(c, w, a)} " for c, w, a in zip(row, widths, aligns)) + b["v"]

    out = [style(title, "bold")] if title else []
    out.append(line(b["tl"], b["tm"], b["tr"]))
    for h in headers:
        out.append(cells([style(c, "bold") for c in h], ["<" if a == "<" else "^" for a in align]))
    out.append(line(b["ml"], b["mm"], b["mr"]))
    out += [cells(r, align) for r in rows]
    out.append(line(b["bl"], b["bm"], b["br"]))
    return out


def _closeness(value: float) -> str:
    """The colour of a C/E (or of a ratio of two runs): how far it is from 1."""
    d = abs(math.log(value))
    return "green" if d <= math.log(1.1) else "yellow" if d <= math.log(1.5) else "red"


def _percent(value: float) -> str:
    return f"{value:.0%}"


def report(rows: list, calculations: list, worst: int = 15, nis_name: str = "NIS", color: bool = False,
           unicode: bool = True) -> str:
    """nis_name: the name of the run in row["nis"], i.e. "NIS (JEFF-4.0)"."""
    style = _Style(color, unicode)
    uncertainties = case_uncertainties(rows)
    first = _short(nis_name)
    lines = []

    for (benchmark, quantity), cases in case_ratios(rows, calculations).items():
        codes = ["nis"] + [c for c in calculations if any(c in r for r in cases.values())]
        names = [first if c == "nis" else _short(c) for c in codes]
        # a letter for each run, so that the tables stay narrow
        letters = [chr(ord("A") + i) for i in range(len(codes))]
        points = sum(len(r["nis"]) for r in cases.values())
        lines.append(style(f"{benchmark} · {quantity.replace('_', ' ')} · {len(cases)} cases · {points} points",
                           "bold", "cyan"))
        lines.append("")
        lines += _table(style, [["", "Run"]], [[style(l, "bold"), n] for l, n in zip(letters, names)], ["^", "<"],
                        "The runs")
        lines.append("")

        # 1. how many cases and points are within 10, 20 and 50% of the experiment
        table_rows = []
        for tolerance in (0.1, 0.2, 0.5):
            values = []
            for code in codes:
                with_code = [r[code] for r in cases.values() if code in r]
                values.append((fraction_within([geometric_mean(v) for v in with_code], tolerance),
                               fraction_within([x for v in with_code for x in v], tolerance)))
            best = max(v[0] for v in values)
            row = [_percent(tolerance)]
            for c, p in values:
                cases_text = style(_percent(c), "bold", "green") if c == best and len(codes) > 1 else _percent(c)
                row.append(f"{_pad(cases_text, 5, '>')}   {_pad(_percent(p), 6, '>')}")
            table_rows.append(row)
        headers = [["C/E within"] + [f"{l}: {n}" for l, n in zip(letters, names)],
                   [""] + [f"{'cases':>5}   {'points':>6}"] * len(codes)]
        lines += _table(style, headers, table_rows, ["<"] + ["^"] * len(codes),
                        "Cases and points with C/E within the limit (the C/E of a case is the geometric mean of its points)")
        lines.append("")

        # 2. how close the other runs and codes are to the first run
        if len(codes) > 1:
            table_rows = []
            for code, letter, name in zip(codes[1:], letters[1:], names[1:]):
                ratio = [geometric_mean(r["nis"]) / geometric_mean(r[code]) for r in cases.values() if code in r]
                row = [f"A / {letter}", f"{first} / {name}"]
                for tolerance in (0.1, 0.2, 0.5):
                    row.append(_percent(fraction_within(ratio, tolerance)))
                table_rows.append(row)
            lines += _table(style, [["", "Runs", "within 10%", "within 20%", "within 50%"]], table_rows,
                            ["^", "<", ">", ">", ">"], "Cases where the C/E of A and of the other runs agree")
            lines.append("")

        # 3. the cases with the largest differences
        def uncertainty_cell(case):
            u = uncertainties.get((benchmark, quantity, case), math.nan)
            if math.isnan(u):
                return "-"
            return style(f"±{u:.0%}", "yellow" if u > UNCERTAIN else "")

        def case_cell(case):
            u = uncertainties.get((benchmark, quantity, case), math.nan)
            if math.isnan(u) or u <= UNCERTAIN:
                return case
            return style(case, "yellow") if color else f"{case} *"

        if worst:
            for code, letter, name in zip(codes[1:], letters[1:], names[1:]):
                pairs = [(case, geometric_mean(r["nis"]), geometric_mean(r[code]))
                         for case, r in cases.items() if code in r]
                pairs.sort(key=lambda p: -abs(math.log(p[1] / p[2])))
                table_rows = [[str(i), case_cell(case), uncertainty_cell(case), style(f"{n:.3f}", _closeness(n)),
                               style(f"{c:.3f}", _closeness(c)), style(f"{n / c:.3f}", _closeness(n / c))]
                              for i, (case, n, c) in enumerate(pairs[:worst], 1)]
                lines += _table(style, [["#", "Case", "Exp. unc.", "C/E A", f"C/E {letter}", f"A / {letter}"]],
                                table_rows, [">", "<", ">", ">", ">", ">"],
                                f"Cases where A ({first}) and {letter} ({name}) differ most")
                lines.append("")
            if len(codes) == 1:
                far = sorted(((case, geometric_mean(r["nis"])) for case, r in cases.items()),
                             key=lambda p: -abs(math.log(p[1])))
                table_rows = [[str(i), case_cell(case), uncertainty_cell(case), style(f"{n:.3f}", _closeness(n))]
                              for i, (case, n) in enumerate(far[:worst], 1)]
                lines += _table(style, [["#", "Case", "Exp. unc.", "C/E A"]], table_rows,
                                [">", "<", ">", ">"], f"Cases where A ({first}) is furthest from the experiment")
                lines.append("")

        if color:
            legend = (f"C/E and ratios: {style('within 10%', 'green')}, {style('within 50%', 'yellow')}, "
                      f"{style('more', 'red')}. The best run of a row is {style('bold', 'bold', 'green')}.\n"
                      f"Exp. unc.: the median relative uncertainty of the measured values of the case. "
                      f"{style(f'Cases over {UNCERTAIN:.0%}', 'yellow')} can not test the codes well.")
        else:
            legend = (f"Exp. unc.: the median relative uncertainty of the measured values of the case. "
                      f"Cases over {UNCERTAIN:.0%} (*) can not test the codes well.")
        lines.append(legend)
        lines.append("")
    return "\n".join(lines)


def strip_colors(text: str) -> str:
    return _ANSI.sub("", text)


def _short(code: str) -> str:
    """NIS, NIS-TENDL-2025 for the run 'NIS (TENDL-2025)', FISPACT-II for 'FISPACT-II (TENDL-2017)'"""
    if code == "nis":
        return "NIS"
    if code.startswith("NIS (") and code.endswith(")"):
        return "NIS-" + code[5:-1]
    return code.split(" ")[0]
