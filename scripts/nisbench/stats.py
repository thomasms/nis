"""
C/E statistics of a results file, for each benchmark and quantity.

For each case, the C/E is the geometric mean over the measured times. Only the times where NIS and all
the calculations of the case have a value are used, so the codes are compared on the same points.
"""

import collections
import math


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


def report(rows: list, calculations: list, worst: int = 25) -> str:
    lines = []
    for (benchmark, quantity), cases in case_ratios(rows, calculations).items():
        codes = ["nis"] + [c for c in calculations if any(c in r for r in cases.values())]
        points = sum(len(r["nis"]) for r in cases.values())
        lines.append(f"{benchmark}, {quantity}: {len(cases)} cases, {points} points")
        lines.append(f"  {'C/E within':>10s}" + "".join(f" {'cases ' + _short(c):>16s} {'points':>7s}" for c in codes))
        for tolerance in (0.1, 0.2, 0.5):
            text = f"  {tolerance:>10.0%}"
            for code in codes:
                with_code = [r[code] for r in cases.values() if code in r]
                means = [geometric_mean(v) for v in with_code]
                all_points = [x for v in with_code for x in v]
                text += f" {fraction_within(means, tolerance):>16.0%} {fraction_within(all_points, tolerance):>7.0%}"
            lines.append(text)

        for code in codes[1:]:
            pairs = [(case, geometric_mean(r["nis"]), geometric_mean(r[code])) for case, r in cases.items() if code in r]
            ratio = [n / c for _, n, c in pairs]
            lines.append(f"  NIS/{code} within 10%: {fraction_within(ratio, 0.1):.0%}, within 20%: "
                         f"{fraction_within(ratio, 0.2):.0%}")
            if worst:
                lines.append(f"  cases where NIS and {code} differ most:")
                for case, n, c in sorted(pairs, key=lambda p: -abs(math.log(p[1] / p[2])))[:worst]:
                    lines.append(f"    {case:24s} NIS C/E {n:7.3f}  {_short(code)} C/E {c:7.3f}  ratio {n / c:7.3f}")
        if len(codes) == 1 and worst:
            lines.append("  cases furthest from the experiment:")
            far = sorted(((case, geometric_mean(r["nis"])) for case, r in cases.items()), key=lambda p: -abs(math.log(p[1])))
            for case, n in far[:worst]:
                lines.append(f"    {case:24s} NIS C/E {n:7.3f}")
        lines.append("")
    return "\n".join(lines)


def _short(code: str) -> str:
    return "NIS" if code == "nis" else code.split(" ")[0]
