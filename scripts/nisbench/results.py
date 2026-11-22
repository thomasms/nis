"""
The standard results file: one CSV for all benchmarks, one row for each measured value.

    benchmark, case, quantity, unit, time_s, experiment, uncertainty, nis, <calculation>, ...

time_s is the cooling time (after the end of the last irradiation step). The values are in the unit
of the measurement. The calculations (for example FISPACT-II) come from the case file, a cell is
empty where a calculation has no value.

A run with a label (bench.py run --label TENDL-2025) has the column "nis (TENDL-2025)" in place of "nis".
Then results files of different runs can be read together, and each run is a column of its own.
"""

import csv
import math
from pathlib import Path

from .format import Case, units_for

BASE_COLUMNS = ["benchmark", "case", "quantity", "unit", "time_s", "experiment", "uncertainty"]
FIXED_COLUMNS = BASE_COLUMNS + ["nis"]


def nis_column(label: str = None) -> str:
    """The column of the NIS values: "nis", or "nis (<label>)" for a run with a label."""
    if not label:
        return "nis"
    if any(c in label for c in "(),\n"):
        raise ValueError(f"the label {label!r} must not have '(', ')', ',' or a new line")
    return f"nis ({label})"


def run_name(column: str) -> str:
    """The name of a run to show: NIS, or NIS (<label>)."""
    return "NIS" + column[3:]


def _is_nis_column(column: str) -> bool:
    return column == "nis" or (column.startswith("nis (") and column.endswith(")"))


def nis_value(case: Case, quantity: str, unit: str, step_result: dict, cooling_time: float) -> float:
    """The NIS value (in the unit of the measurement) from the result of a step (W and Bq)."""
    factor, normalisation, times_time = units_for(quantity)[unit]
    value = step_result[quantity] * factor
    if normalisation == "mass":
        value /= case.mass_kg
    elif normalisation == "fission":
        value /= case.fissions
    elif normalisation == "fission_rate":
        value /= case.last_fission_rate
    if times_time:
        value *= cooling_time
    return value


def case_rows(case: Case, step_results: list, column: str = "nis") -> list:
    rows = []
    for m in case.measurements:
        times = case.cooling_times(m.steps)
        for i, (step, t) in enumerate(zip(m.steps, times)):
            row = {"benchmark": case.benchmark, "case": case.name, "quantity": m.quantity, "unit": m.unit,
                   "time_s": f"{t:.6g}", "experiment": f"{m.values[i]:.6g}", "uncertainty": f"{m.uncertainties[i]:.6g}",
                   column: f"{nis_value(case, m.quantity, m.unit, step_results[step], t):.6e}"}
            for name, values in m.calculations.items():
                row[name] = "" if values[i] is None else f"{values[i]:.6e}"
            rows.append(row)
    return rows


def calculation_names(cases: list) -> list:
    names = []
    for case in cases:
        for m in case.measurements:
            names += [n for n in m.calculations if n not in names]
    return names


def write_results(path: Path, rows: list, calculations: list, column: str = "nis") -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, BASE_COLUMNS + [column] + calculations, restval="")
        writer.writeheader()
        writer.writerows(rows)


def _float(text: str) -> float:
    return math.nan if text in ("", None) else float(text)


def read_runs(paths: list) -> tuple:
    """(rows, calculation names, name of the first run) of one or more results files.

    The NIS values of the first run are in row["nis"]. Each other run (another label) is a calculation, with
    its name (i.e. "NIS (TENDL-2025)"), in front of the calculations of the case files. Files of the same run (the
    same label, or no label) are joined, and a later file replaces the same (benchmark, case, quantity, time) of an
    earlier one. The values are floats, nan where empty."""
    rows, runs, calculations = {}, [], []
    for path in paths:
        with open(path) as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames or []
            missing = [c for c in BASE_COLUMNS if c not in fieldnames]
            nis_columns = [c for c in fieldnames if _is_nis_column(c)]
            if missing or len(nis_columns) != 1:
                what = f"no column {', '.join(missing)}" if missing else "it must have one nis column"
                raise ValueError(f"{path}: not a results file, {what}")
            name = run_name(nis_columns[0])
            if name not in runs:
                runs.append(name)
            key = "nis" if name == runs[0] else name
            calculations += [c for c in fieldnames
                             if c not in BASE_COLUMNS and c != nis_columns[0] and c not in calculations]
            for r in reader:
                row = rows.setdefault((r["benchmark"], r["case"], r["quantity"], r["time_s"]), {})
                row.update({c: r[c] for c in BASE_COLUMNS})
                row[key] = r[nis_columns[0]]
                row.update({c: r[c] for c in fieldnames if c in calculations})

    calculations = runs[1:] + calculations
    for row in rows.values():
        for column in ["time_s", "experiment", "uncertainty", "nis"] + calculations:
            row[column] = _float(row.get(column))
    return list(rows.values()), calculations, runs[0] if runs else "NIS"


def read_results(paths: list) -> tuple:
    """(rows, calculation names) of one or more results files, see read_runs."""
    rows, calculations, _ = read_runs(paths)
    return rows, calculations
