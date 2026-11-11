"""
The standard results file: one CSV for all benchmarks, one row for each measured value.

    benchmark, case, quantity, unit, time_s, experiment, uncertainty, nis, <calculation>, ...

time_s is the cooling time (after the end of the last irradiation step). The values are in the unit
of the measurement. The calculations (for example FISPACT-II) come from the case file, a cell is
empty where a calculation has no value.
"""

import csv
import math
from pathlib import Path

from .format import Case, units_for

FIXED_COLUMNS = ["benchmark", "case", "quantity", "unit", "time_s", "experiment", "uncertainty", "nis"]


def nis_value(case: Case, quantity: str, unit: str, step_result: dict, cooling_time: float) -> float:
    """The NIS value (in the unit of the measurement) from the result of a step (W and Bq)."""
    factor, normalisation, times_time = units_for(quantity)[unit]
    value = step_result[quantity] * factor
    if normalisation == "mass":
        value /= case.mass_kg
    elif normalisation == "fission":
        value /= case.fissions
    if times_time:
        value *= cooling_time
    return value


def case_rows(case: Case, step_results: list) -> list:
    rows = []
    for m in case.measurements:
        times = case.cooling_times(m.steps)
        for i, (step, t) in enumerate(zip(m.steps, times)):
            row = {"benchmark": case.benchmark, "case": case.name, "quantity": m.quantity, "unit": m.unit,
                   "time_s": f"{t:.6g}", "experiment": f"{m.values[i]:.6g}", "uncertainty": f"{m.uncertainties[i]:.6g}",
                   "nis": f"{nis_value(case, m.quantity, m.unit, step_results[step], t):.6e}"}
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


def write_results(path: Path, rows: list, calculations: list) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, FIXED_COLUMNS + calculations, restval="")
        writer.writeheader()
        writer.writerows(rows)


def _float(text: str) -> float:
    return math.nan if text in ("", None) else float(text)


def read_results(paths: list) -> tuple:
    """(rows, calculation names) of one or more results files. A later file replaces the same
    (benchmark, case, quantity, time) of an earlier one. The values are floats, nan where empty."""
    rows, calculations = {}, []
    for path in paths:
        with open(path) as f:
            reader = csv.DictReader(f)
            missing = [c for c in FIXED_COLUMNS if c not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(f"{path}: not a results file, no column {', '.join(missing)}")
            calculations += [c for c in reader.fieldnames if c not in FIXED_COLUMNS and c not in calculations]
            for r in reader:
                row = dict(r)
                for key in ["time_s", "experiment", "uncertainty", "nis"] + calculations:
                    row[key] = _float(r.get(key))
                rows[(r["benchmark"], r["case"], r["quantity"], r["time_s"])] = row
    return list(rows.values()), calculations
