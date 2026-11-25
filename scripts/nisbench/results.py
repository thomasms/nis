"""
The standard results file: one CSV for all benchmarks, one row for each measured value.

    benchmark, case, quantity, unit, time_s, experiment, uncertainty, nis, <calculation>, ...

time_s is the cooling time (after the end of the last irradiation step). The values are in the unit
of the measurement. A file with cross sections at given energies has two more columns after time_s:

    energy_ev   the incident neutron energy of the value (time_s is then empty)
    source      where the measured value comes from, e.g. an EXFOR subentry

Two rows can have the same energy and source (a data set with two values at one energy): they are told apart by
their order in the file. The calculations (for example FISPACT-II) come from the case file, a cell is
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
# only in a file with cross sections at given energies
ENERGY_COLUMNS = ["energy_ev", "source"]


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


def nis_value(case: Case, quantity: str, unit: str, step_result: dict, cooling_time: float,
              reaction: tuple = None, value: float = None) -> float:
    """The NIS value (in the unit of the measurement) from the result of a step (W, Bq, and cm2 for each
    (target, product, all_states) of a cross section). reaction: (target, product, all_states), only for a cross
    section. value: the SI value (cm2 for a cross section) in place of the step result."""
    factor, normalisation, times_time = units_for(quantity)[unit]
    if value is None:
        value = step_result[quantity][reaction] if quantity == "cross_section" else step_result[quantity]
    value *= factor
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
    """step_results: the results of the engine for the case (with at_energy for the cross sections at energies)."""
    rows = []
    for m in case.measurements:
        if m.energies_ev is not None:
            rows += _energy_rows(case, m, step_results.at_energy, column)
            continue
        times = case.cooling_times(m.steps)
        for i, (step, t) in enumerate(zip(m.steps, times)):
            row = {"benchmark": case.benchmark, "case": case.name, "quantity": m.quantity, "unit": m.unit,
                   "time_s": f"{t:.6g}", "experiment": f"{m.values[i]:.6g}", "uncertainty": f"{m.uncertainties[i]:.6g}",
                   column: f"{nis_value(case, m.quantity, m.unit, step_results[step], t, (m.target, m.product, m.all_states)):.6e}"}
            for name, values in m.calculations.items():
                row[name] = "" if values[i] is None else f"{values[i]:.6e}"
            rows.append(row)
    return rows


def _energy_rows(case: Case, m, at_energy: dict, column: str) -> list:
    rows = []
    for i, energy in enumerate(m.energies_ev):
        value = at_energy[(m.target, m.product, m.all_states, energy)]
        row = {"benchmark": case.benchmark, "case": case.name, "quantity": m.quantity, "unit": m.unit, "time_s": "",
               "energy_ev": f"{energy:.6g}", "source": m.source or "", "experiment": f"{m.values[i]:.6g}",
               "uncertainty": f"{m.uncertainties[i]:.6g}",
               column: f"{nis_value(case, m.quantity, m.unit, None, 0.0, value=value):.6e}"}
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
        base = BASE_COLUMNS
        if any("energy_ev" in r for r in rows):
            base = BASE_COLUMNS[:5] + ENERGY_COLUMNS + BASE_COLUMNS[5:]
        writer = csv.DictWriter(f, base + [column] + calculations, restval="")
        writer.writeheader()
        writer.writerows(rows)


def _float(text: str) -> float:
    return math.nan if text in ("", None) else float(text)


def read_runs(paths: list) -> tuple:
    """(rows, calculation names, name of the first run) of one or more results files.

    The NIS values of the first run are in row["nis"]. Each other run (another label) is a calculation, with
    its name (i.e. "NIS (TENDL-2025)"), in front of the calculations of the case files. Files of the same run (the
    same label, or no label) are joined, and a later file replaces the same (benchmark, case, quantity, time, energy,
    source) of an earlier one. The values are floats, nan where empty."""
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
            calculations += [c for c in fieldnames if c not in BASE_COLUMNS and c not in ENERGY_COLUMNS
                             and c != nis_columns[0] and c not in calculations]
            seen = {}
            for r in reader:
                point = (r["benchmark"], r["case"], r["quantity"], r["time_s"], r.get("energy_ev", ""),
                         r.get("source", ""))
                # the n-th row of the same point in this file
                seen[point] = seen.get(point, -1) + 1
                row = rows.setdefault(point + (seen[point],), {})
                row.update({c: r[c] for c in BASE_COLUMNS})
                row["energy_ev"], row["source"] = r.get("energy_ev", ""), r.get("source", "")
                row[key] = r[nis_columns[0]]
                row.update({c: r[c] for c in fieldnames if c in calculations})

    calculations = runs[1:] + calculations
    for row in rows.values():
        for column in ["time_s", "energy_ev", "experiment", "uncertainty", "nis"] + calculations:
            row[column] = _float(row.get(column))
    return list(rows.values()), calculations, runs[0] if runs else "NIS"


def read_results(paths: list) -> tuple:
    """(rows, calculation names) of one or more results files, see read_runs."""
    rows, calculations, _ = read_runs(paths)
    return rows, calculations
