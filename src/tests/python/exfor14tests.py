"""
Cross sections at given energies (energies_ev in the benchmark format), and the EXFOR 14 MeV benchmark
(benchmarks/exfor14):

    - the format: a cross section with energies_ev in place of steps, a case with no steps
    - the results file: the columns energy_ev and source, two values at the same energy of one data set
    - the converter: the rules to select the EXFOR data sets
    - the engine: the cross section at an energy is the value of the group of the reaction data that has it
      (needs data/jeff4.reac and data/jeff311.dec, processing/process.py bootstrap)
    - the benchmark: the case files, and a regression test that runs all the cases with pynis and compares with
      benchmarks/exfor14/results/nis_reference.csv. If a change of the results is intended, make the reference again:

        uv run python scripts/bench.py run exfor14 --decay-data data/jeff311.dec \\
            --reaction-data data/jeff4.reac --jobs 4 --output benchmarks/exfor14/results/nis_reference.csv
"""

import copy
import json
import math
import os
import tempfile
import unittest
from pathlib import Path

from .testdata import BENCHMARKS_DIR, JEFF311_DECAY, JEFF4_REACTION, has_data

from nisbench import FormatError, format, load_benchmark, results
from nisbench.convert import exfor14


def _energy_case() -> dict:
    return {
        "format": "nis-benchmark-case", "version": 1, "benchmark": "test", "name": "Fe56/Mn56",
        "sample": {"atoms": {"Fe56": 1.0}},
        "spectrum": None,
        "measurements": [{"quantity": "cross_section", "unit": "mb", "target": "Fe56", "product": "Mn56",
                          "all_states": True, "source": "EXFOR 41240-012-0, (N,P), A.A.Filatenkov 1999",
                          "energies_ev": [14.1e6, 14.1e6, 14.7e6], "values": [110.0, 112.0, 105.0],
                          "uncertainties": [3.0, 3.0, 4.0]}],
    }


class EnergyFormatTests(unittest.TestCase):
    def test_case(self) -> None:
        case = format.parse_case(_energy_case())
        self.assertEqual(case.steps, [])
        self.assertFalse(case.has_deck)
        m, = case.measurements
        self.assertEqual(m.energies_ev, [14.1e6, 14.1e6, 14.7e6])
        self.assertEqual(m.steps, [])
        self.assertEqual(m.source, "EXFOR 41240-012-0, (N,P), A.A.Filatenkov 1999")

    def test_write_and_read(self) -> None:
        case = format.parse_case(_energy_case())
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Fe56" / "Mn56.json"
            format.write_case(case, path)
            data = json.loads(path.read_text())
            self.assertNotIn("steps", data)
            self.assertNotIn("steps", data["measurements"][0])
            again = format.load_case(path)
        self.assertEqual(again.measurements[0].energies_ev, case.measurements[0].energies_ev)
        self.assertEqual(again.measurements[0].source, case.measurements[0].source)

    def test_bad_cases(self) -> None:
        def changed(change):
            data = copy.deepcopy(_energy_case())
            change(data, data["measurements"][0])
            return data

        bad = {
            "steps and energies": lambda d, m: m.update(steps=[0]),
            "energies of a decay heat": lambda d, m: m.update(quantity="decay_heat", unit="W"),
            "energies go down": lambda d, m: m.update(energies_ev=[14.7e6, 14.1e6, 14.1e6]),
            "energy of 0": lambda d, m: m.update(energies_ev=[0.0, 14.1e6, 14.7e6]),
            "a value too few": lambda d, m: m.update(values=[1.0, 2.0]),
            "steps with no steps in the case": lambda d, m: (m.pop("energies_ev"), m.update(steps=[0])),
            "source not a text": lambda d, m: m.update(source=12),
            "empty steps": lambda d, m: d.update(steps={"dt_s": [], "flux_n_m2_s": []}),
        }
        for name, change in bad.items():
            with self.subTest(name), self.assertRaises(FormatError):
                format.parse_case(changed(change))


class EnergyResultsTests(unittest.TestCase):
    def _results(self, values: dict):
        case = format.parse_case(_energy_case())
        from nisbench.engines import CaseResults
        at_energy = {("Fe56", "Mn56", True, e): v * 1e-27 for e, v in values.items()}
        return case, CaseResults([], at_energy)

    def test_rows(self) -> None:
        case, step_results = self._results({14.1e6: 111.0, 14.7e6: 104.0})
        rows = results.case_rows(case, step_results)
        self.assertEqual([r["energy_ev"] for r in rows], ["1.41e+07", "1.41e+07", "1.47e+07"])
        self.assertEqual([float(r["nis"]) for r in rows], [111.0, 111.0, 104.0])
        self.assertEqual({r["time_s"] for r in rows}, {""})

    def test_write_and_read_runs(self) -> None:
        case, first = self._results({14.1e6: 111.0, 14.7e6: 104.0})
        _, second = self._results({14.1e6: 100.0, 14.7e6: 90.0})
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "a.csv", Path(tmp) / "b.csv"
            results.write_results(a, results.case_rows(case, first), [])
            column = results.nis_column("TENDL")
            results.write_results(b, results.case_rows(case, second, column), [], column)
            self.assertTrue(a.read_text().startswith(
                "benchmark,case,quantity,unit,time_s,energy_ev,source,experiment,uncertainty,nis\n"))
            rows, calculations = results.read_results([a, b])
        # the two values at 14.1 MeV of the data set are two rows
        self.assertEqual(len(rows), 3)
        self.assertEqual(calculations, ["NIS (TENDL)"])
        self.assertEqual(sorted((r["energy_ev"], r["experiment"], r["nis"], r["NIS (TENDL)"]) for r in rows),
                         [(14.1e6, 110.0, 111.0, 100.0), (14.1e6, 112.0, 111.0, 100.0), (14.7e6, 105.0, 104.0, 90.0)])

    def test_files_without_energies(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.csv"
            path.write_text("benchmark,case,quantity,unit,time_s,experiment,uncertainty,nis\n"
                            "fns,Co/x,decay_heat,uW/g,34.98,0.02729,0.001621,2.832815e-02\n")
            rows, _ = results.read_results([path])
        self.assertTrue(math.isnan(rows[0]["energy_ev"]))
        self.assertEqual(rows[0]["source"], "")


def _header(code: str) -> dict:
    return {"Reaction code": code, "Entry id": "41240-012-0 (entry-subentry-pointer)",
            "First author": "A.A.Filatenkov", "Year": "1999"}


# ZAI -> stable, as read_decay_data gives it
_DECAY = {250560: False, 250550: True, 410920: False, 410921: False, 270580: False, 270581: False, 20040: True,
          10010: True, 220500: True, 260560: True, 491160: False, 491161: False, 491162: False}
_ROWS = [(13.4, 0.0, 0.1, 0.01), (14.1, 0.05, 0.11, 0.004), (14.7, 0.05, 0.10, 0.0), (15.0, 0.05, 0.09, 0.005),
         (15.1, 0.0, 0.08, 0.004)]


class Exfor14ConverterTests(unittest.TestCase):
    def select(self, code: str):
        return exfor14.select(_header(code), _ROWS, _DECAY)

    def test_points(self) -> None:
        target, product, all_states, process, points = self.select("(26-FE-56(N,P)25-MN-56,,SIG)")
        self.assertEqual((target, product, all_states, process), ("Fe56", "Mn56", True, "P"))
        # 13.5 to 15 MeV, with an uncertainty, in mb
        self.assertEqual([(e, round(v, 6), round(u, 6)) for e, v, u in points],
                         [(14.1e6, 110.0, 4.0), (15.0e6, 90.0, 5.0)])

    def test_states(self) -> None:
        self.assertEqual(self.select("(41-NB-93(N,2N)41-NB-92-M,,SIG)")[1:3], ("Nb92m", False))
        self.assertEqual(self.select("(41-NB-93(N,2N)41-NB-92-M1,,SIG)")[1:3], ("Nb92m", False))
        self.assertEqual(self.select("(27-CO-59(N,2N)27-CO-58-G,,SIG)")[1:3], ("Co58", False))
        self.assertEqual(self.select("(27-CO-59(N,2N)27-CO-58-T,,SIG)")[1:3], ("Co58", True))
        self.assertEqual(self.select("(50-SN-116(N,P)49-IN-116-M2,,SIG)")[1:3], ("In116n", False))

    def test_not_used(self) -> None:
        for code in ("(26-FE-56(N,P)25-MN-56,,SIG,,SPA)",       # a spectrum average
                     "(26-FE-56(N,P)25-MN-56,PAR,SIG)",         # a partial cross section
                     "(26-FE-56(N,P)25-MN-56,,SIG)/(13-AL-27(N,A)11-NA-24,,SIG)",  # a ratio
                     "(26-FE-56(N,TOT),,SIG)",                  # no product
                     "(50-SN-116(N,P)49-IN-116-M1+M2,,SIG)",    # the sum of two isomers
                     "(42-MO-92(N,INL)42-MO-92-L,,SIG)",        # a level
                     "(26-FE-56(N,EL)26-FE-56,,SIG)",           # the target itself
                     "(26-FE-56(N,D)25-MN-55,,SIG)",            # a stable product of one reaction
                     "(26-FE-56(N,P)25-MN-57,,SIG)"):           # not in the decay data
            with self.subTest(code):
                self.assertIsNone(self.select(code))

    def test_production_of_a_stable_product(self) -> None:
        self.assertEqual(self.select("(26-FE-56(N,X)2-HE-4,,SIG)")[1:4], ("He4", True, "X"))
        self.assertEqual(self.select("(26-FE-56(N,X)25-MN-55,,SIG)")[1:4], ("Mn55", True, "X"))

    def test_read_table_and_source(self) -> None:
        text = ("# Header: \n#   Title: 26-FE-56(N,P)  Cross section \n# EXFOR BIB:\n"
                "#   Entry id: 41240-012-0 (entry-subentry-pointer) \n"
                "#   Reaction code: (26-FE-56(N,P)25-MN-56,,SIG) \n#   First author: A.A.Filatenkov \n"
                "#   Year:  1999 \n#\n#   E_in_lab(MeV)    dE_in_lab(MeV)           data(B)          ddata(B)\n"
                "        1.4100E+01        5.0000E-02        1.1000E-01        4.0000E-03\n")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.txt"
            path.write_text(text)
            header, rows = exfor14.read_table(path)
        self.assertEqual(rows, [(14.1, 0.05, 0.11, 0.004)])
        self.assertEqual(exfor14.source_text(header, "P"), "EXFOR 41240-012-0, (N,P), A.A.Filatenkov 1999")

    def test_make_case(self) -> None:
        case = exfor14.make_case("Fe56", "Mn56", True, [("EXFOR b", [(14.7e6, 1.0, 0.1), (14.1e6, 2.0, 0.1)]),
                                                        ("EXFOR a", [(14.5e6, 3.0, 0.1)])])
        self.assertEqual(case.name, "Fe56/Mn56")
        self.assertEqual([m.source for m in case.measurements], ["EXFOR a", "EXFOR b"])
        self.assertEqual(case.measurements[1].energies_ev, [14.1e6, 14.7e6])
        self.assertEqual(case.measurements[1].values, [2.0, 1.0])
        self.assertEqual(exfor14.make_case("Co59", "Co58", False, [("x", [(14e6, 1.0, 0.1)])]).name, "Co59/Co58g")


EXFOR14_DIR = BENCHMARKS_DIR / "exfor14"
EXFOR14_REFERENCE = EXFOR14_DIR / "results" / "nis_reference.csv"
# the results on another machine or compiler are not the same to the last digit
REFERENCE_REL_TOL = 1e-6


class Exfor14CasesTests(unittest.TestCase):
    def test_all_cases_load(self) -> None:
        benchmark = load_benchmark(EXFOR14_DIR)
        self.assertEqual(len(benchmark.cases), 1162)
        targets = set(exfor14.fns_targets())
        for case in benchmark.cases:
            with self.subTest(case=case.name):
                self.assertFalse(case.has_deck)
                self.assertEqual(case.steps, [])
                target, = case.atoms
                self.assertIn(target, targets)
                for m in case.measurements:
                    self.assertEqual((m.quantity, m.unit, m.target), ("cross_section", "mb", target))
                    self.assertTrue(m.source.startswith("EXFOR "))
                    self.assertTrue(all(13.5e6 <= e <= 15e6 for e in m.energies_ev))
                    self.assertTrue(all(u > 0 for u in m.uncertainties))

    def test_reference_has_all_cases(self) -> None:
        rows, _ = results.read_results([EXFOR14_REFERENCE])
        self.assertEqual({r["case"] for r in rows}, {c.name for c in load_benchmark(EXFOR14_DIR).cases})
        self.assertEqual(len(rows), sum(len(m.values) for c in load_benchmark(EXFOR14_DIR).cases
                                        for m in c.measurements))


@unittest.skipUnless(has_data(JEFF4_REACTION) and has_data(JEFF311_DECAY),
                     "needs data/jeff4.reac and data/jeff311.dec (processing/process.py bootstrap)")
class EnergyEngineTests(unittest.TestCase):
    def test_value_of_the_group(self) -> None:
        # Fe56 makes Mn56 only with (n,p), MT103: the cross section at E is the value of its group that has E
        import bisect
        import pynis as nis
        from nisbench import engines
        reaction_data = nis.core.reaction_data_from_file(str(engines.reaction_cache(JEFF4_REACTION)))
        xs, = reaction_data.find(260560, 103)
        energies = list(xs.energies)
        case = format.parse_case(_energy_case())
        step_results = next(iter(engines.run([case], engines.NuclearData(decay=JEFF311_DECAY,
                                                                          reaction=engines.reaction_cache(JEFF4_REACTION)))))
        for e in (14.1e6, 14.7e6):
            expected = xs.values[bisect.bisect_right(energies, e) - 1] * 1e-24
            self.assertGreater(expected, 0.0)
            self.assertAlmostEqual(step_results.at_energy[("Fe56", "Mn56", True, e)] / expected, 1.0, places=6)


@unittest.skipUnless(has_data(JEFF4_REACTION) and has_data(JEFF311_DECAY),
                     "needs data/jeff4.reac and data/jeff311.dec (processing/process.py bootstrap)")
class Exfor14RegressionTests(unittest.TestCase):
    def test_same_as_reference(self) -> None:
        from nisbench import engines
        cases = load_benchmark(EXFOR14_DIR).cases
        data = engines.NuclearData(decay=JEFF311_DECAY, reaction=engines.reaction_cache(JEFF4_REACTION))
        rows = []
        for case, step_results in zip(cases, engines.run(cases, data, "pynis", jobs=os.cpu_count() or 1)):
            rows += results.case_rows(case, step_results)
        reference, _ = results.read_results([EXFOR14_REFERENCE])
        self.assertEqual(len(rows), len(reference))
        # the rows are in the same order as in the reference file
        self.assertEqual([(r["case"], r["source"], float(r["energy_ev"])) for r in rows],
                         [(e["case"], e["source"], e["energy_ev"]) for e in reference])
        changed = [f"{r['case']} {r['source']} {r['energy_ev']}: {float(r['nis']):.6e}, reference {e['nis']:.6e}"
                   for r, e in zip(rows, reference)
                   if not math.isclose(float(r["nis"]), e["nis"], rel_tol=REFERENCE_REL_TOL, abs_tol=1e-30)]
        self.assertFalse(changed, f"{len(changed)} values are not the same as the reference:\n" + "\n".join(changed[:20]))
