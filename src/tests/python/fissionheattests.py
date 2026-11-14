"""
Fission decay heat benchmark (benchmarks/fission_decay_heat, from the IAEA CoNDERC fission archive):

    - the case files: all the cases load, each is a fission source with no deck, the reference has all of them
    - regression: runs all the cases with pynis and compares with benchmarks/fission_decay_heat/results/
      nis_reference.csv. Any change to the results fails the test. If the change is intended, make the
      reference again:

        uv run python scripts/bench.py run fission_decay_heat --decay-data data/jeff311.dec \\
            --fission-yields data/jeff311.nfy --output benchmarks/fission_decay_heat/results/nis_reference.csv

    - the experiment: NIS must stay as close to the measured heat as it is now.

The regression test needs data/jeff311.dec (git lfs) and data/jeff311.nfy (scripts/get_fission_yields.py).
"""

import math
import os
import unittest

from .testdata import BENCHMARKS_DIR, JEFF311_DECAY, JEFF311_NFY, has_data

from nisbench import load_benchmark
from nisbench.engines import NuclearData, run
from nisbench.results import case_rows, read_results
from nisbench.stats import fraction_within, geometric_mean

FISSION_DIR = BENCHMARKS_DIR / "fission_decay_heat"
FISSION_REFERENCE = FISSION_DIR / "results" / "nis_reference.csv"

# the results on another machine or compiler are not the same to the last digit
REFERENCE_REL_TOL = 1e-4

# quantity -> (fraction of cases with a mean C/E within 10 %, within 20 %), the last numbers minus a margin.
# The beta heat cases out of 20 % are MacMahon (U235, 5 irradiations): NIS is 19-28 % above them.
MIN_CASES_WITHIN = {
    "decay_heat": (0.95, 1.0),
    "beta_heat": (0.75, 0.85),
    "gamma_heat": (0.70, 0.90),
}


class FissionDecayHeatCasesTests(unittest.TestCase):
    def test_all_cases_load(self) -> None:
        benchmark = load_benchmark(FISSION_DIR)
        self.assertEqual(benchmark.category, "fission")
        self.assertEqual(len(benchmark.cases), 44)
        for case in benchmark.cases:
            with self.subTest(case=case.name):
                self.assertIsNotNone(case.fission_source)
                self.assertFalse(case.has_deck)
                self.assertFalse(case.deck_path.exists())
                self.assertEqual(case.atoms, {})
                self.assertTrue(all(flux == 0 for _, flux in case.steps))
                pulse = case.name.endswith("pulse") or "_pulse_series" in case.name
                if pulse:
                    self.assertEqual(case.fissions, 1.0)
                    self.assertIsNone(case.last_fission_rate)
                else:
                    # an irradiation step with a fission rate, then cooling
                    self.assertGreater(case.step_fission_rates()[0], 0)
                    self.assertTrue(all(r == 0 for r in case.step_fission_rates()[1:]))
                for m in case.measurements:
                    self.assertEqual(m.unit, "MeV/fission" if pulse or "Baumung" in case.name else "MeV/(fission/s)")
                    self.assertTrue(all(v > 0 for v in m.values))

    def test_reference_has_all_cases(self) -> None:
        rows, _ = read_results([FISSION_REFERENCE])
        cases = load_benchmark(FISSION_DIR).cases
        self.assertEqual({r["case"] for r in rows}, {c.name for c in cases})
        self.assertEqual(len(rows), sum(len(m.values) for c in cases for m in c.measurements))


@unittest.skipUnless(has_data(JEFF311_DECAY) and has_data(JEFF311_NFY),
                     "needs data/jeff311.dec (git lfs) and data/jeff311.nfy (scripts/get_fission_yields.py)")
class FissionDecayHeatRegressionTests(unittest.TestCase):
    rows = None

    @classmethod
    def setUpClass(cls) -> None:
        cases = load_benchmark(FISSION_DIR).cases
        data = NuclearData(decay=JEFF311_DECAY, fission_yields=JEFF311_NFY)
        cls.rows = []
        for case, step_results in zip(cases, run(cases, data, "pynis", jobs=os.cpu_count() or 1)):
            cls.rows += case_rows(case, step_results)

    def test_same_as_reference(self) -> None:
        reference, _ = read_results([FISSION_REFERENCE])
        # the reference gives the time as a number, a new row as text
        key = lambda r: (r["case"], r["quantity"], float(r["time_s"]))
        expected = {key(r): float(r["nis"]) for r in reference}
        self.assertEqual(set(expected), {key(r) for r in self.rows}, "the values are not the same as in the reference")
        changed = [f"{r['case']} {r['quantity']} t={r['time_s']} s: {float(r['nis']):.6e}, reference {expected[key(r)]:.6e}"
                   for r in self.rows if not math.isclose(float(r["nis"]), expected[key(r)], rel_tol=REFERENCE_REL_TOL)]
        self.assertFalse(changed, f"{len(changed)} values are not the same as the reference, the first ones:\n"
                         + "\n".join(changed[:20]))

    def test_close_to_experiment(self) -> None:
        for quantity, (within_10, within_20) in MIN_CASES_WITHIN.items():
            ce = {}
            for r in self.rows:
                if r["quantity"] == quantity:
                    ce.setdefault(r["case"], []).append(float(r["nis"]) / float(r["experiment"]))
            case_ce = [geometric_mean(v) for v in ce.values()]
            with self.subTest(quantity=quantity):
                self.assertGreater(len(case_ce), 10)
                self.assertGreaterEqual(fraction_within(case_ce, 0.1), within_10)
                self.assertGreaterEqual(fraction_within(case_ce, 0.2), within_20)


__all__ = ["FissionDecayHeatCasesTests", "FissionDecayHeatRegressionTests"]
