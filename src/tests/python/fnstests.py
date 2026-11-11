"""
FNS decay heat regression test: runs all the FNS cases in benchmarks/fns with pynis and compares with

    - benchmarks/fns/results/nis_reference.csv, the NIS results when the reference was made. Any change to
      the results fails the test. If the change is intended, make the reference again (see the README):

        uv run python scripts/bench.py run fns --decay-data data/jeff311.dec \\
            --reaction-data data/jeff4.reac --output benchmarks/fns/results/nis_reference.csv

    - the experiment: NIS must stay as close to the measured heat as it is now.

Needs the nuclear data from git lfs. About 30 s on 4 cores.
"""

import csv
import math
import os
import unittest

import pytest

from .testdata import FNS_DIR, FNS_REFERENCE, JEFF311_DECAY, JEFF4_REACTION, has_data

from nisbench import load_benchmark
from nisbench.engines import NuclearData, reaction_cache, run
from nisbench.results import case_rows
from nisbench.stats import fraction_within, geometric_mean

# the results on another machine or compiler are not the same to the last digit
REFERENCE_REL_TOL = 1e-4

# the fraction of cases with a mean NIS C/E within 20 % and 50 % (the last numbers, minus a margin)
MIN_CASES_WITHIN_20 = 0.58
MIN_CASES_WITHIN_50 = 0.78


@pytest.mark.slow
@unittest.skipUnless(has_data(JEFF311_DECAY) and has_data(JEFF4_REACTION),
                     "needs data/jeff311.dec and data/jeff4.reac from git lfs")
class FnsRegressionTests(unittest.TestCase):
    results = None

    @classmethod
    def setUpClass(cls) -> None:
        cases = load_benchmark(FNS_DIR).cases
        data = NuclearData(decay=JEFF311_DECAY, reaction=reaction_cache(JEFF4_REACTION))
        # case -> [(cooling time, experiment, nis), ...]
        cls.results = {}
        for case, step_results in zip(cases, run(cases, data, "pynis", jobs=os.cpu_count() or 1)):
            cls.results[case.name] = [(float(r["time_s"]), float(r["experiment"]), float(r["nis"]))
                                      for r in case_rows(case, step_results)]

    def test_same_as_reference(self) -> None:
        with open(FNS_REFERENCE) as f:
            reference = {}
            for r in csv.DictReader(f):
                reference.setdefault(r["case"], []).append(float(r["nis"]))

        self.assertEqual(set(reference), set(self.results), "the cases are not the same as in the reference")
        changed = []
        for name, rows in self.results.items():
            self.assertEqual(len(rows), len(reference[name]), f"{name}: not the same number of times")
            for (t, _, heat), expected in zip(rows, reference[name]):
                if not math.isclose(heat, expected, rel_tol=REFERENCE_REL_TOL, abs_tol=1e-30):
                    changed.append(f"{name} t={t:.6g} s: {heat:.6e} uW/g, reference {expected:.6e}")
        self.assertFalse(changed, f"{len(changed)} values are not the same as the reference, the first ones:\n"
                         + "\n".join(changed[:20]))

    def test_close_to_experiment(self) -> None:
        case_ce = []
        for rows in self.results.values():
            ce = [heat / exp for _, exp, heat in rows if exp > 0 and heat > 0]
            if ce:
                case_ce.append(geometric_mean(ce))
        self.assertEqual(len(case_ce), len(self.results))
        self.assertGreaterEqual(fraction_within(case_ce, 0.2), MIN_CASES_WITHIN_20)
        self.assertGreaterEqual(fraction_within(case_ce, 0.5), MIN_CASES_WITHIN_50)


__all__ = ["FnsRegressionTests"]
