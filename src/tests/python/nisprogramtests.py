"""
Tests of the nis program (from the CMake build, or NIS_BIN). Some FNS decks are run with the program and
must give the same decay heat as benchmarks/fns/results/nis_reference.csv (made with pynis). This tests the parts
that only the program has: the input decks, the steps with a flux and the cooling steps, and the output.

Needs the nuclear data from git lfs.
"""

import csv
import math
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

import pytest

from .testdata import FNS_DIR, FNS_REFERENCE, JEFF311_DECAY, JEFF4_REACTION, NIS_BIN, has_data

from nisbench import load_benchmark
from nisbench.engines import NuclearData, parse_nis_output, reaction_cache, run
from nisbench.results import case_rows

CASES = ["Co/2000exp_5min", "Fe/1996exp_7hour", "Ni/2000exp_5min"]


def _nis(*arguments) -> subprocess.CompletedProcess:
    return subprocess.run([str(NIS_BIN), *map(str, arguments)], capture_output=True, text=True)


@pytest.mark.slow
@unittest.skipUnless(NIS_BIN.exists(), f"needs the nis program ({NIS_BIN}), build it with CMake or give NIS_BIN")
@unittest.skipUnless(has_data(JEFF311_DECAY) and has_data(JEFF4_REACTION),
                     "needs data/jeff311.dec and data/jeff4.reac from git lfs")
class NisProgramTests(unittest.TestCase):
    def test_decks_same_as_reference(self) -> None:
        # nis_cache next to nis makes the cache, if it is not there yet
        cache = reaction_cache(JEFF4_REACTION, NIS_BIN.parent / "nis_cache")
        with open(FNS_REFERENCE) as f:
            reference = list(csv.DictReader(f))

        cases = load_benchmark(FNS_DIR, CASES).cases
        self.assertEqual([c.name for c in cases], CASES)
        results = run(cases, NuclearData(decay=JEFF311_DECAY, reaction=cache), engine="nis", nis=NIS_BIN, jobs=3)
        for case, step_results in zip(cases, results):
            expected = [float(r["nis"]) for r in reference if r["case"] == case.name]
            heats = [float(r["nis"]) for r in case_rows(case, step_results)]
            self.assertEqual(len(heats), len(expected), case.name)
            for heat, value in zip(heats, expected):
                self.assertTrue(math.isclose(heat, value, rel_tol=1e-4),
                                f"{case.name}: {heat:.6e} uW/g, reference {value:.6e}")

    def test_decay_only_steps(self) -> None:
        # Co60: 2 steps of one half life with no spectrum, so no reaction data is needed
        half_life = 1.6634e8
        with tempfile.TemporaryDirectory() as directory:
            deck = Path(directory) / "co60.inp"
            deck.write_text(f"dt {half_life}\ndt {half_life} 0\n\nnuclides\n    Co60 1e+20\n")
            result = _nis("-d", JEFF311_DECAY, "-i", deck)
        self.assertEqual(result.returncode, 0, result.stderr)

        # one block for each step
        self.assertEqual(re.findall(r"^step (\d+)$", result.stdout, re.MULTILINE), ["1", "2"])
        co60 = [float(v) for v in re.findall(r"^\s*Co60\s+([0-9.eE+-]+)", result.stdout, re.MULTILINE)]
        self.assertEqual(len(co60), 2)
        self.assertTrue(math.isclose(co60[0], 0.5e20, rel_tol=1e-3), co60)
        self.assertTrue(math.isclose(co60[1], 0.25e20, rel_tol=1e-3), co60)
        heats = [step["decay_heat"] for step in parse_nis_output(result.stdout)]
        self.assertEqual(len(heats), 2)
        self.assertTrue(math.isclose(heats[1] / heats[0], 0.5, rel_tol=1e-3))

    def test_flux_without_spectrum_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            deck = Path(directory) / "bad.inp"
            deck.write_text("dt 10 1e14\n\nnuclides\n    Co59 1e+20\n")
            result = _nis("-d", JEFF311_DECAY, "-i", deck)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("spectrum", result.stderr)

    def test_spectrum_without_reaction_data_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            deck = Path(directory) / "bad.inp"
            deck.write_text("dt 10\n\nnuclides\n    Co59 1e+20\n\nspectrum\n    1e6 2e6\n    1e10\n")
            result = _nis("-d", JEFF311_DECAY, "-i", deck)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--reaction-data", result.stderr)


__all__ = ["NisProgramTests"]
