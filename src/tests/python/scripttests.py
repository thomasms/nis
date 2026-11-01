"""
Fast tests of the benchmark scripts in scripts/. They need no nuclear data and no pynis.
"""

import math
import os
import stat
import struct
import sys
import tempfile
import textwrap
import time
import unittest
from pathlib import Path

from .testdata import FNS_DECKS_DIR

import fns_benchmark
import fns_run_decks
import fns_stats


class InterpolationTests(unittest.TestCase):
    def test_cumulative(self) -> None:
        self.assertEqual(fns_benchmark.cumulative([1.0, 2.0, 3.0]), [1.0, 3.0, 6.0])
        self.assertEqual(fns_benchmark.cumulative([]), [])

    def test_interpolate_log_at_the_points(self) -> None:
        times = [1.0, 2.0, 4.0]
        values = [100.0, 10.0, 1.0]
        for t, v in zip(times, values):
            self.assertAlmostEqual(fns_benchmark.interpolate_log(times, values, t), v)

    def test_interpolate_log_between_the_points(self) -> None:
        # the log is linear in time, so the middle is the geometric mean
        self.assertAlmostEqual(fns_benchmark.interpolate_log([0.0, 2.0], [100.0, 1.0], 1.0), 10.0)

    def test_interpolate_log_outside(self) -> None:
        self.assertTrue(math.isnan(fns_benchmark.interpolate_log([1.0, 2.0], [1.0, 2.0], 0.5)))
        self.assertTrue(math.isnan(fns_benchmark.interpolate_log([1.0, 2.0], [1.0, 2.0], 3.0)))
        # no log of 0
        self.assertTrue(math.isnan(fns_benchmark.interpolate_log([1.0, 2.0], [0.0, 2.0], 1.5)))


class DeckTests(unittest.TestCase):
    deck = FNS_DECKS_DIR / "Co" / "2000exp_5min.inp"

    def test_load_deck(self) -> None:
        case = fns_benchmark.Case(material="Co", name="2000exp_5min")
        fns_benchmark.load_deck(case, self.deck)

        self.assertGreater(case.mass_kg, 0.0)
        self.assertGreater(case.irradiation_time, 0.0)
        self.assertIn("Co59", case.initial_atoms)
        self.assertEqual(len(case.flux_shape), 709)

        # the deck flux is in n/m2/s, the case flux in n/cm2/s
        dt_line = next(line for line in self.deck.read_text().splitlines() if line.startswith("dt "))
        self.assertAlmostEqual(case.flux_amplitude, float(dt_line.split()[2]) * 1e-4)

        # one cooling step to each measured time
        self.assertEqual(len(case.cooling_steps), len(case.exp_time))
        for t, step_end in zip(case.exp_time, fns_benchmark.cumulative(case.cooling_steps)):
            self.assertAlmostEqual(step_end / t, 1.0, places=6)
        self.assertEqual(len(case.exp_heat), len(case.exp_time))
        self.assertEqual(len(case.fispact_at_exp), len(case.exp_time))
        self.assertEqual(fns_run_decks.read_mass(self.deck), case.mass_kg)

    def test_load_cases_reads_all_decks(self) -> None:
        cases = fns_benchmark.load_cases(FNS_DECKS_DIR, [])
        self.assertEqual(len(cases), 132)
        cases = fns_benchmark.load_cases(FNS_DECKS_DIR, ["Co", "Fe"])
        self.assertEqual({c.material for c in cases}, {"Co", "Fe"})

    def test_load_cases_with_no_cases(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(SystemExit):
                fns_benchmark.load_cases(Path(directory), [])

    def test_decay_heats_from_nis_output(self) -> None:
        output = textwrap.dedent("""\
            step 1
                time 300
            decay_heat
                total 2.5e-03
                beta 1e-3
            step 2
            decay_heat
                total 1.25e-03
            """)
        self.assertEqual(fns_run_decks.decay_heats(output), [2.5e-3, 1.25e-3])


class StatsTests(unittest.TestCase):
    def test_fraction_within(self) -> None:
        # symmetric in log: 1.25 and 0.8 are both within 25 %
        self.assertEqual(fns_stats.fraction_within([1.25, 0.8, 2.0, 1.0], 0.25), 0.75)

    def test_geometric_mean(self) -> None:
        self.assertAlmostEqual(fns_stats.geometric_mean([0.5, 2.0]), 1.0)


def _write_cache(path: Path, version: int = fns_benchmark._CACHE_VERSION) -> None:
    path.write_bytes(fns_benchmark._CACHE_MAGIC + struct.pack("=I", version) + b"\0" * 8)


class ReactionCacheTests(unittest.TestCase):
    """reaction_cache makes the cache with a fake nis_cache program, which writes only a header."""

    def setUp(self) -> None:
        self._directory = tempfile.TemporaryDirectory()
        self.directory = Path(self._directory.name)
        self.endf = self.directory / "data.reac"
        self.endf.write_text("not really ENDF\n")
        self.cache = self.directory / "data.nisreac"
        self.calls = self.directory / "calls.txt"

        self.nis_cache = self.directory / "nis_cache"
        self.nis_cache.write_text(textwrap.dedent(f"""\
            #!{sys.executable}
            import struct, sys
            output = sys.argv[sys.argv.index("-o") + 1]
            with open(output, "wb") as f:
                f.write({fns_benchmark._CACHE_MAGIC!r} + struct.pack("=I", {fns_benchmark._CACHE_VERSION}) + bytes(8))
            with open({str(self.calls)!r}, "a") as f:
                f.write("called\\n")
            """))
        self.nis_cache.chmod(self.nis_cache.stat().st_mode | stat.S_IXUSR)

    def tearDown(self) -> None:
        self._directory.cleanup()

    def made(self) -> int:
        return len(self.calls.read_text().splitlines()) if self.calls.exists() else 0

    def test_cache_version(self) -> None:
        _write_cache(self.cache, 7)
        self.assertEqual(fns_benchmark._cache_version(self.cache), 7)
        self.assertIsNone(fns_benchmark._cache_version(self.endf))
        self.assertIsNone(fns_benchmark._cache_version(self.directory / "missing"))

    def test_a_cache_is_used_as_it_is(self) -> None:
        _write_cache(self.cache)
        self.assertEqual(fns_benchmark.reaction_cache(self.cache, self.nis_cache), self.cache)
        self.assertEqual(self.made(), 0)

    def test_made_if_missing(self) -> None:
        self.assertEqual(fns_benchmark.reaction_cache(self.endf, self.nis_cache), self.cache)
        self.assertEqual(self.made(), 1)
        self.assertEqual(fns_benchmark._cache_version(self.cache), fns_benchmark._CACHE_VERSION)
        self.assertFalse(self.cache.with_suffix(".nisreac.tmp").exists())

        # it is up to date now, so it is not made again
        fns_benchmark.reaction_cache(self.endf, self.nis_cache)
        self.assertEqual(self.made(), 1)

    def test_made_if_older_than_the_data(self) -> None:
        _write_cache(self.cache)
        old = time.time() - 3600
        os.utime(self.cache, (old, old))
        fns_benchmark.reaction_cache(self.endf, self.nis_cache)
        self.assertEqual(self.made(), 1)

    def test_made_if_another_version(self) -> None:
        _write_cache(self.cache, fns_benchmark._CACHE_VERSION + 1)
        fns_benchmark.reaction_cache(self.endf, self.nis_cache)
        self.assertEqual(self.made(), 1)
        self.assertEqual(fns_benchmark._cache_version(self.cache), fns_benchmark._CACHE_VERSION)


__all__ = ["InterpolationTests", "DeckTests", "StatsTests", "ReactionCacheTests"]
