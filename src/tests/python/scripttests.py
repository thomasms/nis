"""
Fast tests of the benchmark tools (scripts/nisbench and scripts/bench.py). They need no nuclear data and
no pynis.
"""

import copy
import json
import math
import os
import stat
import struct
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
from pathlib import Path

from .testdata import BENCHMARKS_DIR, FNS_DIR, FNS_REFERENCE, SCRIPTS_DIR

from nisbench import FormatError, format, load_benchmark, load_benchmarks, load_case
from nisbench import engines, results, stats
from nisbench.convert.fns import interpolate_log
from nisbench.deck import deck_is_current, deck_text


def _small_case() -> dict:
    """A case with all the parts of the format: 2 groups, irradiation, then 2 cooling steps."""
    return {
        "format": "nis-benchmark-case", "version": 1, "benchmark": "test", "name": "a/one", "title": "a test",
        "sample": {"mass_kg": 0.002, "atoms": {"Co59": 1e22, "Ag110m": 1e3}},
        "steps": {"dt_s": [100.0, 10.0, 20.0], "flux_n_m2_s": [1e14, 0.0, 0.0]},
        "spectrum": {"energy_bounds_ev": [1.0, 10.0, 100.0], "values": [0.25, 0.75]},
        "measurements": [{"quantity": "decay_heat", "unit": "uW/g", "steps": [1, 2], "values": [3.0, 2.0],
                          "uncertainties": [0.3, 0.2], "calculations": {"Other code": [3.3, None]}}],
    }


class FormatTests(unittest.TestCase):
    def test_parse(self) -> None:
        case = format.parse_case(_small_case())
        self.assertEqual(case.steps, [(100.0, 1e14), (10.0, 0.0), (20.0, 0.0)])
        self.assertEqual(case.spectrum, ([1.0, 10.0, 100.0], [0.25, 0.75]))
        self.assertEqual(case.end_times(), [100.0, 110.0, 130.0])
        # the cooling time is after the end of the last step with a flux
        self.assertEqual(case.cooling_times([1, 2]), [10.0, 30.0])
        self.assertEqual(case.measurements[0].calculations, {"Other code": [3.3, None]})

    def test_round_trip(self) -> None:
        case = format.parse_case(_small_case())
        self.assertEqual(format.case_to_json(format.parse_case(format.case_to_json(case))), format.case_to_json(case))

    def test_write_and_load(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "a" / "one.json"
            format.write_case(format.parse_case(_small_case()), path)
            case = load_case(path)
            self.assertEqual(case.path, path)
            self.assertEqual(case.deck_path, path.with_suffix(".inp"))
            # the lists of numbers are on one line
            self.assertIn('"dt_s": [100.0, 10.0, 20.0]', path.read_text())

    def test_no_irradiation(self) -> None:
        data = _small_case()
        data["steps"]["flux_n_m2_s"] = [0.0, 0.0, 0.0]
        data["spectrum"] = None
        case = format.parse_case(data)
        self.assertEqual(case.cooling_times([1, 2]), [110.0, 130.0])

    def test_group_structure(self) -> None:
        bounds = format.read_group_structure("CCFE-709")
        self.assertEqual(len(bounds), 710)
        data = _small_case()
        data["spectrum"] = {"group_structure": "CCFE-709", "values": [1.0] * 709}
        self.assertEqual(format.parse_case(data).spectrum[0], bounds)

    def errors(self, change) -> str:
        data = copy.deepcopy(_small_case())
        change(data)
        with self.assertRaises(FormatError) as e:
            format.parse_case(data)
        return str(e.exception)

    def test_bad_cases(self) -> None:
        def set_(path, value):
            def change(data):
                target = data
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
            return change

        checks = [
            (set_(["version"], 2), '"version" must be 1'),
            (set_(["format"], "x"), '"format" must be'),
            (set_(["sample", "atoms"], {"Cobalt": 1.0}), "not a nuclide name"),
            (set_(["sample", "atoms"], {}), "not empty"),
            (set_(["steps", "dt_s"], [100.0, -1.0, 20.0]), "dt_s must be more than 0"),
            (set_(["steps", "flux_n_m2_s"], [1e14, 0.0]), "one flux for each step"),
            (set_(["spectrum"], None), "needs a spectrum"),
            (set_(["spectrum", "energy_bounds_ev"], [1.0, 10.0]), "needs 3 energy bounds"),
            (set_(["spectrum", "energy_bounds_ev"], [1.0, 100.0, 10.0]), "must go up"),
            (set_(["measurements", 0, "quantity"], "heat"), "unknown quantity"),
            (set_(["measurements", 0, "unit"], "W/cm3"), "unknown unit"),
            (set_(["measurements", 0, "steps"], [1, 3]), "step index from 0 to 2"),
            (set_(["measurements", 0, "steps"], [2, 1]), "must go up"),
            (set_(["measurements", 0, "values"], [1.0]), '"values" has 1 items'),
            (set_(["measurements", 0, "calculations"], {"x": [1.0]}), "one value (or null) for each step"),
            (set_(["measurements", 0, "unit"], "MeV/fission"), 'needs "fissions"'),
            (set_(["sample", "mass_kg"], None), "needs sample.mass_kg"),
            (set_(["measurements", 0, "values"], [1.0, float("inf")]), "must be a number"),
        ]
        for change, message in checks:
            with self.subTest(message=message):
                self.assertIn(message, self.errors(change))


class FnsBenchmarkTests(unittest.TestCase):
    def test_all_cases_load(self) -> None:
        benchmark = load_benchmark(FNS_DIR)
        self.assertEqual(benchmark.category, "fusion")
        self.assertEqual(len(benchmark.cases), 132)
        for case in benchmark.cases:
            self.assertEqual(len(case.spectrum[1]), 709)
            self.assertGreater(case.steps[0][1], 0.0)
            self.assertTrue(all(flux == 0.0 for _, flux in case.steps[1:]))
            # one cooling step to each measured time
            self.assertEqual(case.measurements[0].steps, list(range(1, len(case.steps))))

    def test_decks_are_made_from_the_cases(self) -> None:
        for case in load_benchmark(FNS_DIR).cases:
            self.assertTrue(deck_is_current(case), f"{case.deck_path} is not made from {case.path}")

    def test_select_cases(self) -> None:
        names = [c.name for c in load_benchmark(FNS_DIR, ["Co", "Fe/1996exp_7hour"]).cases]
        self.assertEqual(names, ["Co/1996exp_5min", "Co/1996exp_7hour", "Co/2000exp_5min", "Fe/1996exp_7hour"])
        # a prefix is a full part of the name: "C" is not "Co"
        self.assertEqual(load_benchmark(FNS_DIR, ["C"]).cases, [])

    def test_all_benchmarks(self) -> None:
        self.assertIn("fns", [b.name for b in load_benchmarks()])
        with self.assertRaises(FormatError):
            load_benchmarks(["no_such_benchmark"])

    def test_reference_has_all_cases(self) -> None:
        rows, calculations = results.read_results([FNS_REFERENCE])
        self.assertEqual(calculations, ["FISPACT-II (TENDL-2017)"])
        cases = load_benchmark(FNS_DIR).cases
        self.assertEqual({r["case"] for r in rows}, {c.name for c in cases})
        self.assertEqual(len(rows), sum(len(c.measurements[0].values) for c in cases))


class DeckTests(unittest.TestCase):
    def test_deck_text(self) -> None:
        case = format.parse_case(_small_case())
        text = deck_text(case)
        lines = [line for line in text.splitlines() if line and not line.startswith("#")]
        self.assertEqual(lines[:3], ["dt 1.000000000e+02 1.000000000e+14", "dt 1.000000000e+01 0",
                                     "dt 2.000000000e+01 0"])
        self.assertIn("    Ag110m 1.000000000e+03", lines)
        self.assertEqual(lines[-3:], ["spectrum", "    1.000000000e+00 1.000000000e+01 1.000000000e+02",
                                      "    2.500000000e-01 7.500000000e-01"])

    def test_no_spectrum(self) -> None:
        data = _small_case()
        data["steps"]["flux_n_m2_s"] = [0.0, 0.0, 0.0]
        data["spectrum"] = None
        lines = [line for line in deck_text(format.parse_case(data)).splitlines() if not line.startswith("#")]
        self.assertNotIn("spectrum", lines)


class ResultsTests(unittest.TestCase):
    def test_units(self) -> None:
        case = format.parse_case(_small_case())
        step = {"decay_heat": 2.0, "beta_heat": 1.0, "gamma_heat": 1.0, "alpha_heat": 0.0, "activity": 1e6}
        # 2 W in 2 g is 1 W/g = 1e6 uW/g
        self.assertAlmostEqual(results.nis_value(case, "decay_heat", "uW/g", step, 10.0), 1e6)
        self.assertAlmostEqual(results.nis_value(case, "decay_heat", "W", step, 10.0), 2.0)
        self.assertAlmostEqual(results.nis_value(case, "activity", "Bq/g", step, 10.0), 5e5)
        case.fissions = 1e10
        mev = 1.602176634e-13
        self.assertAlmostEqual(results.nis_value(case, "beta_heat", "MeV/s/fission", step, 10.0), 1.0 / mev / 1e10)
        # t * P(t)
        self.assertAlmostEqual(results.nis_value(case, "beta_heat", "MeV/fission", step, 10.0), 10.0 / mev / 1e10)

    def test_rows_and_file(self) -> None:
        case = format.parse_case(_small_case())
        steps = [{"decay_heat": h, "beta_heat": 0, "gamma_heat": 0, "alpha_heat": 0, "activity": 0}
                 for h in (9e-6, 6e-6, 4e-6)]
        rows = results.case_rows(case, steps)
        self.assertEqual([r["time_s"] for r in rows], ["10", "30"])
        self.assertEqual([float(r["nis"]) for r in rows], [3.0, 2.0])
        self.assertEqual([r["Other code"] for r in rows], ["3.300000e+00", ""])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r.csv"
            results.write_results(path, rows, results.calculation_names([case]))
            read, calculations = results.read_results([path])
        self.assertEqual(calculations, ["Other code"])
        self.assertEqual(read[0]["nis"], 3.0)
        self.assertTrue(math.isnan(read[1]["Other code"]))

    def test_nis_output(self) -> None:
        output = textwrap.dedent("""\
            step 1
                time 300
                flux 1e+14

            nuclides
                Co60 1e10

            activity
                total 7.0

            decay_heat
                total 2.5e-03
                light_particle 1e-3
                electromagnetic 1.5e-3
                heavy_particle 0

            step 2
            decay_heat
                total 1.25e-03
            """)
        steps = engines.parse_nis_output(output)
        self.assertEqual(len(steps), 2)
        self.assertEqual(steps[0], {"activity": 7.0, "decay_heat": 2.5e-3, "beta_heat": 1e-3,
                                    "gamma_heat": 1.5e-3, "alpha_heat": 0.0})
        self.assertEqual(steps[1], {"decay_heat": 1.25e-3})


class StatsTests(unittest.TestCase):
    def test_fraction_within(self) -> None:
        # symmetric in log: 1.25 and 0.8 are both within 25 %
        self.assertEqual(stats.fraction_within([1.25, 0.8, 2.0, 1.0], 0.25), 0.75)

    def test_geometric_mean(self) -> None:
        self.assertAlmostEqual(stats.geometric_mean([0.5, 2.0]), 1.0)

    def test_same_points_for_all_codes(self) -> None:
        nan = math.nan
        rows = [{"benchmark": "b", "case": "c", "quantity": "decay_heat", "experiment": e, "nis": n, "X": x}
                for e, n, x in ((1.0, 2.0, 1.0), (1.0, 4.0, nan), (0.0, 1.0, 1.0))]
        ratios = stats.case_ratios(rows, ["X"])
        # only the first point has the experiment, NIS and X
        self.assertEqual(ratios[("b", "decay_heat")]["c"], {"nis": [2.0], "X": [1.0]})

    def test_report(self) -> None:
        rows, calculations = results.read_results([FNS_REFERENCE])
        text = stats.report(rows, calculations, worst=3)
        self.assertIn("fns, decay_heat: 132 cases", text)


class InterpolationTests(unittest.TestCase):
    def test_interpolate_log_at_the_points(self) -> None:
        times = [1.0, 2.0, 4.0]
        values = [100.0, 10.0, 1.0]
        for t, v in zip(times, values):
            self.assertAlmostEqual(interpolate_log(times, values, t), v)

    def test_interpolate_log_between_the_points(self) -> None:
        # the log is linear in time, so the middle is the geometric mean
        self.assertAlmostEqual(interpolate_log([0.0, 2.0], [100.0, 1.0], 1.0), 10.0)

    def test_interpolate_log_outside(self) -> None:
        self.assertTrue(math.isnan(interpolate_log([1.0, 2.0], [1.0, 2.0], 0.5)))
        self.assertTrue(math.isnan(interpolate_log([1.0, 2.0], [1.0, 2.0], 3.0)))
        # no log of 0
        self.assertTrue(math.isnan(interpolate_log([1.0, 2.0], [0.0, 2.0], 1.5)))


class CommandLineTests(unittest.TestCase):
    def bench(self, *arguments) -> subprocess.CompletedProcess:
        return subprocess.run([sys.executable, str(SCRIPTS_DIR / "bench.py"), *map(str, arguments)],
                              capture_output=True, text=True)

    def test_list_and_check(self) -> None:
        result = self.bench("list")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("fns", result.stdout)
        result = self.bench("check", "fns")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_check_finds_a_changed_deck(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "test").mkdir()
            (root / "test" / "benchmark.json").write_text(json.dumps({
                "format": "nis-benchmark", "version": 1, "name": "test", "title": "Test", "category": "other",
                "description": "a test", "source": "made by the test"}))
            format.write_case(format.parse_case(_small_case() | {"benchmark": "test"}), root / "test" / "a" / "one.json")
            self.assertEqual(self.bench("check", "--root", root).returncode, 1)
            self.assertEqual(self.bench("decks", "--root", root).returncode, 0)
            self.assertEqual(self.bench("check", "--root", root).returncode, 0)
            # a bad case file gives an error, not a traceback
            (root / "test" / "a" / "one.json").write_text("{}")
            result = self.bench("list", "--root", root)
            self.assertEqual(result.returncode, 2)
            self.assertIn("error:", result.stderr)


def _write_cache(path: Path, version: int = engines._CACHE_VERSION) -> None:
    path.write_bytes(engines._CACHE_MAGIC + struct.pack("=I", version) + b"\0" * 8)


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
                f.write({engines._CACHE_MAGIC!r} + struct.pack("=I", {engines._CACHE_VERSION}) + bytes(8))
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
        self.assertEqual(engines._cache_version(self.cache), 7)
        self.assertIsNone(engines._cache_version(self.endf))
        self.assertIsNone(engines._cache_version(self.directory / "missing"))

    def test_a_cache_is_used_as_it_is(self) -> None:
        _write_cache(self.cache)
        self.assertEqual(engines.reaction_cache(self.cache, self.nis_cache), self.cache)
        self.assertEqual(self.made(), 0)

    def test_made_if_missing(self) -> None:
        self.assertEqual(engines.reaction_cache(self.endf, self.nis_cache), self.cache)
        self.assertEqual(self.made(), 1)
        self.assertEqual(engines._cache_version(self.cache), engines._CACHE_VERSION)
        self.assertFalse(self.cache.with_suffix(".nisreac.tmp").exists())

        # it is up to date now, so it is not made again
        engines.reaction_cache(self.endf, self.nis_cache)
        self.assertEqual(self.made(), 1)

    def test_made_if_older_than_the_data(self) -> None:
        _write_cache(self.cache)
        old = time.time() - 3600
        os.utime(self.cache, (old, old))
        engines.reaction_cache(self.endf, self.nis_cache)
        self.assertEqual(self.made(), 1)

    def test_made_if_another_version(self) -> None:
        _write_cache(self.cache, engines._CACHE_VERSION + 1)
        engines.reaction_cache(self.endf, self.nis_cache)
        self.assertEqual(self.made(), 1)
        self.assertEqual(engines._cache_version(self.cache), engines._CACHE_VERSION)


__all__ = ["FormatTests", "FnsBenchmarkTests", "DeckTests", "ResultsTests", "StatsTests",
           "InterpolationTests", "CommandLineTests", "ReactionCacheTests"]
