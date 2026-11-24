"""
The cross_section quantity of the benchmark format, and the IRDFF-II converter (benchmarks/irdff2):

    - the format: a cross section needs a target and a product, the target is in the sample, it is measured at a
      step with a flux, and a case with one has no deck
    - the converter: the spectrum files (pointwise or groups), the RR_UNC listings of IRDFF-II and measured.csv
    - the engine: the cross section of a case is the same as the collapsed cross sections of the reactions of the
      target that make the product (needs data/jeff4.reac and data/jeff311.dec, processing/process.py bootstrap)
    - the IRDFF-II benchmark: the case files, and a regression test that runs all the cases with pynis and compares
      with benchmarks/irdff2/results/nis_reference.csv. If a change of the results is intended, make the reference
      again:

        uv run python scripts/bench.py run irdff2 --decay-data data/jeff311.dec \\
            --reaction-data data/jeff4.reac --output benchmarks/irdff2/results/nis_reference.csv
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
from nisbench.convert import irdff2


def _cross_section_case() -> dict:
    return {
        "format": "nis-benchmark-case", "version": 1, "benchmark": "test", "name": "Cf252/Al27-Na24",
        "sample": {"atoms": {"Al27": 1.0}},
        "steps": {"dt_s": [1.0], "flux_n_m2_s": [1.0]},
        "spectrum": {"energy_bounds_ev": [1.0e6, 1.0e7, 2.0e7], "values": [1.0, 1.0]},
        "measurements": [{"quantity": "cross_section", "unit": "mb", "target": "Al27", "product": "Na24",
                          "all_states": True, "steps": [0], "values": [1.0], "uncertainties": [0.02]}],
    }


class CrossSectionFormatTests(unittest.TestCase):
    def errors(self, change) -> str:
        data = copy.deepcopy(_cross_section_case())
        change(data)
        with self.assertRaises(FormatError) as context:
            format.parse_case(data)
        return str(context.exception)

    def test_case(self) -> None:
        case = format.parse_case(_cross_section_case())
        m = case.measurements[0]
        self.assertEqual((m.target, m.product, m.all_states), ("Al27", "Na24", True))
        # the nis program does not give cross sections
        self.assertFalse(case.has_deck)
        self.assertEqual(format.parse_case(json.loads(json.dumps(format.case_to_json(case)))), case)

    def test_all_states_is_false_by_default(self) -> None:
        data = _cross_section_case()
        del data["measurements"][0]["all_states"]
        self.assertFalse(format.parse_case(data).measurements[0].all_states)

    def test_units(self) -> None:
        case = format.parse_case(_cross_section_case())
        step = {"cross_section": {("Al27", "Na24", True): 1.0e-27}}
        self.assertAlmostEqual(results.nis_value(case, "cross_section", "mb", step, 0.0, ("Al27", "Na24", True)), 1.0)
        self.assertAlmostEqual(results.nis_value(case, "cross_section", "b", step, 0.0, ("Al27", "Na24", True)), 1e-3)

    def test_bad_cases(self) -> None:
        m = lambda d: d["measurements"][0]
        checks = [
            (lambda d: m(d).pop("target"), 'needs "target"'),
            (lambda d: m(d).update(product="sodium"), 'needs "product"'),
            (lambda d: m(d).update(target="Al28"), "must be in sample.atoms"),
            (lambda d: m(d).update(product="Al26m"), '"all_states" needs the ground state'),
            (lambda d: m(d).update(all_states="yes"), "must be true or false"),
            (lambda d: m(d).update(unit="Bq"), "unknown unit"),
            (lambda d: d["steps"].update(flux_n_m2_s=[0.0]), "at the end of a step with a flux"),
            (lambda d: m(d).update(quantity="activity", unit="Bq"), "only for a cross_section"),
        ]
        for change, expected in checks:
            with self.subTest(expected=expected):
                self.assertIn(expected, self.errors(change))


def _endf_lines(mat: int, energies: list, values: list, law: int) -> str:
    """The MF3/MT261 of a spectrum, as in the IRDFF-II files."""
    def line(fields: str, seq: int) -> str:
        return f"{fields:66s}{mat:4d}{3:2d}{261:3d}{seq:5d}\n"

    def number(v: float) -> str:
        return f"{v:11.4E}"

    rows = [line(f"{number(0.0)}{number(0.0)}{0:11d}{0:11d}{1:11d}{0:11d}", 1),
            line(f"{number(0.0)}{number(0.0)}{0:11d}{0:11d}{1:11d}{len(energies):11d}", 2),
            line(f"{len(energies):11d}{law:11d}", 3)]
    pairs = [number(v) for xy in zip(energies, values) for v in xy]
    for k in range(0, len(pairs), 6):
        rows.append(line("".join(pairs[k:k + 6]), 4 + k // 6))
    return "".join(rows) + f"{'':66s}{mat:4d}{3:2d}{0:3d}{99999:5d}\n"


_LISTING = """
   12    Spectrum MAT No.         :   9861
   57      36 12000      11024g 8.26095 1.6608E+00 +/- 2.958E-02      0.80   1.59   1.78
   60      39 13027         103 5.84274 4.7445E+00 +/- 1.116E-01      2.05   1.15   2.35
   61      40 13027         107 8.66807 1.0167E+00 +/- 1.804E-02      0.71   1.63   1.77
   63      42 13027      13026g 15.9694 6.8607E-03 +/- 6.958E-04      3.81   9.40 10.14
   64      43 13027      13026m 17.9578 8.6409E-04 +/- 2.211E-04      9.44  23.79 25.59
   22       1 3000            1 1.87387 1.8353E+03 +/- 2.084E+00      0.00   0.11   0.11
   28        7 3006     105 0.65237 3.2063E+02 +/- 3.341E+00      0.89   0.55   1.04
   29        8 3006     205 0.65237 3.2063E+02 +/- 3.341E+00      0.89   0.55   1.04
   30        9 3006     207 1.70993 4.8215E+02 +/- 3.236E+01      6.70   0.40   6.71
"""


class Irdff2ConverterTests(unittest.TestCase):
    def test_products(self) -> None:
        self.assertEqual(irdff2.product_state("Na24"), ("Na24", True))
        self.assertEqual(irdff2.product_state("Al26g"), ("Al26", False))
        self.assertEqual(irdff2.product_state("Al26m"), ("Al26m", False))
        for token in ("Na24", "Al26g", "Al26m"):
            self.assertEqual(irdff2.product_token(*irdff2.product_state(token)), token)
        with self.assertRaises(ValueError):
            irdff2.product_state("Al26x")

    def test_group_spectrum(self) -> None:
        # lin-lin: the trapezoids, the interval from 0 eV and the zeros at the end are left out
        bounds, groups = irdff2.group_spectrum([0.0, 1.0, 3.0, 4.0, 5.0], [0.0, 2.0, 4.0, 0.0, 0.0], [2], [5])
        self.assertEqual(bounds, [1.0, 3.0, 4.0])
        self.assertEqual(groups, [6.0, 2.0])
        # histogram: the value of the left point
        bounds, groups = irdff2.group_spectrum([1.0, 2.0, 4.0], [3.0, 1.0, 0.0], [1], [3])
        self.assertEqual(groups, [3.0, 2.0])
        # a step (two points at the same energy) has no width
        bounds, groups = irdff2.group_spectrum([1.0, 2.0, 2.0, 3.0], [1.0, 1.0, 3.0, 3.0], [2], [4])
        self.assertEqual((bounds, groups), ([1.0, 2.0, 3.0], [1.0, 3.0]))
        with self.assertRaises(ValueError):
            irdff2.group_spectrum([1.0, 2.0], [1.0, 1.0], [5], [2])

    def test_read_spectrum(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "irdff_sp_test.endf"
            energies = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
            path.write_text("a title line\n" + _endf_lines(9861, energies, [1.0] * 7, 2))
            mat, x, y, laws, breaks = irdff2.read_spectrum(path)
            self.assertEqual((mat, x, y, laws, breaks), (9861, energies, [1.0] * 7, [2], [7]))

    def test_read_listing(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "irdff_sp_test_RR_UNC.txt"
            path.write_text(_LISTING)
            mat, values = irdff2.read_rr_unc(path)
        self.assertEqual(mat, 9861)
        # the MT gives the product in all its states, g and m one state; elements (12000) are left out. He4 is the
        # gas production (MT 207), not the He4 of Li6(n,t) (MT 105)
        self.assertEqual(values, {("Al27", "Mg27", True): 4.7445, ("Al27", "Na24", True): 1.0167,
                                  ("Al27", "Al26", False): 6.8607e-3, ("Al27", "Al26m", False): 8.6409e-4,
                                  ("Li6", "H3", True): 320.63, ("Li6", "He4", True): 482.15})

    def test_case_spectrum(self) -> None:
        # a histogram on the IRDFF-725 bounds keeps all the groups, and gives the name
        bounds = format.read_group_structure("IRDFF-725")
        values = [1.0] * (len(bounds) - 2) + [0.0, 0.0]
        result_bounds, groups, name = irdff2.case_spectrum(bounds, values, [1], [len(bounds)])
        self.assertEqual(name, "IRDFF-725")
        self.assertEqual(len(groups), len(bounds) - 1)
        self.assertAlmostEqual(groups[0], bounds[1] - bounds[0])
        self.assertEqual(groups[-1], 0.0)
        # other bounds: the groups of group_spectrum
        self.assertEqual(irdff2.case_spectrum([1.0, 2.0, 4.0], [3.0, 1.0, 0.0], [1], [3]), ([1.0, 2.0, 4.0], [3.0, 2.0], None))

    def test_convert(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder)
            # the pointwise form is not used when there is a groups form
            (source / "irdff_sp_Cf252_ln.endf").write_text(_endf_lines(9861, [1e5, 2e7], [5.0, 0.0], 2))
            (source / "irdff_sp_Cf252_ln.g").write_text(_endf_lines(9861, [1e5, 1e6, 1e7, 2e7], [1.0] * 4, 1))
            (source / "irdff_sp_Cf252_ln_RR_UNC.txt").write_text(_LISTING)
            (source / "measured.csv").write_text(
                "# a comment\nspectrum,table,notation,target,product,value_mb,uncertainty_percent,reference\n"
                'Cf252,18,Al27a,Al27,Na24,1.006,2.0,"[1, 2]"\nCf252,,,Al27,Al26m,0.0009,20,\n')
            cases = irdff2.convert(source)
            self.assertEqual([c.name for c in cases], ["Cf252/Al27-Na24", "Cf252/Al27-Al26m"])
            na24 = cases[0]
            self.assertEqual(na24.spectrum, ([1e5, 1e6, 1e7, 2e7], [9e5, 9e6, 1e7]))
            self.assertEqual(na24.measurements[0].values, [1.006])
            self.assertAlmostEqual(na24.measurements[0].uncertainties[0], 0.02012)
            self.assertEqual(na24.measurements[0].calculations, {"IRDFF-II": [1.0167]})
            self.assertEqual([c.name for c in irdff2.convert(source, ["Cf252/Al27-Al26m"])], ["Cf252/Al27-Al26m"])
            # each case is a valid case file
            for case in cases:
                format.write_case(case, source / "out" / f"{case.name}.json")
            (source / "measured.csv").write_text("spectrum,target,product,value_mb,uncertainty_percent\nU235,Al27,Na24,1,1\n")
            with self.assertRaises(ValueError):
                irdff2.convert(source)


@unittest.skipUnless(has_data(JEFF4_REACTION) and has_data(JEFF311_DECAY),
                     "needs data/jeff4.reac and data/jeff311.dec (processing/process.py bootstrap)")
class CrossSectionEngineTests(unittest.TestCase):
    def run_case(self, target: str, product: str, all_states: bool) -> float:
        from nisbench import engines
        data = _cross_section_case()
        data["sample"]["atoms"] = {target: 1.0}
        data["measurements"][0].update(target=target, product=product, all_states=all_states)
        case = format.parse_case(data)
        nuclear_data = engines.NuclearData(decay=JEFF311_DECAY, reaction=JEFF4_REACTION)
        step = next(iter(engines.run([case], nuclear_data)))[0]
        return results.nis_value(case, "cross_section", "mb", step, 0.0, (target, product, all_states))

    def test_same_as_the_collapsed_cross_section(self) -> None:
        import pynis as nis
        bounds, values = _cross_section_case()["spectrum"].values()
        spectrum = nis.core.Spectrum(bounds, values)
        reaction_data = nis.core.reaction_data_from_file(str(engines_cache()))
        # Al27 makes Na24 only with (n,a), MT107
        expected = sum(nis.core.collapse_cross_section(spectrum, list(xs.energies), list(xs.values))
                       for xs in reaction_data.find(130270, 107)) / sum(values) * 1e27
        self.assertGreater(expected, 1.0)
        self.assertAlmostEqual(self.run_case("Al27", "Na24", True) / expected, 1.0, places=10)

    def test_all_states(self) -> None:
        ground = self.run_case("Nb93", "Nb92", False)
        isomer = self.run_case("Nb93", "Nb92m", False)
        self.assertGreater(isomer, 0.0)
        self.assertAlmostEqual(self.run_case("Nb93", "Nb92", True), ground + isomer, delta=1e-12 * (ground + isomer))


def engines_cache() -> Path:
    from nisbench import engines
    return engines.reaction_cache(JEFF4_REACTION)


IRDFF2_DIR = BENCHMARKS_DIR / "irdff2"
IRDFF2_REFERENCE = IRDFF2_DIR / "results" / "nis_reference.csv"

# the results on another machine or compiler are not the same to the last digit
REFERENCE_REL_TOL = 1e-6


class Irdff2CasesTests(unittest.TestCase):
    def test_all_cases_load(self) -> None:
        benchmark = load_benchmark(IRDFF2_DIR)
        self.assertEqual(benchmark.category, "other")
        rows = irdff2.read_measured(irdff2.MEASURED)
        self.assertEqual(len(benchmark.cases), len(rows))
        self.assertEqual(len(rows), 87)
        for case in benchmark.cases:
            with self.subTest(case=case.name):
                self.assertFalse(case.has_deck)
                self.assertEqual(case.group_structure, "IRDFF-725")
                m, = case.measurements
                self.assertEqual((m.quantity, m.unit, list(case.atoms)), ("cross_section", "mb", [m.target]))
                # each reaction has the value that IRDFF-II calculates, so the product is the right one
                self.assertIsNotNone(m.calculations.get("IRDFF-II", [None])[0])

    def test_reference_has_all_cases(self) -> None:
        rows, _ = results.read_results([IRDFF2_REFERENCE])
        self.assertEqual({r["case"] for r in rows}, {c.name for c in load_benchmark(IRDFF2_DIR).cases})


@unittest.skipUnless(has_data(JEFF4_REACTION) and has_data(JEFF311_DECAY),
                     "needs data/jeff4.reac and data/jeff311.dec (processing/process.py bootstrap)")
class Irdff2RegressionTests(unittest.TestCase):
    def test_same_as_reference(self) -> None:
        from nisbench import engines
        cases = load_benchmark(IRDFF2_DIR).cases
        data = engines.NuclearData(decay=JEFF311_DECAY, reaction=engines.reaction_cache(JEFF4_REACTION))
        rows = []
        for case, step_results in zip(cases, engines.run(cases, data, "pynis", jobs=os.cpu_count() or 1)):
            rows += results.case_rows(case, step_results)
        reference, _ = results.read_results([IRDFF2_REFERENCE])
        expected = {r["case"]: float(r["nis"]) for r in reference}
        self.assertEqual(set(expected), {r["case"] for r in rows})
        changed = [f"{r['case']}: {float(r['nis']):.6e}, reference {expected[r['case']]:.6e}"
                   for r in rows if not math.isclose(float(r["nis"]), expected[r["case"]], rel_tol=REFERENCE_REL_TOL)]
        self.assertFalse(changed, f"{len(changed)} values are not the same as the reference:\n" + "\n".join(changed))
