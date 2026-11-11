"""
Tests of the data processing (processing/nisdata and processing/process.py).

The fast tests need nothing. The pipeline tests need NJOY and PREPRO (processing/process.py tools, or
NJOY_EXE and PREPRO_EXE_DIR) and the NJOY2016 test files, which 'tools' clones (or NISDATA_TEST_ENDF_DIR).
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from .testdata import PROCESSING_DIR, ROOT

from nisdata import decks, endf
from nisdata.pipeline import Tools, check_output, error_text, process_file
from nisdata.runner import check_unique, process_all
from nisdata.tools import default_tools

PROCESS = PROCESSING_DIR / "process.py"


def _line(data: str, mat: int, mf: int, mt: int, n: int = 1) -> str:
    return f"{data:66s}{mat:4d}{mf:2d}{mt:3d}{n:5d}\n"


def _fake_endf(za: int, mat: int, liso: int = 0, nsub: int = 10, tpid: str = "a test tape") -> str:
    """The start of an ENDF file: TPID, the first 3 lines of MF1/MT451, a cross section and the end lines."""
    return (
        _line(tpid, 1, 0, 0, 0)
        + _line(f"{za:11.1f} 1.000000+0          0          0          0          0", mat, 1, 451)
        + _line(" 0.000000+0 0.000000+0          0" + f"{liso:11d}" + "          0          6", mat, 1, 451, 2)
        + _line(" 1.000000+0 2.000000+7          0          0" + f"{nsub:11d}" + "          0", mat, 1, 451, 3)
        + _line(" 1.0", mat, 3, 1, 1)
        + _line("", mat, 3, 0, 99999)
        + _line("", mat, 0, 0, 0)
        + _line("", 0, 0, 0, 0)
        + _line("", -1, 0, 0, 0)
    )


class EndfTextTests(unittest.TestCase):
    def setUp(self) -> None:
        self.dir = Path(tempfile.mkdtemp())

    def write(self, name: str, text: str) -> Path:
        path = self.dir / name
        path.write_text(text)
        return path

    def test_numbers(self) -> None:
        self.assertEqual(endf.endf_float(" 2.605600+4"), 26056.0)
        self.assertEqual(endf.endf_float("-1.5-3"), -0.0015)
        self.assertEqual(endf.endf_float("1.0E+2"), 100.0)
        self.assertEqual(endf.endf_float("   "), 0.0)
        self.assertEqual(endf.endf_int("    12"), 12)

    def test_read_material(self) -> None:
        m = endf.read_material(self.write("fe", _fake_endf(26056, 2631)))
        self.assertEqual((m.mat, m.za, m.liso, m.nsub), (2631, 26056, 0, 10))
        self.assertEqual((m.z, m.a, m.zai, m.name), (26, 56, 260560, "Fe56"))
        isomer = endf.read_material(self.write("ag", _fake_endf(47110, 4746, liso=1)))
        self.assertEqual((isomer.name, isomer.zai), ("Ag110m", 471101))
        self.assertEqual(endf.Material(0, 72178, 2, 10).name, "Hf178n")

    def test_not_endf(self) -> None:
        with self.assertRaises(ValueError):
            endf.read_material(self.write("readme.txt", "not an ENDF file\n"))

    def test_find_neutron_files(self) -> None:
        self.write("n-Fe056.tendl", _fake_endf(26056, 2631))
        (self.dir / "sub").mkdir()
        self.write("sub/n-Co059.tendl", _fake_endf(27059, 2725))
        self.write("tsl-HinH2O.endf", _fake_endf(1, 1, nsub=12))
        self.write("README", "text\n")
        names = [p.name for p in endf.find_endf_files(self.dir)]
        self.assertEqual(sorted(names), ["n-Co059.tendl", "n-Fe056.tendl"])

    def test_combine(self) -> None:
        files = [self.write("b", _fake_endf(27059, 2725)), self.write("a", _fake_endf(26056, 2631)),
                 self.write("c", _fake_endf(27058, 2722, liso=1))]
        out = self.dir / "all.reac"
        self.assertEqual(endf.combine(files, out, "My library"), 3)
        lines = out.read_text().splitlines()
        self.assertTrue(lines[0].startswith("My library"))
        self.assertEqual(endf.mat_mf_mt(lines[0]), (1, 0, 0))
        self.assertEqual(endf.mat_mf_mt(lines[-1]), (-1, 0, 0))
        # one TPID and one TEND, the materials ordered by ZA and isomer state
        self.assertEqual(sum(endf.mat_mf_mt(line)[0] == -1 for line in lines), 1)
        mats = [endf.mat_mf_mt(line)[0] for line in lines if endf.mat_mf_mt(line)[1:] == (1, 451)]
        self.assertEqual(sorted(set(mats), key=mats.index), [2631, 2722, 2725])
        # each material ends with its MEND line (MAT 0)
        self.assertEqual(sum(endf.mat_mf_mt(line)[0] == 0 for line in lines), 3)
        self.assertTrue(all(len(line) == 80 for line in lines))

    def test_unique_nuclides(self) -> None:
        a = self.write("a", _fake_endf(26056, 2631))
        b = self.write("b", _fake_endf(26056, 2631, tpid="another library"))
        c = self.write("c", _fake_endf(27059, 2725))
        check_unique([a, c])
        with self.assertRaisesRegex(ValueError, "Fe56"):
            check_unique([a, b, c])


class DeckTests(unittest.TestCase):
    def test_fixed_format(self) -> None:
        # the same as the UKAEA inputs: the codes read these lines by column
        self.assertEqual(decks.sigma1("A", "B").splitlines()[0],
                         "          0          0 293.6      1.00000-10          1          0")
        self.assertEqual(decks.recent("A", "B").splitlines()[0],
                         "          0 1.00000-10          1          1          1          1")
        self.assertEqual(decks.groupie("A", "B").splitlines()[0],
                         "         0       1102         0          -1 1.00000-03           0   1")

    def test_groups(self) -> None:
        bounds = [float(v) for v in decks.GROUPS_FILE.read_text().split()]
        self.assertEqual(len(bounds), decks.NGROUPS + 1)
        self.assertEqual((bounds[0], bounds[-1]), (1e-5, 1e9))
        self.assertTrue(all(b < c for b, c in zip(bounds, bounds[1:])))

    def test_njoy(self) -> None:
        deck = decks.njoy(2631, 300.0)
        self.assertIn("'pendf for material 2631'/", deck)
        self.assertIn("2631 1/", deck)
        self.assertEqual(deck.count("300.0/"), 2)
        self.assertTrue(deck.rstrip().endswith("stop"))

    def test_merger(self) -> None:
        lines = decks.merger("OUT", ["IN1", "IN2"], decks.MERGE_BOTH).splitlines()
        self.assertEqual(lines[:5], ["OUT", " " * 69 + "0 0", "IN1", "IN2", "END"])


class ErrorTextTests(unittest.TestCase):
    def test_first_error(self) -> None:
        log = "start\n ERROR - At least one Spin-Pair# MUST = 2\n   Check and Correct Data\n more\n"
        self.assertEqual(error_text(log), "ERROR - At least one Spin-Pair# MUST = 2 Check and Correct Data")
        self.assertEqual(error_text("a\nb\nlast line\n"), "last line")
        self.assertEqual(error_text(""), "no output")


class MissingToolsTests(unittest.TestCase):
    def test_missing_programs(self) -> None:
        tools = Tools(njoy=Path("/nonexistent/njoy"), prepro_dir=Path("/nonexistent"))
        with self.assertRaisesRegex(FileNotFoundError, "njoy"):
            tools.check()
        with self.assertRaises(FileNotFoundError):
            process_all([], Path(tempfile.mkdtemp()), tools)


def _tools() -> Tools:
    if os.environ.get("NJOY_EXE") and os.environ.get("PREPRO_EXE_DIR"):
        return Tools(Path(os.environ["NJOY_EXE"]), Path(os.environ["PREPRO_EXE_DIR"]))
    return default_tools(PROCESSING_DIR / "tools")


def _test_files_dir() -> Path:
    if os.environ.get("NISDATA_TEST_ENDF_DIR"):
        return Path(os.environ["NISDATA_TEST_ENDF_DIR"])
    return PROCESSING_DIR / "tools" / "src" / "NJOY2016" / "tests" / "resources"


def _have_tools() -> bool:
    try:
        _tools().check()
    except FileNotFoundError:
        return False
    return _test_files_dir().is_dir()


# two small files: a TENDL-2023 file and an isomer
TEST_FILES = ("n-018_Ar_37-tendl2023.endf", "n-027_Co_058m1-ENDF8.0.endf")


@unittest.skipUnless(_have_tools(), "needs NJOY and PREPRO (processing/process.py tools) and the NJOY2016 test files")
class PipelineTests(unittest.TestCase):
    """Process real ENDF files with NJOY and PREPRO."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.dir = Path(tempfile.mkdtemp())
        cls.raw = cls.dir / "raw"
        cls.raw.mkdir()
        for name in TEST_FILES:
            (cls.raw / name).symlink_to(_test_files_dir() / name)
        cls.out = cls.dir / "out"
        cls.results = process_all(sorted(cls.raw.iterdir()), cls.out, _tools(), jobs=2, quiet=True)

    def test_all_done(self) -> None:
        self.assertEqual(sorted((r.name, r.ok) for r in self.results), [("Ar37", True), ("Co58m", True)],
                         [r.error for r in self.results])
        summary = json.loads((self.out / "summary.json").read_text())
        self.assertEqual((summary["files"], summary["processed"], summary["failed"]), (2, 2, 0))

    def test_output(self) -> None:
        for name, za, liso in (("Ar37", 18037, 0), ("Co58m", 27058, 1)):
            path = self.out / f"{name}.asc"
            check_output(path, endf.read_material(path))
            material = endf.read_material(path)
            self.assertEqual((material.za, material.liso), (za, liso))
            sections = {endf.mat_mf_mt(line)[1:] for line in path.read_text().splitlines()}
            # the total, elastic and capture, kerma and damage from NJOY HEATR, gas production from GASPR
            for section in [(3, 1), (3, 2), (3, 102), (3, 301), (3, 444), (3, 203)]:
                self.assertIn(section, sections, f"{name} has no MF{section[0]}/MT{section[1]}")

    def test_done_nuclides_are_skipped(self) -> None:
        again = process_all(sorted(self.raw.iterdir()), self.out, _tools(), jobs=2, quiet=True)
        self.assertTrue(all(r.ok and r.skipped for r in again))

    def test_failure_is_reported(self) -> None:
        bad = self.dir / "bad"
        bad.mkdir()
        # the header is good but there is no data
        (bad / "n-bad.endf").write_text(_fake_endf(26056, 2631))
        out = self.dir / "bad_out"
        result = process_file(bad / "n-bad.endf", out, _tools())
        self.assertFalse(result.ok)
        self.assertTrue(result.error)
        self.assertTrue((out / "logs" / "Fe56.log").is_file())
        self.assertFalse((out / "Fe56.asc").exists())

    def test_combine_and_read_with_nis(self) -> None:
        output = self.dir / "test.reac"
        run = subprocess.run([sys.executable, str(PROCESS), "combine", str(self.out), "--output", str(output),
                              "--title", "test", "--nis-cache", ""], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        try:
            import pynis
        except ImportError:
            self.skipTest("needs pynis")
        data = pynis.core.reaction_data_from_file(str(output))
        # thermal capture of Ar37 gives Ar38
        index = {zai: i for i, zai in enumerate(pynis.data.ALL_ZAIS())}
        spectrum = pynis.core.Spectrum([0.0252, 0.0254], [1e14])
        mat = pynis.core.reaction_mat(data, spectrum)
        ar37, ar38 = pynis.zai.get_from_name("Ar37"), pynis.zai.get_from_name("Ar38")
        self.assertGreater(mat[index[ar37], index[ar38]], 0.0)


class CommandLineTests(unittest.TestCase):
    def run_process(self, *args) -> subprocess.CompletedProcess:
        return subprocess.run([sys.executable, str(PROCESS), *map(str, args)], capture_output=True, text=True,
                              cwd=ROOT)

    def test_help(self) -> None:
        run = self.run_process("--help")
        self.assertEqual(run.returncode, 0)
        for command in ("tools", "fetch", "run", "combine", "tendl"):
            self.assertIn(command, run.stdout)

    def test_duplicates_are_an_error(self) -> None:
        folder = Path(tempfile.mkdtemp())
        (folder / "a.endf").write_text(_fake_endf(26056, 2631))
        (folder / "b.endf").write_text(_fake_endf(26056, 2631))
        tools = Path(tempfile.mkdtemp())
        (tools / "njoy").touch(mode=0o755)
        for code in ("linear", "recent", "sigma1", "sixpak", "merger", "activate", "fixup", "dictin", "groupie"):
            (tools / code).touch(mode=0o755)
        run = self.run_process("run", folder, "--output", folder / "out", "--njoy", tools / "njoy",
                               "--prepro-dir", tools)
        self.assertEqual(run.returncode, 2)
        self.assertIn("more than one file for the same nuclide", run.stderr)

    def test_combine_needs_files(self) -> None:
        run = self.run_process("combine", tempfile.mkdtemp(), "--output", Path(tempfile.mkdtemp()) / "x.reac")
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("no processed files", run.stderr)
