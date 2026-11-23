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

    def test_from_archive(self) -> None:
        import tarfile
        from nisdata.runner import process_archive
        archive = self.dir / "library.tgz"
        with tarfile.open(archive, "w:gz") as tar:
            for name in TEST_FILES:
                tar.add(_test_files_dir() / name, arcname=f"lib/{name}")
            # not a neutron file, and not an ENDF file: both ignored
            tar.add(_test_files_dir() / "tsl-HinZrH-ENDF8.0.endf", arcname="lib/tsl-HinZrH.endf")
            readme = self.dir / "README"
            readme.write_text("text\n")
            tar.add(readme, arcname="lib/README")
        out = self.dir / "from_archive"
        results = process_archive(archive, out, _tools(), jobs=2, quiet=True, only=("Ar37",))
        self.assertEqual([(r.name, r.ok) for r in results], [("Ar37", True)], [r.error for r in results])
        self.assertTrue(str(results[0].source).endswith("library.tgz:lib/n-018_Ar_37-tendl2023.endf"))
        # the same output as from the file
        self.assertEqual((out / "Ar37.asc").read_text(), (self.out / "Ar37.asc").read_text())

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
        for command in ("bootstrap", "libraries", "decay", "yields", "reactions", "tools", "fetch", "run", "combine", "tendl"):
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


class DownloadTests(unittest.TestCase):
    """fetch.download continues a stopped download with a range request."""

    def test_continue_after_a_drop(self) -> None:
        import http.server
        import threading
        from nisdata import fetch

        content = bytes(range(256)) * 4000
        calls = []

        class Handler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                start = int(self.headers.get("Range", "bytes=0-")[6:].rstrip("-") or 0)
                calls.append(start)
                self.send_response(206 if start else 200)
                if start:
                    self.send_header("Content-Range", f"bytes {start}-{len(content) - 1}/{len(content)}")
                self.send_header("Content-Length", str(len(content) - start))
                self.end_headers()
                # the first time, stop after a part of the file
                end = len(content) // 3 if len(calls) == 1 else len(content)
                self.wfile.write(content[start:end])

            def log_message(self, *args):
                pass

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        threading.Thread(target=server.serve_forever, daemon=True).start()
        self.addCleanup(server.shutdown)
        original_sleep = fetch.time.sleep
        fetch.time.sleep = lambda s: None
        self.addCleanup(setattr, fetch.time, "sleep", original_sleep)

        target = Path(tempfile.mkdtemp()) / "file.tgz"
        fetch.download(f"http://127.0.0.1:{server.server_port}/file.tgz", target, tries=3)
        self.assertEqual(target.read_bytes(), content)
        self.assertEqual(calls, [0, len(content) // 3])
        self.assertFalse(target.with_suffix(".tgz.part").exists())


def _serve(content: bytes, seen_headers: list = None):
    """A local web server that gives the content. Returns its address and a function to stop it."""
    import http.server
    import threading

    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            if seen_headers is not None:
                seen_headers.append(dict(self.headers))
            self.send_response(200)
            self.send_header("Content-Length", str(len(content)))
            self.end_headers()
            self.wfile.write(content)

        def log_message(self, *args):
            pass

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return f"http://127.0.0.1:{server.server_port}", server.shutdown


class LibraryTests(unittest.TestCase):
    """The downloads of nisdata/libraries.py, without the network."""

    def test_the_list(self) -> None:
        from nisdata import libraries
        self.assertIn("jeff311", libraries.LIBRARIES)
        self.assertEqual(libraries.get("jeff4").parts(), ["decay", "nfy", "sfy", "neutron"])
        for lib in libraries.LIBRARIES.values():
            for part in lib.parts():
                source = getattr(lib, part)
                self.assertRegex(source.md5, "^[0-9a-f]{32}$")
                self.assertTrue(source.url.startswith("https://"))
        with self.assertRaises(ValueError):
            libraries.get("jeff5")

    def test_the_first_line_is_added(self) -> None:
        from nisdata import libraries
        tape = _fake_endf(92235, 9228, nsub=5)
        without = "".join(tape.splitlines(keepends=True)[1:])
        fixed = libraries.fix_tape(without, "JEFF-4.0 yields")
        self.assertTrue(fixed.startswith("JEFF-4.0 yields"))
        self.assertEqual(fixed.splitlines()[1:], without.splitlines())
        # a tape with its first line does not change
        self.assertEqual(libraries.fix_tape(tape, "x"), tape)

    def test_a_list_of_8_items_is_made_12(self) -> None:
        from nisdata import libraries
        head = _line(" 8.400000+4 3.000000+3          0          0          8          0", 128, 8, 457, 8)
        first = _line(" 1.000000+0 0.000000+0 8.420000+1 1.500000+0 0.000000+0 0.000000+0", 128, 8, 457, 9)
        second = _line(" 1.280000-1 2.100000-2", 128, 8, 457, 10)
        # a decay mode list has N2 > 0 and does not change
        modes = _line(" 0.000000+0 0.000000+0          0          0          8          1", 128, 8, 457, 5)
        tape = _fake_endf(2004, 128, nsub=4)
        fixed = libraries.fix_tape(tape + head + first + second, "x").splitlines(keepends=True)
        self.assertEqual(fixed[-3][44:55], "         12")
        self.assertEqual(fixed[-3][:44], head[:44])
        self.assertEqual(fixed[-2], first)
        self.assertEqual(fixed[-1], second[:22] + " 0.000000+0" * 4 + second[66:])
        self.assertEqual(len(fixed[-1]), len(second))
        # the same head with something in the second line after 2 items is an error
        bad = _line(" 1.280000-1 2.100000-2 1.000000+0", 128, 8, 457, 10)
        with self.assertRaises(ValueError):
            libraries.fix_tape(tape + head + first + bad, "x")
        self.assertEqual(libraries.fix_tape(modes, "x"), modes)

    def test_one_tape_from_a_zip(self) -> None:
        import zipfile
        from nisdata import libraries
        folder = Path(tempfile.mkdtemp())
        archive = folder / "decay.zip"
        with zipfile.ZipFile(archive, "w") as z:
            z.writestr("decay/dec-027_Co_060.endf", _fake_endf(27060, 2725, nsub=4))
            z.writestr("decay/dec-001_H_003.endf", _fake_endf(1003, 3, nsub=4))
            z.writestr("decay/CHANGELOG.txt", "not ENDF\n")
            z.writestr("decay/n-026_Fe_056.endf", _fake_endf(26056, 2631, nsub=10))
        tape = libraries.tape_from(archive, libraries.DECAY_SUBLIBRARY, "Decay")
        lines = tape.splitlines()
        self.assertTrue(lines[0].startswith("Decay"))
        mats = [endf.mat_mf_mt(line)[0] for line in lines if line[70:80] == " 1451    1"]
        self.assertEqual(mats, [3, 2725])     # H3 then Co60, and not the neutron file
        self.assertEqual(endf.mat_mf_mt(lines[-1])[0], -1)
        self.assertEqual(sum(1 for line in lines if endf.mat_mf_mt(line)[0] == -1), 1)

    def test_the_md5_is_checked(self) -> None:
        import hashlib
        from nisdata import libraries
        content = b"some data\n" * 1000
        headers = []
        address, stop = _serve(content, headers)
        self.addCleanup(stop)
        folder = Path(tempfile.mkdtemp())

        good = libraries.Source(f"{address}/a.txt", "a.txt", hashlib.md5(content).hexdigest())
        self.assertEqual(libraries.download(good, folder).read_bytes(), content)
        # the data bank refuses the default python name
        self.assertTrue(headers[0]["User-Agent"].startswith("nisdata"))

        bad = libraries.Source(f"{address}/b.txt", "b.txt", "0" * 32)
        with self.assertRaises(ValueError):
            libraries.download(bad, folder)
        self.assertFalse((folder / "b.txt").exists())
        self.assertTrue((folder / "b.txt.bad").exists())

    def test_make_a_decay_file(self) -> None:
        import hashlib
        import zipfile
        from nisdata import libraries
        buffer = Path(tempfile.mkdtemp()) / "d.zip"
        with zipfile.ZipFile(buffer, "w") as z:
            z.writestr("dec-001_H_003.endf", _fake_endf(1003, 3, nsub=4))
        content = buffer.read_bytes()
        address, stop = _serve(content)
        self.addCleanup(stop)
        lib = libraries.Library("test", "Test-1.0", "none",
                                decay=libraries.Source(f"{address}/d.zip", "d.zip", hashlib.md5(content).hexdigest()))
        raw, data = Path(tempfile.mkdtemp()), Path(tempfile.mkdtemp())
        path = libraries.make(lib, "decay", raw, data)
        self.assertEqual(path, data / "test.dec")
        self.assertTrue(path.read_text().startswith("Test-1.0 Radioactive Decay Data"))
        # it is not made again
        path.write_text("changed")
        libraries.make(lib, "decay", raw, data)
        self.assertEqual(path.read_text(), "changed")
        libraries.make(lib, "decay", raw, data, force=True)
        self.assertNotEqual(path.read_text(), "changed")
        with self.assertRaises(ValueError):
            libraries.make(lib, "nfy", raw, data)
        # a git lfs pointer file is made again
        path.write_text("version https://git-lfs.github.com/spec/v1\noid sha256:642e\nsize 67543551\n")
        self.assertFalse(libraries.is_there(path))
        libraries.make(lib, "decay", raw, data)
        self.assertTrue(path.read_text().startswith("Test-1.0"))

    def test_libraries_command(self) -> None:
        run = subprocess.run([sys.executable, str(PROCESS), "libraries"], capture_output=True, text=True, cwd=ROOT)
        self.assertEqual(run.returncode, 0)
        for key in ("jeff311", "jeff4", "endfbvii1", "endfbviii0"):
            self.assertIn(key, run.stdout)


def _numbers(values: list) -> list:
    """Lines of 6 ENDF numbers (the MAT, MF, MT are added later)."""
    from nisdata.mt5 import _number
    return ["".join(v if isinstance(v, str) else _number(v) for v in values[k:k + 6])
            for k in range(0, len(values), 6)]


def _tab1(c1, c2, l1, l2, xs: list, ys: list, law: int = 2) -> list:
    head = ["".join([*_numbers([c1, c2]), f"{l1:11d}{l2:11d}{1:11d}{len(xs):11d}"])]
    pairs = []
    for x, y in zip(xs, ys):
        pairs += [x, y]
    return head + [f"{len(xs):11d}{law:11d}"] + _numbers(pairs)


def _section(mat: int, mf: int, mt: int, rows: list) -> str:
    return "".join(_line(r, mat, mf, mt, i + 1) for i, r in enumerate(rows)) + _line("", mat, mf, 0, 99999)


def _material_with_mt5(mat: int = 2525) -> str:
    """Mn55 with MF3/MT5 and MF6/MT5: V52 (LAW=0), Cr54 in its 1st isomer (LAW=1, with angular data) and alpha."""
    head = "".join([*_numbers([25055.0, 54.466]), f"{0:11d}{0:11d}{0:11d}{0:11d}"])
    mf3 = _section(mat, 3, 5, [head] + _tab1(0.0, 0.0, 0, 0, [1e6, 1e7, 2e7], [0.0, 0.2, 0.4]))
    mf6_head = "".join([*_numbers([25055.0, 54.466]), f"{0:11d}{1:11d}{3:11d}{0:11d}"])
    v52 = _tab1(23052.0, 51.0, 0, 0, [1e6, 2e7], [1.0, 1.0])
    # LAW=1: a TAB2 with 1 energy, then a LIST of 3 numbers
    cr54 = _tab1(24054.0, 53.0, 1, 1, [1e6, 1e7, 2e7], [0.5, 0.5, 0.0])
    cr54 += [f"{'':22s}{1:11d}{2:11d}{1:11d}{1:11d}", f"{1:11d}{2:11d}"]
    cr54 += ["".join([*_numbers([0.0, 1e7]), f"{0:11d}{1:11d}{3:11d}{1:11d}"])] + _numbers([1.0, 2.0, 3.0])
    alpha = _tab1(2004.0, 3.97, 0, 0, [1e6, 2e7], [2.0, 2.0])
    mf6 = _section(mat, 6, 5, [mf6_head] + v52 + cr54 + alpha)
    mf10_other = _section(mat, 10, 37, ["".join([*_numbers([25055.0, 54.466]), f"{0:11d}{0:11d}{1:11d}{0:11d}"])]
                          + _tab1(0.0, 0.0, 25052, 0, [1e7, 2e7], [0.1, 0.1]))
    start = "".join(_fake_endf(25055, mat, nsub=10).splitlines(keepends=True)[:4])
    return (start + _line("", mat, 1, 0, 99999) + _line("", mat, 0, 0, 0) + mf3 + _line("", mat, 0, 0, 0)
            + mf6 + _line("", mat, 0, 0, 0) + mf10_other + _line("", mat, 0, 0, 0)
            + _line("", 0, 0, 0, 0) + _line("", -1, 0, 0, 0))


class ResidualProductionTests(unittest.TestCase):
    """The products of MT=5 from MF6 into MF10 (nisdata/mt5.py)."""

    def products(self, text: str) -> dict:
        from nisdata import mt5
        lines = [l for l in text.splitlines(keepends=True) if endf.mat_mf_mt(l)[1:] == (10, 5)]
        reader = mt5._Reader(lines)
        out = {}
        for _ in range(reader.cont()[4]):
            head, x, y, _, _ = reader.tab1()
            out[(head[2], head[3])] = dict(zip(x, y))
        return out

    def test_products(self) -> None:
        from nisdata import mt5
        folder = Path(tempfile.mkdtemp())
        source = folder / "in.endf"
        source.write_text(_material_with_mt5())
        self.assertEqual(mt5.add_residual_production(source, folder / "out.endf"), 3)
        out = (folder / "out.endf").read_text()
        products = self.products(out)
        self.assertEqual(set(products), {(23052, 0), (24054, 1), (2004, 0)})
        # sigma_5 times the yield, on the energies of both
        self.assertAlmostEqual(products[(23052, 0)][1e7], 0.2)
        self.assertAlmostEqual(products[(24054, 1)][1e7], 0.1)
        self.assertAlmostEqual(products[(24054, 1)][2e7], 0.0)
        self.assertEqual(sorted(products[(23052, 0)]), [1e6, 1e7, 2e7])
        self.assertAlmostEqual(products[(2004, 0)][2e7], 0.8)
        # MF10/MT5 is before MF10/MT37, in the same file
        keys = [endf.mat_mf_mt(l) for l in out.splitlines() if len(l) >= 75]
        sections = [k for i, k in enumerate(keys) if k[2] > 0 and (i == 0 or keys[i - 1] != k)]
        self.assertEqual([k for k in sections if k[1] == 10], [(2525, 10, 5), (2525, 10, 37)])
        self.assertEqual(sum(1 for k in keys if k == (2525, 0, 0)), 4)
        # nothing to add a second time
        self.assertEqual(mt5.add_residual_production(folder / "out.endf", folder / "again.endf"), 0)

    def test_a_step_stays_a_step(self) -> None:
        from nisdata import mt5
        # as after SIGMA1: the cross section starts with a step at its threshold (two points at 1e7)
        grid, values = mt5._product([1e7, 1e7, 2e7], [0.0, 0.2, 0.4], [3], [2], [1e-5, 2e7], [1.0, 1.0], [2], [2])
        self.assertEqual(grid, [1e-5, 1e7, 1e7, 2e7])
        self.assertEqual(values, [0.0, 0.0, 0.2, 0.4])
        # a step in the yield, and a histogram (INT=1)
        grid, values = mt5._product([1e6, 2e7], [1.0, 1.0], [2], [2], [1e6, 1e7, 1e7, 2e7], [0.5, 0.5, 1.0, 1.0],
                                    [4], [2])
        self.assertEqual(list(zip(grid, values)), [(1e6, 0.5), (1e7, 0.5), (1e7, 1.0), (2e7, 1.0)])
        self.assertEqual(mt5._value([1e6, 1e7, 2e7], [1.0, 3.0, 3.0], [3], [1], 5e6), 1.0)

    def test_numbers_fit_in_11_characters(self) -> None:
        from nisdata import mt5
        for v in (1.0, 1.18e7, 3.878602e-10, -5.536369e-26, 1.5e-120, -1.0e100, 2.0e-99, 123456789.0):
            text = mt5._number(v)
            self.assertEqual(len(text), 11, text)
            if abs(v) >= 1e-99:
                self.assertAlmostEqual(endf.endf_float(text) / v, 1.0, places=3)
        self.assertEqual(mt5._number(1.5e-120), " 0.000000+0")
        # no negative cross sections
        _, values = mt5._product([1e6, 2e7], [-1e-26, 1.0], [2], [2], [1e6, 2e7], [1.0, 1.0], [2], [2])
        self.assertEqual(values[0], 0.0)

    def test_new_mf10(self) -> None:
        from nisdata import mt5
        text = _material_with_mt5()
        lines = [l for l in text.splitlines(keepends=True) if endf.mat_mf_mt(l)[1] != 10]
        # one FEND line too many is left where MF10 was
        fends = [i for i, l in enumerate(lines) if endf.mat_mf_mt(l) == (2525, 0, 0)]
        del lines[fends[-1]]
        folder = Path(tempfile.mkdtemp())
        (folder / "in.endf").write_text("".join(lines))
        self.assertEqual(mt5.add_residual_production(folder / "in.endf", folder / "out.endf"), 3)
        out = (folder / "out.endf").read_text().splitlines()
        keys = [endf.mat_mf_mt(l) for l in out if len(l) >= 75]
        last10 = max(i for i, k in enumerate(keys) if k[1] == 10)
        self.assertEqual(keys[last10 + 1], (2525, 0, 0))
        self.assertEqual(keys[last10 + 2], (0, 0, 0))

    def test_no_mt5(self) -> None:
        from nisdata import mt5
        folder = Path(tempfile.mkdtemp())
        (folder / "in.endf").write_text(_fake_endf(26056, 2631))
        self.assertEqual(mt5.add_residual_production(folder / "in.endf", folder / "out.endf"), 0)
        self.assertEqual((folder / "out.endf").read_text(), _fake_endf(26056, 2631))


class VersionTests(unittest.TestCase):
    """Reaction data of an older version of the processing is made again."""

    def test_up_to_date(self) -> None:
        sys.path.insert(0, str(PROCESSING_DIR))
        import importlib
        process = importlib.import_module("process")
        folder = Path(tempfile.mkdtemp())
        new, old = folder / "new.reac", folder / "old.reac"
        new.write_text(endf.tpid_line(process._title("JEFF-4.0 Incident Neutron File")) + "x" * 2000)
        old.write_text(endf.tpid_line("JEFF-4.0 Incident Neutron File (1102 groups)") + "x" * 2000)
        self.assertTrue(process._up_to_date(new))
        self.assertFalse(process._up_to_date(old))
        self.assertFalse(process._up_to_date(folder / "none.reac"))

    def test_old_processed_files_are_moved(self) -> None:
        import importlib
        sys.path.insert(0, str(PROCESSING_DIR))
        process = importlib.import_module("process")
        folder = Path(tempfile.mkdtemp())
        (folder / "Fe56.asc").write_text("old")
        process._check_version(folder)
        self.assertFalse((folder / "Fe56.asc").exists())
        self.assertEqual((folder / "version-1" / "Fe56.asc").read_text(), "old")
        self.assertEqual((folder / "VERSION").read_text().strip(), str(process.VERSION))
        # a file of this version stays (a stopped run continues)
        (folder / "Co59.asc").write_text("new")
        process._check_version(folder)
        self.assertTrue((folder / "Co59.asc").exists())

    def test_allow_old(self) -> None:
        import importlib
        sys.path.insert(0, str(PROCESSING_DIR))
        folder = Path(tempfile.mkdtemp())
        old = folder / "tendl2025.reac"
        old.write_text(endf.tpid_line("TENDL-2025 Incident Neutron File (1102 groups)") + "x" * 2000)
        run = subprocess.run([sys.executable, str(PROCESS), "tendl", "--output", old, "--allow-old"],
                             capture_output=True, text=True, cwd=ROOT)
        self.assertEqual(run.returncode, 0)
        self.assertIn("it is used as it is", run.stderr)
        # a file of a bad version is never used
        bad = folder / "bad.reac"
        bad.write_text(endf.tpid_line("TENDL-2025 Incident Neutron File (1102 groups, nisdata 2)") + "x" * 2000)
        process = importlib.import_module("process")
        args = type("Args", (), {"force": False, "allow_old": True})()
        self.assertFalse(process._use_existing(bad, args))
        incomplete = folder / "incomplete.reac"
        incomplete.write_text(endf.tpid_line("TENDL-2025 Incident Neutron File (1102 groups, INCOMPLETE)") + "x" * 2000)
        self.assertFalse(process._use_existing(incomplete, args))
        # a file with a version, but nuclides failed in the last run of its folder (made by an older NIS)
        new = folder / "new.reac"
        new.write_text(endf.tpid_line(process._title("JEFF-4.0 Incident Neutron File")) + "x" * 2000)
        out = folder / "out"
        out.mkdir()
        self.assertTrue(process._use_existing(new, args, out))
        (out / "summary.json").write_text('{"files": 593, "failed": 1}')
        self.assertFalse(process._use_existing(new, args, out))
        self.assertTrue(process._use_existing(old, args))
