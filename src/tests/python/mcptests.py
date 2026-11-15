"""
Tests of the MCP server (nismcp): times, materials, API keys, the runs, the plots and the tools, in-process and
over HTTP.

Most tests use a small fake data set (H3 decay, and Co59 (n,g) Co60 with a fixed cross section), so they need no
nuclear data. The last test runs 1 kg of U235 in a 14 MeV flux with the real data, if it is there.
"""

import asyncio
import math
import socket
import tempfile
import threading
import time
import unittest
from pathlib import Path

import numpy as np

try:
    import mcp  # noqa: F401
    import periodictable  # noqa: F401
    _has_mcp = True
except ImportError:  # python 3.9, or without the mcp group
    _has_mcp = False

import pynis as nis

from .testdata import JEFF311_DECAY, JEFF311_NFY, JEFF4_REACTION, ROOT, has_data

if str(ROOT) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(ROOT))

AVOGADRO = 6.02214076e23
YEAR = 365.25 * 86400
H3_HALF_LIFE = 12.32 * YEAR
CO60_HALF_LIFE = 5.2714 * YEAR
CO59_CAPTURE_BARN = 37.2
EV = 1.602176634e-19


class FakeData:
    """The interface of nismcp.data.NuclearData, with a few nuclides: H3 -> He3, Co60 -> Ni60, and a reaction
    Co59 (n,g) Co60 of 37.2 b at all energies."""

    def __init__(self, files=None):
        import scipy.sparse as sp
        self.nis = nis
        self.files = files
        size = len(nis.core.vector_from_inventory({}))
        index = {name: int(np.argmax(nis.core.vector_from_inventory({name: 1.0})))
                 for name in ("H3", "He3", "Co59", "Co60", "Ni60")}
        self.index = index
        lam = np.zeros(size)
        lam[index["H3"]] = math.log(2) / H3_HALF_LIFE
        lam[index["Co60"]] = math.log(2) / CO60_HALF_LIFE
        # the matrix has the daughters of a parent in its row, and -lambda on the diagonal
        rows, cols, values = [], [], []
        for parent, daughter in (("H3", "He3"), ("Co60", "Ni60")):
            p, d = index[parent], index[daughter]
            rows += [p, p]
            cols += [p, d]
            values += [-lam[p], lam[p]]
        self.decay_mat = sp.csc_matrix((values, (rows, cols)), shape=(size, size))
        energy = nis.core.DecayEnergyData()
        energy.decay_constants = lam
        light, em, heavy = np.zeros(size), np.zeros(size), np.zeros(size)
        light[index["H3"]] = 5.7e3
        light[index["Co60"]] = 9.6e4
        em[index["Co60"]] = 2.5e6
        energy.light_particle, energy.electromagnetic, energy.heavy_particle = light, em, heavy
        self.energy_data = energy
        self.has_reactions = True
        self.has_fission_yields = False

    def reaction_matrix(self, bounds, values):
        import scipy.sparse as sp
        size = self.decay_mat.shape[0]
        rate = CO59_CAPTURE_BARN * 1e-24 * float(sum(values))
        p, d = self.index["Co59"], self.index["Co60"]
        return sp.csc_matrix(([-rate, rate], ([p, p], [p, d])), shape=(size, size))

    def describe(self):
        return {"decay": "fake"}


def _loader():
    from nismcp.data import DataFiles, DataLoader
    loader = DataLoader(DataFiles(decay=Path("fake")), factory=FakeData).start()
    loader.get()
    return loader


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@unittest.skipUnless(_has_mcp, "the mcp group is not installed (uv sync --group mcp)")
class MCPUnitsTest(unittest.TestCase):

    def test_parse_time(self) -> None:
        from nismcp.units import parse_time
        self.assertEqual(parse_time("1 d"), 86400)
        self.assertEqual(parse_time("30min"), 1800)
        self.assertEqual(parse_time("2.5 h"), 9000)
        self.assertEqual(parse_time("1e3 s"), 1000)
        self.assertEqual(parse_time(60), 60)
        self.assertEqual(parse_time("5 m"), 300)
        self.assertAlmostEqual(parse_time("1 y"), YEAR)
        self.assertEqual(parse_time("10"), 10)
        for bad in ("abc", "1 parsec", "-1 d", True):
            with self.assertRaises(ValueError):
                parse_time(bad)

    def test_format_time(self) -> None:
        from nismcp.units import format_time
        self.assertEqual(format_time(86400), "1 d")
        self.assertEqual(format_time(300), "5 min")
        self.assertEqual(format_time(0), "0 s")
        self.assertEqual(format_time(0.5), "0.5 s")


@unittest.skipUnless(_has_mcp, "the mcp group is not installed (uv sync --group mcp)")
class MCPMaterialsTest(unittest.TestCase):

    def setUp(self) -> None:
        zais = nis.data.ALL_ZAIS()
        self.known = {nis.zai.get_name_from(z): i for i, z in enumerate(zais)}

    def test_one_nuclide(self) -> None:
        from nismcp.materials import atoms_from_mass
        atoms, notes = atoms_from_mass(1.0, {"U235": 1.0}, self.known)
        self.assertEqual(list(atoms), ["U235"])
        self.assertAlmostEqual(atoms["U235"] / (1000 / 235.0439 * AVOGADRO), 1.0, places=5)
        self.assertEqual(notes, [])

    def test_names(self) -> None:
        from nismcp.materials import canonical_nuclide
        self.assertEqual(canonical_nuclide("u235", self.known), "U235")
        self.assertEqual(canonical_nuclide("Tc-99m", self.known), "Tc99m")
        self.assertEqual(canonical_nuclide("Am242m1", self.known), "Am242m")
        for bad in ("U", "Xx12", "U999", "235U"):
            with self.assertRaises(ValueError):
                canonical_nuclide(bad, self.known)

    def test_natural_elements(self) -> None:
        from nismcp.materials import atoms_from_mass, natural_abundances
        atoms, _ = atoms_from_mass(1.0, {"Fe": 1.0}, self.known)
        self.assertEqual(set(atoms), {"Fe54", "Fe56", "Fe57", "Fe58"})
        # 1 kg of iron is 1000 / 55.845 mol
        self.assertAlmostEqual(sum(atoms.values()) / (1000 / 55.845 * AVOGADRO), 1.0, places=3)
        self.assertAlmostEqual(atoms["Fe56"] / sum(atoms.values()), 0.91754, places=4)
        u = natural_abundances("U")
        self.assertAlmostEqual(u[235], 0.007204, places=6)
        self.assertAlmostEqual(sum(u.values()), 1.0)
        with self.assertRaises(ValueError):
            natural_abundances("Pu")
        # the natural Ta180 is the metastable state (almost stable), not the 8 h ground state
        atoms, _ = atoms_from_mass(1.0, {"Ta": 1.0}, self.known)
        self.assertEqual(set(atoms), {"Ta180m", "Ta181"})
        self.assertAlmostEqual(atoms["Ta180m"] / sum(atoms.values()), 1.201e-4, places=7)
        # natural Th is all Th232
        atoms, _ = atoms_from_mass(1.0, {"Th": 1.0}, self.known)
        self.assertEqual(set(atoms), {"Th232"})

    def test_mixture_and_normalisation(self) -> None:
        from nismcp.materials import atoms_from_mass
        atoms, notes = atoms_from_mass(2.0, {"Fe": 70, "Cr": 18, "Ni": 12}, self.known)
        self.assertEqual(notes, [])
        cr = sum(v for k, v in atoms.items() if k.startswith("Cr"))
        self.assertAlmostEqual(cr / (2000 * 0.18 / 51.9961 * AVOGADRO), 1.0, places=3)
        _, notes = atoms_from_mass(1.0, {"Fe": 0.5, "Cr": 0.2}, self.known)
        self.assertEqual(len(notes), 1)
        for bad in ({}, {"Fe": -1}, {"Fe": 0}):
            with self.assertRaises(ValueError):
                atoms_from_mass(1.0, bad, self.known)
        with self.assertRaises(ValueError):
            atoms_from_mass(0.0, {"Fe": 1}, self.known)


@unittest.skipUnless(_has_mcp, "the mcp group is not installed (uv sync --group mcp)")
class MCPAuthTest(unittest.TestCase):

    def setUp(self) -> None:
        self.folder = tempfile.TemporaryDirectory()
        self.keys = Path(self.folder.name) / "keys.txt"

    def tearDown(self) -> None:
        self.folder.cleanup()

    def test_keys_file(self) -> None:
        from nismcp.auth import KeyStore, add_key, hash_key
        alice = add_key(self.keys, "alice")
        bob = add_key(self.keys, "bob")
        self.assertNotEqual(alice, bob)
        text = self.keys.read_text()
        # only the hashes are in the file
        self.assertNotIn(alice, text)
        self.assertIn(hash_key(alice), text)
        self.assertEqual(self.keys.stat().st_mode & 0o777, 0o600)
        with self.assertRaises(ValueError):
            add_key(self.keys, "alice")
        with self.assertRaises(ValueError):
            add_key(self.keys, "bad name")

        store = KeyStore(self.keys)
        self.assertEqual(store.user(alice), "alice")
        self.assertEqual(store.user(bob), "bob")
        self.assertIsNone(store.user("nis_wrong"))

        # remove bob: the store reads the file again
        time.sleep(0.01)
        self.keys.write_text("\n".join(l for l in text.splitlines() if not l.startswith("bob")) + "\n# a comment\n")
        self.assertIsNone(store.user(bob))
        self.assertEqual(store.user(alice), "alice")

        # fail closed: a bad edit, or no file, and no key works (else a user that was just removed keeps access)
        good = self.keys.read_text()
        time.sleep(0.01)
        self.keys.write_text("this is not a key line\n")
        self.assertIsNone(store.user(alice))
        self.keys.unlink()
        self.assertIsNone(store.user(alice))
        self.keys.write_text(good)
        self.assertEqual(store.user(alice), "alice")

    def test_bad_file(self) -> None:
        from nismcp.auth import parse_keys
        with self.assertRaises(ValueError):
            parse_keys("alice 1234")
        with self.assertRaises(ValueError):
            parse_keys(f"alice {'a' * 64}\nbob {'a' * 64}")

    def test_headers(self) -> None:
        from nismcp.auth import key_from_headers
        self.assertEqual(key_from_headers({"authorization": "Bearer abc"}), "abc")
        self.assertEqual(key_from_headers({"authorization": "bearer abc "}), "abc")
        self.assertEqual(key_from_headers({"x-api-key": "abc"}), "abc")
        self.assertIsNone(key_from_headers({"authorization": "Basic abc"}))
        self.assertIsNone(key_from_headers({}))

    def test_middleware(self) -> None:
        from nismcp.auth import ApiKeyMiddleware, KeyStore, add_key
        key = add_key(self.keys, "alice")
        seen = {}

        async def app(scope, receive, send):
            seen["user"] = scope["state"]["nis_user"]
            await send({"type": "http.response.start", "status": 200, "headers": []})
            await send({"type": "http.response.body", "body": b"ok"})

        middleware = ApiKeyMiddleware(app, KeyStore(self.keys))

        def call(headers, path="/mcp"):
            sent = []

            async def send(message):
                sent.append(message)

            async def receive():
                return {"type": "http.request", "body": b""}

            scope = {"type": "http", "path": path, "headers": [(k.encode(), v.encode()) for k, v in headers.items()]}
            asyncio.run(middleware(scope, receive, send))
            return sent[0]["status"]

        self.assertEqual(call({}), 401)
        self.assertEqual(call({"Authorization": "Bearer nis_wrong"}), 401)
        self.assertEqual(call({"Authorization": f"Bearer {key}"}), 200)
        self.assertEqual(seen["user"], "alice")
        self.assertEqual(call({"X-API-Key": key}), 200)


@unittest.skipUnless(_has_mcp, "the mcp group is not installed (uv sync --group mcp)")
class MCPModelTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        from nismcp.model import Runner
        cls.loader = _loader()
        cls.data, cls.tables = cls.loader.get()
        cls.runner = Runner(cls.loader, max_runs=3)

    def test_steps(self) -> None:
        from nismcp.model import DEFAULT_FLUX, steps_from_input
        steps = steps_from_input("1 d", ["1 h", "1 d", "1 h"], None, None)
        self.assertEqual([(s.duration, s.flux) for s in steps],
                         [(86400, DEFAULT_FLUX), (3600, 0.0), (82800, 0.0)])
        steps = steps_from_input(None, None, [{"duration": "1 h", "flux": 1e12}, {"duration": 60, "flux": 0},
                                              {"duration": "1 h"}], 5e10)
        self.assertEqual([s.flux for s in steps], [1e12, 0.0, 5e10])
        with self.assertRaises(ValueError):
            steps_from_input(None, None, None, None)
        with self.assertRaises(ValueError):
            steps_from_input("1 d", None, [{"duration": "1 h"}], None)
        with self.assertRaises(ValueError):
            steps_from_input("1 d", None, None, -1)
        with self.assertRaises(ValueError):
            steps_from_input("20000 y", None, None, None)

    def test_spectrum(self) -> None:
        from nismcp.model import spectrum_from_input
        bounds, values, text = spectrum_from_input("14 MeV")
        self.assertTrue(bounds[0] < 14.1e6 < bounds[1])
        self.assertTrue(text.startswith("14MeV"))
        self.assertTrue(spectrum_from_input(None)[2].startswith("14MeV"))
        self.assertTrue(spectrum_from_input("Thermal")[2].startswith("thermal"))
        bounds, values, _ = spectrum_from_input({"bounds_ev": [1, 10, 100], "flux": [1, 2]})
        self.assertEqual(values, [1.0, 2.0])
        for bad in ("fission", {"bounds_ev": [10, 1], "flux": [1]}, {"bounds_ev": [1, 10], "flux": [1, 2]}, 3):
            with self.assertRaises(ValueError):
                spectrum_from_input(bad)

    def test_decay_half_life(self) -> None:
        from nismcp.model import Step, spectrum_from_input
        result = self.runner.run("t", {"H3": 1e20}, [Step(H3_HALF_LIFE, 0.0)], spectrum_from_input(None), {}, [])
        atoms = result.points[-1].atoms
        self.assertAlmostEqual(atoms[self.tables.index["H3"]] / 1e20, 0.5, places=6)
        self.assertAlmostEqual(atoms[self.tables.index["He3"]] / 1e20, 0.5, places=6)
        # the sub-steps go up in time, the last is the end of the step
        times = [p.time for p in result.points]
        self.assertEqual(times, sorted(times))
        self.assertAlmostEqual(times[-1], H3_HALF_LIFE)

    def test_activation_bateman(self) -> None:
        """Co59 (n,g) Co60 in a flux, then cooling: the analytic solution."""
        from nismcp.model import Step, spectrum_from_input, totals
        flux, t_irr, t_cool, n0 = 1e14, 30 * 86400, YEAR, 1e24
        steps = [Step(t_irr, flux), Step(t_cool, 0.0)]
        result = self.runner.run("t", {"Co59": n0}, steps, spectrum_from_input("thermal"), {}, [])
        rate = CO59_CAPTURE_BARN * 1e-24 * flux
        lam = math.log(2) / CO60_HALF_LIFE
        co60_end = n0 * rate / (lam - rate) * (math.exp(-rate * t_irr) - math.exp(-lam * t_irr))
        end = result.point_at(t_irr)
        self.assertAlmostEqual(end.atoms[self.tables.index["Co60"]] / co60_end, 1.0, places=6)
        last = result.points[-1]
        self.assertAlmostEqual(last.atoms[self.tables.index["Co60"]] / (co60_end * math.exp(-lam * t_cool)), 1.0,
                               places=6)
        activity = totals(self.tables, last.atoms)["activity_Bq"]
        self.assertAlmostEqual(activity / (lam * co60_end * math.exp(-lam * t_cool)), 1.0, places=3)
        self.assertAlmostEqual(result.irradiation_end, t_irr)
        with self.assertRaises(ValueError):
            result.point_at(1.2345)

    def test_short_cooling_after_long_irradiation(self) -> None:
        """After 1 y of irradiation, a cooling time of 1 min is not the point at 44 s (the nearest sub-step)."""
        from nismcp.model import spectrum_from_input, steps_from_input
        steps = steps_from_input("1 y", ["1 h"], None, 1e14)
        result = self.runner.run("t", {"Co59": 1e24}, steps, spectrum_from_input("thermal"), {}, [])
        end = result.irradiation_end
        self.assertEqual(result.point_at(0, after=end).time, end)
        self.assertAlmostEqual(result.point_at(3600, after=end).time - end, 3600, places=3)
        for t in (60, 10, 1):
            with self.assertRaises(ValueError):
                result.point_at(t, after=end)

    def test_point_storage(self) -> None:
        import numpy as np
        from nismcp.model import Point
        atoms = np.zeros(4131)
        atoms[[3, 100, 4000]] = [1.5, 2e20, 3e-5]
        point = Point(1.0, 0, True, atoms)
        np.testing.assert_array_equal(point.atoms, atoms)
        self.assertLess(point.nbytes, 1000)

    def test_memory_of_all_runs(self) -> None:
        """At most max_results_bytes for all the users: the oldest runs go first, the last run stays."""
        from nismcp.model import Runner, Step, spectrum_from_input
        spec = spectrum_from_input(None)
        runner = Runner(self.loader, max_results_bytes=1)
        first = runner.run("alice", {"H3": 1.0}, [Step(1.0, 0.0)], spec, {}, []).run_id
        last = runner.run("bob", {"H3": 1.0}, [Step(1.0, 0.0)], spec, {}, []).run_id
        self.assertEqual(runner.list("alice"), [])
        self.assertEqual(runner.get("bob", last).run_id, last)
        with self.assertRaises(ValueError):
            runner.get("alice", first)

    def test_runs_per_user(self) -> None:
        from nismcp.model import Runner, Step, spectrum_from_input
        runner = Runner(self.loader, max_runs=2)
        spec = spectrum_from_input(None)
        ids = [runner.run("alice", {"H3": 1.0}, [Step(1.0, 0.0)], spec, {}, []).run_id for _ in range(3)]
        bob = runner.run("bob", {"H3": 1.0}, [Step(1.0, 0.0)], spec, {}, []).run_id
        self.assertEqual([r.run_id for r in runner.list("alice")], ids[1:])
        with self.assertRaises(ValueError):
            runner.get("alice", ids[0])
        with self.assertRaises(ValueError):
            runner.get("alice", bob)
        self.assertEqual(runner.get("bob", bob).run_id, bob)

    def test_plot(self) -> None:
        from nismcp.model import Step, spectrum_from_input
        from nismcp.plot import plot_png
        result = self.runner.run("t", {"Co59": 1e24, "H3": 1e20}, [Step(86400, 1e14), Step(YEAR, 0.0)],
                                 spectrum_from_input(None), {}, [])
        png = plot_png(result, self.tables, ["Co60", "H3"], "activity")
        self.assertTrue(png.startswith(b"\x89PNG"))
        png = plot_png(result, self.tables, ["Co60"], "decay_heat", log_time=False, log_value=False)
        self.assertTrue(png.startswith(b"\x89PNG"))
        with self.assertRaises(ValueError):
            plot_png(result, self.tables, ["Ni60"], "activity")
        with self.assertRaises(ValueError):
            plot_png(result, self.tables, ["Co60"] * 9, "activity")

    def test_plots_in_threads(self) -> None:
        """The tools run in threads, so plots can be made at the same time."""
        from concurrent.futures import ThreadPoolExecutor
        from nismcp.model import Step, spectrum_from_input
        from nismcp.plot import plot_png
        result = self.runner.run("t", {"Co59": 1e24, "H3": 1e20}, [Step(86400, 1e14), Step(YEAR, 0.0)],
                                 spectrum_from_input(None), {}, [])
        with ThreadPoolExecutor(8) as pool:
            pngs = list(pool.map(lambda i: plot_png(result, self.tables, ["Co60", "H3"][: 1 + i % 2]), range(16)))
        self.assertTrue(all(png.startswith(b"\x89PNG") for png in pngs))


def _text(result) -> str:
    return "\n".join(c.text for c in result.content if c.type == "text")


@unittest.skipUnless(_has_mcp, "the mcp group is not installed (uv sync --group mcp)")
class MCPToolsTest(unittest.TestCase):
    """The tools through an MCP client, in-process (as over stdio: no keys)."""

    @classmethod
    def setUpClass(cls) -> None:
        from nismcp.model import Runner
        from nismcp.server import build_server
        cls.loader = _loader()
        cls.server = build_server(cls.loader, Runner(cls.loader))

    def call(self, tool, arguments):
        from mcp import Client

        async def run():
            async with Client(self.server) as client:
                return await client.call_tool(tool, arguments)
        return asyncio.run(run())

    def test_list_tools(self) -> None:
        from mcp import Client

        async def run():
            async with Client(self.server) as client:
                return await client.list_tools()
        names = {t.name for t in asyncio.run(run()).tools}
        self.assertEqual(names, {"nuclear_data_info", "irradiate", "get_results", "plot_nuclides", "list_runs"})

    def test_irradiate_results_plot(self) -> None:
        result = self.call("irradiate", {"material": {"Co59": 1.0}, "mass_kg": 1.0, "irradiation_time": "30 d",
                                         "cooling_times": ["1 d", "1 y"], "spectrum": "thermal"})
        self.assertFalse(result.is_error, _text(result))
        out = result.structured_content
        self.assertEqual(len(out["results_at_step_ends"]), 4)
        last = out["results_at_step_ends"][-1]
        self.assertEqual(last["cooling_time"], "1 y")
        self.assertEqual(last["top_by_activity"][0]["nuclide"], "Co60")
        self.assertTrue(any("default" in a for a in out["assumptions"]))
        self.assertTrue(any("0-D" in a for a in out["assumptions"]))
        run_id = out["run_id"]

        result = self.call("get_results", {"run_id": run_id, "cooling_time": "1 d", "top_n": 3})
        self.assertFalse(result.is_error, _text(result))
        self.assertEqual(result.structured_content["nuclides"][0]["nuclide"], "Co60")
        self.assertEqual(result.structured_content["nuclides"][0]["half_life"], "5.271 y")

        result = self.call("get_results", {"run_id": run_id, "nuclides": ["co60", "Co59"], "quantity": "atoms"})
        self.assertFalse(result.is_error, _text(result))
        self.assertEqual([r["nuclide"] for r in result.structured_content["nuclides"]], ["Co60", "Co59"])

        result = self.call("plot_nuclides", {"run_id": run_id, "top_n": 3})
        self.assertFalse(result.is_error, _text(result))
        self.assertEqual([c.type for c in result.content], ["image", "text"])
        self.assertEqual(result.content[0].mime_type, "image/png")
        self.assertIn("Co60", _text(result))

        result = self.call("list_runs", {})
        self.assertIn(run_id, _text(result))

    def test_errors(self) -> None:
        for arguments in ({"material": {"Co59": 1}, "irradiation_time": "1 d"},
                          {"material": {"Xx99": 1}, "mass_kg": 1, "irradiation_time": "1 d"},
                          {"material": {"Co59": 1}, "mass_kg": 1, "irradiation_time": "1 parsec"},
                          {"atoms": {"Co59": 1}, "mass_kg": 1, "irradiation_time": "1 d"}):
            result = self.call("irradiate", arguments)
            self.assertTrue(result.is_error, arguments)
        result = self.call("get_results", {"run_id": "nope"})
        self.assertTrue(result.is_error)
        self.assertIn("no run", _text(result))

    def test_info(self) -> None:
        result = self.call("nuclear_data_info", {})
        self.assertFalse(result.is_error, _text(result))
        self.assertEqual(result.structured_content["status"], "ready")
        self.assertIn("14MeV", result.structured_content["spectra"])


@unittest.skipUnless(_has_mcp, "the mcp group is not installed (uv sync --group mcp)")
class MCPHttpTest(unittest.TestCase):
    """The HTTP server with API keys, with uvicorn, and an MCP client over HTTP."""

    @classmethod
    def setUpClass(cls) -> None:
        import uvicorn
        from nismcp.auth import ApiKeyMiddleware, KeyStore, add_key
        from nismcp.model import Runner
        from nismcp.server import build_server
        cls.folder = tempfile.TemporaryDirectory()
        keys = Path(cls.folder.name) / "keys.txt"
        cls.alice = add_key(keys, "alice")
        cls.bob = add_key(keys, "bob")
        store = KeyStore(keys)
        loader = _loader()
        server = build_server(loader, Runner(loader), key_store=store)
        app = ApiKeyMiddleware(server.streamable_http_app(), store)
        cls.port = _free_port()
        cls.uvicorn = uvicorn.Server(uvicorn.Config(app, host="127.0.0.1", port=cls.port, log_level="warning"))
        cls.thread = threading.Thread(target=cls.uvicorn.run, daemon=True)
        cls.thread.start()
        for _ in range(100):
            if cls.uvicorn.started:
                break
            time.sleep(0.05)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.uvicorn.should_exit = True
        cls.thread.join(5)
        cls.folder.cleanup()

    def url(self, path="/mcp") -> str:
        return f"http://127.0.0.1:{self.port}{path}"

    def call(self, key, tool, arguments):
        from mcp import Client
        from mcp.client.streamable_http import streamable_http_client
        from mcp.shared._httpx_utils import create_mcp_http_client

        async def run():
            http = create_mcp_http_client(headers={"Authorization": f"Bearer {key}"})
            async with http:
                async with Client(streamable_http_client(self.url(), http_client=http)) as client:
                    return await client.call_tool(tool, arguments)
        return asyncio.run(run())

    def test_no_key(self) -> None:
        import httpx2
        body = {"jsonrpc": "2.0", "id": 1, "method": "tools/list"}
        response = httpx2.post(self.url(), json=body)
        self.assertEqual(response.status_code, 401)
        response = httpx2.post(self.url(), json=body, headers={"Authorization": "Bearer nis_wrong"})
        self.assertEqual(response.status_code, 401)
        # the health check needs no key
        response = httpx2.get(self.url("/health"))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["status"], "ready")

    def test_runs_of_each_user(self) -> None:
        result = self.call(self.alice, "irradiate", {"atoms": {"H3": 1e20}, "irradiation_time": "1 y",
                                                     "flux_n_cm2_s": 0})
        self.assertFalse(result.is_error, _text(result))
        run_id = result.structured_content["run_id"]
        result = self.call(self.alice, "get_results", {"run_id": run_id})
        self.assertFalse(result.is_error, _text(result))
        self.assertEqual(result.structured_content["nuclides"][0]["nuclide"], "H3")
        # bob cannot see the run of alice
        result = self.call(self.bob, "get_results", {"run_id": run_id})
        self.assertTrue(result.is_error)


class FakeGitHub:
    """GitHub for the tests: the code "<login>-code" logs in the user <login>. The account ids are in ids."""

    def __init__(self):
        from nismcp.oauth import GitHub
        self._real = GitHub("fake-id", "fake-secret")
        self.redirects = []
        self.ids = {}

    def login_url(self, state, redirect_uri):
        self.redirects.append(redirect_uri)
        return self._real.login_url(state, redirect_uri)

    def user_login(self, code, redirect_uri):
        if not code.endswith("-code"):
            raise ValueError("bad code")
        login = code[:-len("-code")]
        return login, self.ids.setdefault(login.lower(), 1000 + len(self.ids))


@unittest.skipUnless(_has_mcp, "the mcp group is not installed (uv sync --group mcp)")
class MCPAllowListTest(unittest.TestCase):

    def test_allow_list(self) -> None:
        from nismcp.oauth import AllowList
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "users.txt"
            path.write_text("# the users\nAlice\n  bob  # a comment\n\n")
            allow = AllowList(path)
            self.assertTrue(allow.allows("alice"))
            self.assertTrue(allow.allows("ALICE"))
            self.assertTrue(allow.allows("bob"))
            self.assertFalse(allow.allows("mallory"))
            time.sleep(0.01)
            path.write_text("bob\n")
            self.assertFalse(allow.allows("alice"))
            self.assertTrue(allow.allows("bob"))
            # a login with its account id: only that account
            time.sleep(0.01)
            path.write_text("bob 42\n")
            self.assertTrue(allow.allows("bob", 42))
            self.assertFalse(allow.allows("bob", 43))
            # fail closed: a bad line, or no file, lets no one in
            time.sleep(0.01)
            path.write_text("bob 42 extra\n")
            self.assertFalse(allow.allows("bob", 42))
            path.unlink()
            self.assertFalse(allow.allows("bob", 42))
            path.write_text("bob\n")
            self.assertTrue(allow.allows("bob", 42))


@unittest.skipUnless(_has_mcp, "the mcp group is not installed (uv sync --group mcp)")
class MCPOAuthTest(unittest.TestCase):
    """OAuth as claude.ai does it: metadata, registration, authorize (PKCE), the GitHub login (fake), token,
    refresh, then the tools with the token. The API keys work on the same server."""

    REDIRECT = "https://claude.ai/api/mcp/auth_callback"

    @classmethod
    def setUpClass(cls) -> None:
        import uvicorn
        from nismcp.auth import KeyStore, add_key
        from nismcp.model import Runner
        from nismcp.oauth import AllowList, ApiKeyHeader, GitHubOAuthProvider, auth_settings
        from nismcp.server import build_server
        cls.folder = tempfile.TemporaryDirectory()
        folder = Path(cls.folder.name)
        keys = folder / "keys.txt"
        cls.key = add_key(keys, "scripts")
        cls.users = folder / "users.txt"
        cls.users.write_text("alice\nbob\n")
        cls.state_file = folder / "oauth.json"
        cls.port = _free_port()
        cls.public = f"http://127.0.0.1:{cls.port}"
        cls.github = FakeGitHub()
        store = KeyStore(keys)
        cls.provider = GitHubOAuthProvider(cls.github, AllowList(cls.users), cls.public, state_file=cls.state_file,
                                           key_store=store)
        loader = _loader()
        server = build_server(loader, Runner(loader), key_store=store,
                              oauth=(cls.provider, auth_settings(cls.public)))
        app = ApiKeyHeader(server.streamable_http_app())
        cls.uvicorn = uvicorn.Server(uvicorn.Config(app, host="127.0.0.1", port=cls.port, log_level="warning"))
        cls.thread = threading.Thread(target=cls.uvicorn.run, daemon=True)
        cls.thread.start()
        for _ in range(100):
            if cls.uvicorn.started:
                break
            time.sleep(0.05)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.uvicorn.should_exit = True
        cls.thread.join(5)
        cls.folder.cleanup()

    def http(self):
        import httpx2
        return httpx2.Client(base_url=self.public, follow_redirects=False)

    def register(self, http) -> str:
        response = http.post("/register", json={
            "redirect_uris": [self.REDIRECT], "token_endpoint_auth_method": "none", "client_name": "test",
            "grant_types": ["authorization_code", "refresh_token"], "response_types": ["code"]})
        self.assertEqual(response.status_code, 201, response.text)
        return response.json()["client_id"]

    def login(self, http, client_id: str, login: str):
        """The browser part: /authorize, GitHub, the callback. The redirect to the client, as a dict."""
        import base64
        import hashlib
        import secrets
        from urllib.parse import parse_qs, urlparse
        verifier = secrets.token_urlsafe(48)
        challenge = base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest()).rstrip(b"=").decode()
        response = http.get("/authorize", params={
            "response_type": "code", "client_id": client_id, "redirect_uri": self.REDIRECT, "state": "s1",
            "code_challenge": challenge, "code_challenge_method": "S256", "resource": f"{self.public}/mcp"})
        self.assertEqual(response.status_code, 302, response.text)
        consent = urlparse(response.headers["location"])
        self.assertEqual(consent.path, "/oauth/consent")
        response = self.consent(http, f"{consent.path}?{consent.query}", "allow")
        self.assertEqual(response.status_code, 303, response.text)
        github = urlparse(response.headers["location"])
        self.assertEqual(github.netloc, "github.com")
        query = parse_qs(github.query)
        self.assertEqual(query["redirect_uri"], [f"{self.public}/oauth/github/callback"])
        response = http.get("/oauth/github/callback", params={"code": f"{login}-code", "state": query["state"][0]})
        self.assertEqual(response.status_code, 302, response.text)
        back = urlparse(response.headers["location"])
        self.assertEqual(f"{back.scheme}://{back.netloc}{back.path}", self.REDIRECT)
        result = {k: v[0] for k, v in parse_qs(back.query).items()}
        result["verifier"] = verifier
        return result

    def consent(self, http, url: str, action: str, page_check=None):
        """Open the consent page, and click a button."""
        import re
        page = http.get(url)
        self.assertEqual(page.status_code, 200, page.text)
        self.assertEqual(page.headers["x-frame-options"], "DENY")
        self.assertIn("frame-ancestors 'none'", page.headers["content-security-policy"])
        if page_check:
            page_check(page.text)
        fields = dict(re.findall(r"name=(state|csrf) value='([^']*)'", page.text))
        return http.post("/oauth/consent", data={**fields, "action": action})

    def start(self, http, client_id: str, redirect=None) -> str:
        """/authorize: the URL of the consent page."""
        from urllib.parse import urlparse
        response = http.get("/authorize", params={
            "response_type": "code", "client_id": client_id, "redirect_uri": redirect or self.REDIRECT,
            "state": "s1", "code_challenge": "x" * 43, "code_challenge_method": "S256"})
        self.assertEqual(response.status_code, 302, response.text)
        location = urlparse(response.headers["location"])
        return f"{location.path}?{location.query}"

    def tokens(self, http, client_id: str, login: str) -> dict:
        back = self.login(http, client_id, login)
        self.assertEqual(back["state"], "s1")
        response = http.post("/token", data={
            "grant_type": "authorization_code", "code": back["code"], "redirect_uri": self.REDIRECT,
            "client_id": client_id, "code_verifier": back["verifier"], "resource": f"{self.public}/mcp"})
        self.assertEqual(response.status_code, 200, response.text)
        return response.json()

    def call(self, token, tool, arguments, header="authorization"):
        from mcp import Client
        from mcp.client.streamable_http import streamable_http_client
        from mcp.shared._httpx_utils import create_mcp_http_client
        headers = {"Authorization": f"Bearer {token}"} if header == "authorization" else {"X-API-Key": token}

        async def run():
            http = create_mcp_http_client(headers=headers)
            async with http:
                async with Client(streamable_http_client(f"{self.public}/mcp", http_client=http)) as client:
                    return await client.call_tool(tool, arguments)
        return asyncio.run(run())

    def test_metadata(self) -> None:
        with self.http() as http:
            response = http.post("/mcp", json={"jsonrpc": "2.0", "id": 1, "method": "tools/list"})
            self.assertEqual(response.status_code, 401)
            self.assertIn("resource_metadata=", response.headers["www-authenticate"])
            resource = http.get("/.well-known/oauth-protected-resource/mcp").json()
            self.assertEqual(resource["resource"], f"{self.public}/mcp")
            self.assertEqual(resource["authorization_servers"], [self.public])
            server = http.get("/.well-known/oauth-authorization-server").json()
            self.assertEqual(server["registration_endpoint"], f"{self.public}/register")
            self.assertIn("S256", server["code_challenge_methods_supported"])
            self.assertEqual(http.get("/health").status_code, 200)

    def test_login_and_tools(self) -> None:
        with self.http() as http:
            client_id = self.register(http)
            alice = self.tokens(http, client_id, "alice")
            bob = self.tokens(http, client_id, "Bob")
        result = self.call(alice["access_token"], "irradiate",
                           {"atoms": {"H3": 1e20}, "irradiation_time": "1 y", "flux_n_cm2_s": 0})
        self.assertFalse(result.is_error, _text(result))
        run_id = result.structured_content["run_id"]
        subject = self.provider.access[next(iter(self.provider.access))].subject
        # the user is the GitHub account id, not the login
        self.assertTrue(subject.startswith("github:") and subject[7:].isdigit(), subject)
        result = self.call(alice["access_token"], "get_results", {"run_id": run_id})
        self.assertFalse(result.is_error, _text(result))
        # bob cannot see the run of alice
        result = self.call(bob["access_token"], "get_results", {"run_id": run_id})
        self.assertTrue(result.is_error)
        # a wrong token
        with self.http() as http:
            response = http.post("/mcp", json={"jsonrpc": "2.0", "id": 1, "method": "tools/list"},
                                 headers={"Authorization": "Bearer not-a-token"})
            self.assertEqual(response.status_code, 401)

    def test_not_allowed(self) -> None:
        with self.http() as http:
            client_id = self.register(http)
            back = self.login(http, client_id, "mallory")
        self.assertEqual(back["error"], "access_denied")
        self.assertNotIn("code", back)

    def test_consent_page(self) -> None:
        from urllib.parse import parse_qs, urlparse
        with self.http() as http:
            client_id = self.register(http)
            # the page names the client and the place the code goes
            def check(text):
                self.assertIn("test", text)
                self.assertIn("claude.ai", text)
                self.assertNotIn("does not send you back to claude.ai", text)
            response = self.consent(http, self.start(http, client_id), "deny", check)
            self.assertEqual(response.status_code, 303)
            back = parse_qs(urlparse(response.headers["location"]).query)
            self.assertEqual(back["error"], ["access_denied"])

        # a client that sends the code somewhere else: the page warns
        with self.http() as http:
            response = http.post("/register", json={"redirect_uris": ["https://evil.example.com/cb"],
                                                    "token_endpoint_auth_method": "none", "client_name": "evil"})
            evil = response.json()["client_id"]
            url = self.start(http, evil, "https://evil.example.com/cb")
            response = self.consent(http, url, "allow", lambda text: self.assertIn(
                "does not send you back to claude.ai", text))
            self.assertEqual(response.status_code, 303)

        # a form post from another site has no cookie: refused, and it does not go to GitHub
        import re
        with self.http() as http:
            url = self.start(http, client_id)
            page = http.get(url).text
        fields = dict(re.findall(r"name=(state|csrf) value='([^']*)'", page))
        with self.http() as other:
            response = other.post("/oauth/consent", data={**fields, "action": "allow"})
            self.assertEqual(response.status_code, 400)
            # and a form with no csrf field, with the cookie, is refused too
        with self.http() as http:
            url = self.start(http, client_id)
            http.get(url)
            state = parse_qs(urlparse(url).query)["state"][0]
            response = http.post("/oauth/consent", data={"state": state, "action": "allow"})
            self.assertEqual(response.status_code, 400)
            # the GitHub callback of a login that the user did not allow is refused
            response = http.get("/oauth/github/callback", params={"code": "alice-code", "state": state})
            self.assertEqual(response.status_code, 400)

    def test_login_in_another_browser(self) -> None:
        """An attacker allows the consent page in their own browser, then sends the GitHub link to a user. GitHub
        sends the user back with a code at once (the user allowed the app before), but the user's browser does
        not have the cookie of the consent page: refused, and no code goes to the attacker."""
        import base64
        import hashlib
        import secrets
        from urllib.parse import parse_qs, urlparse
        evil = "https://evil.example.com/cb"
        with self.http() as attacker:
            response = attacker.post("/register", json={"redirect_uris": [evil], "token_endpoint_auth_method": "none",
                                                        "client_name": "Claude"})
            client_id = response.json()["client_id"]
            verifier = secrets.token_urlsafe(48)
            challenge = base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest()).rstrip(b"=").decode()
            response = attacker.get("/authorize", params={
                "response_type": "code", "client_id": client_id, "redirect_uri": evil, "state": "x",
                "code_challenge": challenge, "code_challenge_method": "S256"})
            consent = urlparse(response.headers["location"])
            response = self.consent(attacker, f"{consent.path}?{consent.query}", "allow")
            state = parse_qs(urlparse(response.headers["location"]).query)["state"][0]
        with self.http() as victim:
            response = victim.get("/oauth/github/callback", params={"code": "alice-code", "state": state})
            self.assertEqual(response.status_code, 400)
            self.assertNotIn("location", response.headers)
        # and the attacker cannot complete the login with the code of the user either: the server refused
        # before it called GitHub, but the login is still there for the attacker's own browser only
        self.assertTrue(state in self.provider.logins)

    def test_renamed_github_login(self) -> None:
        """A login that is taken by another GitHub account (after a rename) does not get in."""
        with self.http() as http:
            client_id = self.register(http)
            first = self.tokens(http, client_id, "alice")
            old = self.github.ids["alice"]
            self.github.ids["alice"] = 999999
            try:
                back = self.login(http, client_id, "alice")
                self.assertEqual(back["error"], "access_denied")
            finally:
                self.github.ids["alice"] = old
        # the tokens of the first account still work
        self.assertFalse(self.call(first["access_token"], "list_runs", {}).is_error)

    def test_bad_callback(self) -> None:
        with self.http() as http:
            response = http.get("/oauth/github/callback", params={"code": "alice-code", "state": "unknown"})
            self.assertEqual(response.status_code, 400)
            response = http.get("/oauth/github/callback", params={"error": "access_denied"})
            self.assertEqual(response.status_code, 400)

    def test_code_once_and_pkce(self) -> None:
        with self.http() as http:
            client_id = self.register(http)
            back = self.login(http, client_id, "alice")
            data = {"grant_type": "authorization_code", "code": back["code"], "redirect_uri": self.REDIRECT,
                    "client_id": client_id, "code_verifier": "wrong" * 10}
            self.assertEqual(http.post("/token", data=data).status_code, 400)
            data["code_verifier"] = back["verifier"]
            self.assertEqual(http.post("/token", data=data).status_code, 200)
            # a code is used one time only
            self.assertEqual(http.post("/token", data=data).status_code, 400)

    def test_refresh_and_remove_user(self) -> None:
        with self.http() as http:
            client_id = self.register(http)
            first = self.tokens(http, client_id, "bob")
            data = {"grant_type": "refresh_token", "refresh_token": first["refresh_token"], "client_id": client_id}
            response = http.post("/token", data=data)
            self.assertEqual(response.status_code, 200, response.text)
            second = response.json()
            self.assertNotEqual(second["access_token"], first["access_token"])
            # the refresh token is rotated
            self.assertEqual(http.post("/token", data=data).status_code, 400)
            self.assertFalse(self.call(second["access_token"], "list_runs", {}).is_error)

            # the state file has no tokens in clear, and a new provider knows the client and the refresh token
            text = self.state_file.read_text()
            self.assertNotIn(second["refresh_token"], text)
            self.assertNotIn(second["access_token"], text)
            from nismcp.oauth import AllowList, GitHubOAuthProvider
            again = GitHubOAuthProvider(FakeGitHub(), AllowList(self.users), self.public, state_file=self.state_file)
            client = asyncio.run(again.get_client(client_id))
            self.assertIsNotNone(client)
            self.assertIsNotNone(asyncio.run(again.load_refresh_token(client, second["refresh_token"])))
            self.assertIsNotNone(asyncio.run(again.load_access_token(second["access_token"])))

            # remove bob: his token stops at once, and he cannot refresh it
            text = self.users.read_text()
            time.sleep(0.01)
            self.users.write_text("alice\n")
            try:
                response = http.post("/mcp", json={"jsonrpc": "2.0", "id": 1, "method": "tools/list"},
                                     headers={"Authorization": f"Bearer {second['access_token']}"})
                self.assertEqual(response.status_code, 401)
                data["refresh_token"] = second["refresh_token"]
                self.assertEqual(http.post("/token", data=data).status_code, 400)
            finally:
                time.sleep(0.01)
                self.users.write_text(text)

    def test_api_key_too(self) -> None:
        for header in ("authorization", "x-api-key"):
            result = self.call(self.key, "irradiate", {"atoms": {"H3": 1.0}, "irradiation_time": 60,
                                                       "flux_n_cm2_s": 0}, header=header)
            self.assertFalse(result.is_error, _text(result))


@unittest.skipUnless(_has_mcp and has_data(JEFF311_DECAY) and has_data(JEFF4_REACTION),
                     "needs the mcp group and the decay and reaction data")
class MCPRealDataTest(unittest.TestCase):
    """1 kg of U235 in a 14 MeV flux of 1e14 n/cm2/s for 1 day, with the real data."""

    def test_u235_14mev(self) -> None:
        from nismcp.data import DataFiles, DataLoader
        from nismcp.materials import atoms_from_mass
        from nismcp.model import Runner, spectrum_from_input, steps_from_input, top_nuclides, totals
        from nismcp.plot import plot_png
        from nisbench.engines import reaction_cache
        files = DataFiles(decay=JEFF311_DECAY, reaction=reaction_cache(JEFF4_REACTION),
                          fission_yields=JEFF311_NFY if has_data(JEFF311_NFY, 1000) else None)
        loader = DataLoader(files).start()
        data, tables = loader.get()
        atoms, _ = atoms_from_mass(1.0, {"U235": 1.0}, tables.index)
        steps = steps_from_input("1 d", ["1 h", "1 d"], None, 1e14)
        result = Runner(loader).run("t", atoms, steps, spectrum_from_input("14MeV"), {}, [])
        end = result.point_at(86400)
        self.assertGreater(totals(tables, end.atoms)["activity_Bq"], 1e12)
        products = top_nuclides(tables, end.atoms, "activity", 3, exclude={"U235"})
        self.assertEqual(len(products), 3)
        self.assertTrue(plot_png(result, tables, products).startswith(b"\x89PNG"))


if __name__ == "__main__":
    unittest.main()
