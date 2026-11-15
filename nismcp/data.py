"""
The nuclear data of the server: read one time at the start, then used by all the runs.

The reaction data can be an ENDF file or a reaction cache (.nisreac, made with nis_cache). Use the cache: it is
read in a few seconds, and it needs less memory than the ENDF file.
"""

import logging
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

log = logging.getLogger("nismcp.data")


@dataclass
class DataFiles:
    decay: Path
    reaction: Path = None
    fission_yields: Path = None
    sf_yields: Path = None

    def describe(self) -> dict:
        def name(path):
            return None if path is None else Path(path).name
        return {"decay": name(self.decay), "reaction": name(self.reaction),
                "fission_yields": name(self.fission_yields), "sf_yields": name(self.sf_yields)}


def _check_file(path: Path, what: str) -> Path:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"the {what} file {path} does not exist")
    # without git lfs, the data file is a small text pointer
    with open(path, "rb") as f:
        if f.read(40).startswith(b"version https://git-lfs"):
            raise ValueError(f"the {what} file {path} is a git lfs pointer, get the data with 'git lfs pull'")
    return path


class NuclearData:
    """The decay matrix, the decay energies, the reaction data and the fission yields.

    For the tests, a subclass can give the matrices without files (see reaction_matrix)."""

    def __init__(self, files: DataFiles):
        import pynis as nis
        self.nis = nis
        self.files = files
        start = time.time()
        sf_yields = None
        if files.sf_yields:
            sf_yields = nis.core.fission_yields_from_file(str(_check_file(files.sf_yields, "spontaneous fission yields")))
        decay = str(_check_file(files.decay, "decay data"))
        self.decay_mat = nis.core.decay_mat_from_file(decay, sf_yields)
        self.energy_data = nis.core.decay_energy_data_from_file(decay)
        self.fission_yields = None
        if files.fission_yields:
            self.fission_yields = nis.core.fission_yields_from_file(
                str(_check_file(files.fission_yields, "fission yields")))
        self.reaction_data = None
        if files.reaction:
            reaction = str(_check_file(files.reaction, "reaction data"))
            if not nis.core.is_reaction_cache(reaction):
                log.warning("%s is an ENDF file: make a cache with nis_cache, it is faster and needs less memory",
                            reaction)
            self.reaction_data = nis.core.reaction_data_from_file(reaction)
        log.info("nuclear data read in %.1f s", time.time() - start)

    @property
    def has_reactions(self) -> bool:
        return self.reaction_data is not None

    @property
    def has_fission_yields(self) -> bool:
        return self.fission_yields is not None

    def reaction_matrix(self, bounds: list, values: list):
        """The reaction matrix for the spectrum (group bounds in eV, ascending, and the flux in n/cm2/s per group)."""
        if self.reaction_data is None:
            raise ValueError("the server has no reaction data, so it can only do decay (flux 0)")
        spectrum = self.nis.core.Spectrum(list(bounds), list(values))
        return self.nis.core.reaction_mat(self.reaction_data, spectrum, self.fission_yields)

    def describe(self) -> dict:
        return self.files.describe()


@dataclass
class _Tables:
    """Tables from the decay data that the results need: names, decay constants and decay energies."""
    names: list
    index: dict
    decay_constants: np.ndarray
    # the mass number of each nuclide, for the mass of a nuclide (g, approximate: A g/mol)
    mass_numbers: np.ndarray
    # decay energy per decay in J: total, light particle (beta), electromagnetic (gamma), heavy particle (alpha)
    energy: dict = field(default_factory=dict)


def tables(data: NuclearData) -> _Tables:
    nis = data.nis
    zais = nis.data.ALL_ZAIS()
    names = [nis.zai.get_name_from(z) for z in zais]
    ev_to_j = 1.602176634e-19
    e = data.energy_data
    energy = {"light_particle": np.asarray(e.light_particle) * ev_to_j,
              "electromagnetic": np.asarray(e.electromagnetic) * ev_to_j,
              "heavy_particle": np.asarray(e.heavy_particle) * ev_to_j}
    energy["total"] = energy["light_particle"] + energy["electromagnetic"] + energy["heavy_particle"]
    mass_numbers = np.array([nis.zai.get_a(z) for z in zais], dtype=float)
    return _Tables(names=names, index={n: i for i, n in enumerate(names)},
                   decay_constants=np.asarray(e.decay_constants), mass_numbers=mass_numbers, energy=energy)


class DataLoader:
    """Reads the nuclear data in a background thread, so the server can answer (for example the list of tools)
    while the data is read. A tool that needs the data waits for it."""

    def __init__(self, files: DataFiles, factory=NuclearData):
        self.files = files
        self._factory = factory
        self._ready = threading.Event()
        self._data = None
        self._tables = None
        self._error = None
        self.started = None

    def start(self) -> "DataLoader":
        """Start to read the data. A second call does nothing."""
        if self.started is not None:
            return self
        self.started = time.time()
        threading.Thread(target=self._load, name="nis-data", daemon=True).start()
        return self

    def _load(self) -> None:
        try:
            data = self._factory(self.files)
            self._tables = tables(data)
            self._data = data
        except Exception as e:  # noqa: BLE001, the error goes to the tools
            log.exception("cannot read the nuclear data")
            self._error = e
        finally:
            self._ready.set()

    @property
    def ready(self) -> bool:
        return self._ready.is_set() and self._error is None

    def status(self) -> str:
        if self.started is None:
            return "not started"
        if not self._ready.is_set():
            return f"loading ({time.time() - self.started:.0f} s)"
        return f"error: {self._error}" if self._error else "ready"

    def get(self, timeout: float = 600.0):
        """(data, tables), after the data is read."""
        if not self._ready.wait(timeout):
            raise TimeoutError("the nuclear data is still loading, try again in a minute")
        if self._error is not None:
            raise RuntimeError(f"the server cannot read its nuclear data: {self._error}")
        return self._data, self._tables
