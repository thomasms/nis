"""
Run the cases of a benchmark with NIS. Two engines give the same results:

    pynis   the python bindings, in this process. The nuclear data is read one time for all the cases.
    nis     the nis program, one run of the deck (<case>.inp) for each case.

Each run gives, at the end of each step, the decay heat (W, total and by type) and the activity (Bq).
"""

import multiprocessing
import os
import struct
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from .format import Case

# NIS result names -> quantities of the format
_HEAT_PARTS = {"total": "decay_heat", "light_particle": "beta_heat",
               "electromagnetic": "gamma_heat", "heavy_particle": "alpha_heat"}


@dataclass
class NuclearData:
    decay: Path
    reaction: Path = None
    fission_yields: Path = None
    sf_yields: Path = None


# the header of a reaction cache, as src/nis/endf/reaction_cache.cpp writes it
_CACHE_MAGIC = b"NISREAC\0"
_CACHE_VERSION = 1


def _cache_version(path: Path):
    """The version of the reaction cache, or None if the file is not a cache."""
    try:
        with open(path, "rb") as f:
            header = f.read(len(_CACHE_MAGIC) + 4)
    except OSError:
        return None
    if len(header) < len(_CACHE_MAGIC) + 4 or not header.startswith(_CACHE_MAGIC):
        return None
    return struct.unpack("=I", header[len(_CACHE_MAGIC):])[0]


def reaction_cache(reaction_data: Path, nis_cache: Path = None) -> Path:
    """The reaction cache to use for the reaction data.

    If the reaction data is a cache, it is used as it is. Else the cache is <reaction data>.nisreac,
    and it is made first if it does not exist, is older than the reaction data, or is of another version.
    It is made with the nis_cache program if one is given, else with the python bindings."""
    reaction_data = Path(reaction_data)
    if _cache_version(reaction_data) is not None:
        return reaction_data
    cache = reaction_data.with_suffix(".nisreac")
    if (_cache_version(cache) == _CACHE_VERSION
            and cache.stat().st_mtime >= reaction_data.stat().st_mtime):
        return cache

    print(f"making the reaction cache {cache}", file=sys.stderr)
    start = time.time()
    # write to a temporary file first, so a stopped run does not leave half a cache
    temporary = cache.with_suffix(".nisreac.tmp")
    if nis_cache is not None and Path(nis_cache).exists():
        subprocess.run([str(nis_cache), "-r", str(reaction_data), "-o", str(temporary)], check=True,
                       stdout=subprocess.DEVNULL)
    else:
        import pynis as nis
        nis.core.write_reaction_cache(nis.core.reaction_data_from_file(str(reaction_data)), str(temporary))
    os.replace(temporary, cache)
    print(f"reaction cache made in {time.time() - start:.1f} s", file=sys.stderr)
    return cache


# ---------------------------------------------------------------- pynis

class _PynisData:
    def __init__(self, data: NuclearData):
        import pynis as nis
        self.nis = nis
        self.fission_yields = nis.core.fission_yields_from_file(str(data.fission_yields)) if data.fission_yields else None
        sf_yields = nis.core.fission_yields_from_file(str(data.sf_yields)) if data.sf_yields else None
        self.decay_mat = nis.core.decay_mat_from_file(str(data.decay), sf_yields)
        self.energy_data = nis.core.decay_energy_data_from_file(str(data.decay))
        self.reaction_data = nis.core.reaction_data_from_file(str(data.reaction)) if data.reaction else None


def _run_pynis(case: Case, d: _PynisData) -> list:
    nis = d.nis
    reaction_mat, spectrum_total = None, None
    if case.spectrum is not None and any(flux > 0 for _, flux in case.steps):
        if d.reaction_data is None:
            raise RuntimeError(f"{case.name}: the case has a flux, so it needs reaction data")
        bounds, values = case.spectrum
        reaction_mat = nis.core.reaction_mat(d.reaction_data, nis.core.Spectrum(bounds, values), d.fission_yields)
        spectrum_total = sum(values)

    matrices = {}
    nuclides = dict(case.atoms)
    results = []
    for dt, flux in case.steps:
        if flux not in matrices:
            # the same as solver::flux_scale: the spectrum is in n/cm2/s and its total is the step flux
            matrices[flux] = d.decay_mat if flux == 0 else d.decay_mat + reaction_mat * (flux * 1e-4 / spectrum_total)
        model = nis.core.BasicModel()
        model.nuclides = nuclides
        model.dt = dt
        nuclides = nis.core.compute_inventory(model, matrices[flux], nis.core.solver_cram48).nuclides
        inventory = nis.core.vector_from_inventory(nuclides)
        heat = nis.core.compute_decay_heat(d.energy_data, inventory)
        result = {quantity: getattr(heat, part) for part, quantity in _HEAT_PARTS.items()}
        result["activity"] = nis.core.compute_activity(d.energy_data, inventory)
        results.append(result)
    return results


# set before the worker processes are forked, so they share the data and do not read it again
_SHARED = {}


def _run_shared(index: int) -> list:
    return _run_pynis(_SHARED["cases"][index], _SHARED["data"])


def run_pynis(cases: list, data: NuclearData, jobs: int = 1):
    """Yields the results of each case, in the same order as the cases."""
    start = time.time()
    _SHARED.update(cases=cases, data=_PynisData(data))
    print(f"nuclear data read in {time.time() - start:.1f} s", file=sys.stderr)
    if jobs <= 1 or len(cases) <= 1:
        yield from (_run_shared(i) for i in range(len(cases)))
        return
    # fork, so the workers get the data that is already read (Linux and macOS only)
    with multiprocessing.get_context("fork").Pool(jobs) as pool:
        yield from pool.imap(_run_shared, range(len(cases)))


# ---------------------------------------------------------------- the nis program

def parse_nis_output(output: str) -> list:
    """The decay heat and activity at the end of each step of the nis output."""
    results = []
    section = None
    for line in output.splitlines():
        words = line.split()
        if line.startswith("step "):
            results.append({})
            section = None
        elif len(words) == 1 and not line.startswith(" "):
            section = words[0]
        elif results and len(words) == 2:
            if section == "decay_heat" and words[0] in _HEAT_PARTS:
                results[-1][_HEAT_PARTS[words[0]]] = float(words[1])
            elif section == "activity" and words[0] == "total":
                results[-1]["activity"] = float(words[1])
    return results


def _run_program(case: Case, nis: Path, data: NuclearData) -> list:
    if not case.deck_path.exists():
        raise RuntimeError(f"{case.name}: no deck {case.deck_path}, make it with scripts/bench.py decks")
    # a large threshold, so the output has no nuclide list
    command = [str(nis), "-d", str(data.decay), "-i", str(case.deck_path), "-t", "1e300"]
    if data.reaction:
        command += ["-r", str(data.reaction)]
    if data.fission_yields:
        command += ["-y", str(data.fission_yields)]
    if data.sf_yields:
        command += ["--sf-yields", str(data.sf_yields)]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{case.name}: {result.stderr.strip()}")
    results = parse_nis_output(result.stdout)
    if len(results) != len(case.steps):
        raise RuntimeError(f"{case.name}: the nis output has {len(results)} steps, the case has {len(case.steps)}")
    return results


def run_program(cases: list, nis: Path, data: NuclearData, jobs: int = 1):
    """Yields the results of each case, in the same order as the cases."""
    with ThreadPoolExecutor(max(jobs, 1)) as pool:
        yield from pool.map(lambda c: _run_program(c, nis, data), cases)


def run(cases: list, data: NuclearData, engine: str = "pynis", nis: Path = None, jobs: int = 1):
    if engine == "pynis":
        return run_pynis(cases, data, jobs)
    if engine == "nis":
        if nis is None or not Path(nis).exists():
            raise SystemExit(f"the nis program {nis} does not exist, build it or give --nis")
        return run_program(cases, nis, data, jobs)
    raise ValueError(f"unknown engine {engine!r}")
