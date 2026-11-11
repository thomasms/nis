"""
Process one ENDF neutron file to 1102 group cross sections (INT=1), the form that NIS reads.

    NJOY:    RECONR, BROADR, UNRESR, HEATR (kerma, damage), GASPR  ->  NJOY.OUT (PENDF)
    PREPRO:  LINEAR, RECENT, SIGMA1, SIXPAK (+ LINEAR of MF9), MERGER, ACTIVATE  ->  ACTIVATE.OUT
    then:    MERGER (the NJOY heating, damage and gas data with the PREPRO data), MERGER, MERGER,
             FIXUP, DICTIN, GROUPIE (1102 groups)  ->  the output

All the steps of a nuclide run in one work folder. NJOY runs at the same time as the first PREPRO
steps (they do not need each other), so a nuclide uses up to 2 cores for part of its time.
"""

import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

from . import decks
from .endf import Material, mat_mf_mt, read_material

PREPRO_CODES = ("linear", "recent", "sigma1", "sixpak", "merger", "activate", "fixup", "dictin", "groupie")


class StepError(RuntimeError):
    pass


@dataclass
class Tools:
    njoy: Path
    prepro_dir: Path

    def prepro(self, code: str) -> Path:
        return self.prepro_dir / code

    def check(self) -> None:
        missing = [str(p) for p in [self.njoy] + [self.prepro(c) for c in PREPRO_CODES]
                   if not (p.is_file() and os.access(p, os.X_OK))]
        if missing:
            raise FileNotFoundError("these programs are missing (build them with: process.py tools): "
                                    + ", ".join(missing))


@dataclass
class Result:
    source: Path
    name: str = ""
    output: Path = None
    ok: bool = False
    skipped: bool = False
    seconds: float = 0.0
    error: str = ""
    # time of each step (s)
    steps: dict = field(default_factory=dict)


def default_work_root() -> Path:
    """A memory file system if there is one (PREPRO reads and writes many scratch files)."""
    shm = Path("/dev/shm")
    if shm.is_dir() and os.access(shm, os.W_OK):
        return shm
    return Path(tempfile.gettempdir())


def _log(work: Path, text: str) -> None:
    with open(work / "steps.log", "a") as f:
        f.write(text)


def _prepro(tools: Tools, code: str, deck: str, work: Path, outputs: list, timings: dict, label: str = "") -> None:
    """Run one PREPRO code with its input deck (<CODE>.INP) in the work folder."""
    start = time.perf_counter()
    (work / f"{code.upper()}.INP").write_text(deck)
    run = subprocess.run([str(tools.prepro(code))], cwd=work, stdin=subprocess.DEVNULL,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log = run.stdout.decode(errors="replace")
    _log(work, f"==== {label or code} (exit {run.returncode})\n{log}\n")
    timings[label or code] = time.perf_counter() - start
    missing = [o for o in outputs if not (work / o).is_file() or (work / o).stat().st_size == 0]
    # a PREPRO code can stop on an error in the data with exit code 0
    if run.returncode != 0 or missing or "Execution Terminated" in log:
        raise StepError(f"{label or code} failed (exit {run.returncode}"
                        f"{', no ' + ', '.join(missing) if missing else ''}): {error_text(log)}")


def error_text(log: str) -> str:
    """The first error message in the output of a program, else its last line."""
    lines = [line.strip() for line in log.splitlines() if line.strip()]
    for i, line in enumerate(lines):
        if line.startswith(("ERROR", "Fortran runtime error", "***error", "*** error")) or "error in" in line:
            return " ".join(lines[i:i + 2])
    return lines[-1] if lines else "no output"


def _start_njoy(tools: Tools, material: Material, source: Path, work: Path, temperature: float):
    folder = work / "njoy"
    folder.mkdir()
    (folder / "tape20").symlink_to(source.resolve())
    deck = decks.njoy(material.mat, temperature)
    (folder / "input").write_text(deck)
    log = open(folder / "stdout", "wb")
    with open(folder / "input") as stdin:
        process = subprocess.Popen([str(tools.njoy)], cwd=folder, stdin=stdin, stdout=log, stderr=subprocess.STDOUT)
    return process, log, time.perf_counter()


def _finish_njoy(process, log, start, work: Path, timings: dict) -> Path:
    code = process.wait()
    log.close()
    timings["njoy"] = time.perf_counter() - start
    folder = work / "njoy"
    output = folder / "tape35"
    text = (folder / "stdout").read_text(errors="replace")
    _log(work, f"==== njoy (exit {code})\n{text}\n")
    if code != 0 or not output.is_file() or output.stat().st_size == 0:
        raise StepError(f"njoy failed (exit {code}): {error_text(text)}")
    target = work / "NJOY.OUT"
    output.rename(target)
    return target


def check_output(path: Path, material: Material) -> None:
    """The output has the material and the 1102 group total cross section (MF3/MT1)."""
    with open(path) as f:
        f.readline()
        first = f.readline()
        if mat_mf_mt(first)[0] != material.mat:
            raise StepError(f"the output has MAT {mat_mf_mt(first)[0]}, not {material.mat}")
        for line in f:
            if mat_mf_mt(line)[1:] == (3, 1):
                points = int(f.readline()[55:66])
                if points != decks.NGROUPS + 1:
                    raise StepError(f"MF3/MT1 has {points} points, not {decks.NGROUPS + 1}")
                return
    raise StepError("the output has no MF3/MT1")


def process_file(source: Path, output_dir: Path, tools: Tools, temperature: float = decks.DEFAULT_TEMPERATURE,
                 work_root: Path = None, keep: bool = False, force: bool = False) -> Result:
    """Process one ENDF file into output_dir/<name>.asc. The logs of a failure are in output_dir/logs."""
    source = Path(source)
    result = Result(source=source)
    start = time.perf_counter()
    try:
        material = read_material(source)
    except (ValueError, OSError) as e:
        result.error = str(e)
        return result
    result.name = material.name
    result.output = Path(output_dir) / f"{material.name}.asc"
    if result.output.is_file() and not force:
        result.ok = result.skipped = True
        return result

    work = Path(tempfile.mkdtemp(prefix=f"nisdata_{material.name}_", dir=work_root or default_work_root()))
    timings = result.steps
    njoy = None
    try:
        njoy = _start_njoy(tools, material, source, work, temperature)
        (work / "ENDF.IN").symlink_to(source.resolve())

        _prepro(tools, "linear", decks.linear("ENDF.IN", "LINEAR.OUT", "1.00000-04"), work, ["LINEAR.OUT"], timings)
        _prepro(tools, "recent", decks.recent("LINEAR.OUT", "RECENT.OUT"), work, ["RECENT.OUT"], timings)
        _prepro(tools, "sigma1", decks.sigma1("RECENT.OUT", "SIGMA1.OUT", temperature), work, ["SIGMA1.OUT"],
                timings)
        _prepro(tools, "sixpak", decks.sixpak("SIGMA1.OUT"), work, [], timings)
        sources = ["SIGMA1.OUT"]
        # the isomer production (MF9) from the MF6 data, if the file has it
        if (work / "ENDFB.MF9").is_file() and (work / "ENDFB.MF9").stat().st_size > 0:
            _prepro(tools, "linear", decks.linear("ENDFB.MF9", "MF9.LINEAR", "1.00000-02"), work, ["MF9.LINEAR"],
                    timings, "linear (MF9)")
            sources.append("MF9.LINEAR")
        _prepro(tools, "merger", decks.merger("MERGER.OUT0", sources, decks.MERGE_PREPRO), work, ["MERGER.OUT0"],
                timings, "merger 0")
        _prepro(tools, "activate", decks.activate("MERGER.OUT0", "ACTIVATE.OUT"), work, ["ACTIVATE.OUT"], timings)

        process, log, njoy_start = njoy
        njoy = None
        _finish_njoy(process, log, njoy_start, work, timings)
        _prepro(tools, "merger", decks.merger("MERGER.OUT1", ["NJOY.OUT"], decks.MERGE_NJOY), work,
                ["MERGER.OUT1"], timings, "merger 1")
        _prepro(tools, "merger", decks.merger("MERGER.OUT2", ["ACTIVATE.OUT", "MERGER.OUT1"], decks.MERGE_BOTH),
                work, ["MERGER.OUT2"], timings, "merger 2")
        _prepro(tools, "merger", decks.merger("MERGER.OUT3", ["MERGER.OUT2"], decks.MERGE_SELECT), work,
                ["MERGER.OUT3"], timings, "merger 3")
        _prepro(tools, "merger", decks.merger("MERGER.OUT4", ["MERGER.OUT3"], decks.MERGE_FINAL), work,
                ["MERGER.OUT4"], timings, "merger 4")
        _prepro(tools, "fixup", decks.fixup("MERGER.OUT4", "FIXUP.OUT"), work, ["FIXUP.OUT"], timings)
        _prepro(tools, "dictin", decks.dictin("FIXUP.OUT", "DICTIN.OUT"), work, ["DICTIN.OUT"], timings)
        _prepro(tools, "groupie", decks.groupie("DICTIN.OUT", "GROUPIE.OUT"), work, ["GROUPIE.OUT"], timings)
        check_output(work / "GROUPIE.OUT", material)

        result.output.parent.mkdir(parents=True, exist_ok=True)
        # write to a temporary name first, so a stopped run never leaves a part file
        partial = result.output.with_suffix(".part")
        shutil.copyfile(work / "GROUPIE.OUT", partial)
        partial.replace(result.output)
        result.ok = True
    except (StepError, OSError) as e:
        result.error = str(e)
    finally:
        if njoy is not None:
            process, log, _ = njoy
            process.kill()
            process.wait()
            log.close()
        if not result.ok:
            logs = Path(output_dir) / "logs"
            logs.mkdir(parents=True, exist_ok=True)
            if (work / "steps.log").is_file():
                shutil.copyfile(work / "steps.log", logs / f"{material.name}.log")
        if not keep:
            shutil.rmtree(work, ignore_errors=True)
        result.seconds = time.perf_counter() - start
    return result
