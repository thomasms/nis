"""Get and build NJOY2016 and PREPRO (needs git, cmake, make and gfortran)."""

import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .pipeline import PREPRO_CODES, Tools

NJOY_REPO = "https://github.com/njoy/NJOY2016"
PREPRO_REPO = "https://github.com/IAEA-NDS/PREPRO"
# the PREPRO files that each code needs with its own source
PREPRO_COMMON = ("endfio.f", "scratchb.f")


def default_tools(prefix: Path) -> Tools:
    return Tools(njoy=Path(prefix) / "bin" / "njoy", prepro_dir=Path(prefix) / "prepro")


def _run(command: list, cwd: Path = None) -> None:
    print("+", " ".join(str(c) for c in command), file=sys.stderr, flush=True)
    subprocess.run([str(c) for c in command], cwd=cwd, check=True)


def _clone(repo: str, folder: Path, ref: str = None) -> None:
    if folder.is_dir():
        print(f"{folder} is there, not cloned again", file=sys.stderr)
        return
    command = ["git", "clone", "--depth", "1"] + (["--branch", ref] if ref else []) + [repo, folder]
    _run(command)


def build_njoy(source: Path, prefix: Path, jobs: int) -> Path:
    build = source / "build"
    _run(["cmake", "-S", source, "-B", build, "-DCMAKE_BUILD_TYPE=Release"])
    _run(["cmake", "--build", build, "-j", jobs])
    target = Path(prefix) / "bin" / "njoy"
    target.parent.mkdir(parents=True, exist_ok=True)
    # NJOY links to its own shared library (libnjoy), keep the build and link the program
    if target.exists() or target.is_symlink():
        target.unlink()
    target.symlink_to((build / "njoy").resolve())
    return target


def build_prepro(source: Path, prefix: Path, jobs: int, fortran: str = "gfortran", flags: str = "-O2") -> Path:
    """
    Compile the PREPRO codes that the processing uses, at the same time. The PREPRO makefile uses -O0;
    -O2 makes the codes about 2 times faster and gives the same output (checked for the test nuclides).
    """
    src = source / "source"
    target = Path(prefix) / "prepro"
    target.mkdir(parents=True, exist_ok=True)

    def compile_one(code: str) -> None:
        _run([fortran] + flags.split() + ["-o", (target / code).resolve(), f"{code}.f", *PREPRO_COMMON], cwd=src)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        list(pool.map(compile_one, PREPRO_CODES))
    return target


def build(prefix: Path, jobs: int = None, njoy_ref: str = None, prepro_ref: str = None,
          fortran: str = "gfortran", flags: str = "-O2") -> Tools:
    """Clone into <prefix>/src and build. Run it again to build again (the clones are kept)."""
    for program in ("git", "cmake", fortran):
        if shutil.which(program) is None:
            raise SystemExit(f"{program} is needed to build NJOY and PREPRO")
    jobs = jobs or os.cpu_count() or 1
    prefix = Path(prefix)
    sources = prefix / "src"
    sources.mkdir(parents=True, exist_ok=True)
    _clone(NJOY_REPO, sources / "NJOY2016", njoy_ref)
    _clone(PREPRO_REPO, sources / "PREPRO", prepro_ref)
    build_njoy(sources / "NJOY2016", prefix, jobs)
    build_prepro(sources / "PREPRO", prefix, jobs, fortran, flags)
    tools = default_tools(prefix)
    tools.check()
    return tools
