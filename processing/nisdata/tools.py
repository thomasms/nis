"""
Get and build rustjoy and PREPRO (needs git, cargo and gfortran).

rustjoy (a Rust port of NJOY2016) does the NJOY step: it gives the same tapes as NJOY2016, byte for byte, and is
faster. The NJOY2016 repository is cloned too (it is not built), because the processing tests use its test files.
"""

import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .pipeline import PREPRO_CODES, Tools

RUSTJOY_REPO = "https://github.com/thomasms/rustjoy"
NJOY_REPO = "https://github.com/njoy/NJOY2016"
PREPRO_REPO = "https://github.com/IAEA-NDS/PREPRO"
# the PREPRO files that each code needs with its own source
PREPRO_COMMON = ("endfio.f", "scratchb.f")


def default_tools(prefix: Path) -> Tools:
    return Tools(njoy=Path(prefix) / "bin" / "rustjoy", prepro_dir=Path(prefix) / "prepro")


def _run(command: list, cwd: Path = None) -> None:
    print("+", " ".join(str(c) for c in command), file=sys.stderr, flush=True)
    subprocess.run([str(c) for c in command], cwd=cwd, check=True)


def _clone(repo: str, folder: Path, ref: str = None) -> None:
    if folder.is_dir():
        print(f"{folder} is there, not cloned again", file=sys.stderr)
        return
    command = ["git", "clone", "--depth", "1"] + (["--branch", ref] if ref else []) + [repo, folder]
    _run(command)


def build_rustjoy(source: Path, prefix: Path, jobs: int) -> Path:
    _run(["cargo", "build", "--release", "--locked", "-j", jobs], cwd=source)
    target = Path(prefix) / "bin" / "rustjoy"
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        target.unlink()
    target.symlink_to((source / "target" / "release" / "rustjoy").resolve())
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


def build(prefix: Path, jobs: int = None, rustjoy_ref: str = None, prepro_ref: str = None,
          fortran: str = "gfortran", flags: str = "-O2") -> Tools:
    """Clone into <prefix>/src and build. Run it again to build again (the clones are kept)."""
    for program in ("git", "cargo", fortran):
        if shutil.which(program) is None:
            raise SystemExit(f"{program} is needed to build rustjoy and PREPRO")
    jobs = jobs or os.cpu_count() or 1
    prefix = Path(prefix)
    sources = prefix / "src"
    sources.mkdir(parents=True, exist_ok=True)
    _clone(RUSTJOY_REPO, sources / "rustjoy", rustjoy_ref)
    _clone(NJOY_REPO, sources / "NJOY2016")  # only for the test files
    _clone(PREPRO_REPO, sources / "PREPRO", prepro_ref)
    build_rustjoy(sources / "rustjoy", prefix, jobs)
    build_prepro(sources / "PREPRO", prefix, jobs, fortran, flags)
    tools = default_tools(prefix)
    tools.check()
    return tools
