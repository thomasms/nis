"""
Make NIS reaction data from raw ENDF neutron files, with NJOY and PREPRO (see processing/README.md).

Commands:
    tools      get and build NJOY2016 and PREPRO
    fetch      download the raw ENDF files of a library (TENDL)
    run        process ENDF files to 1102 group cross sections, one file for each nuclide (in parallel)
    combine    join the processed files into one reaction data file (and make the NIS cache)
    tendl      all of it for a TENDL release: fetch, run, combine

Examples:
    uv run python processing/process.py tools
    uv run python processing/process.py tendl --year 2025 --output data/tendl2025.reac
    uv run python processing/process.py run path/to/endf --output processing/out/mylib -j 16
    uv run python processing/process.py combine processing/out/mylib --output data/mylib.reac --title "My library"
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

from nisdata import endf, fetch, tools as build_tools
from nisdata.pipeline import Tools
from nisdata.runner import is_archive, process_all, process_archive

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_TOOLS = HERE / "tools"
DEFAULT_RAW = HERE / "raw"
DEFAULT_OUT = HERE / "out"
DEFAULT_NIS_CACHE = ROOT / "build" / "bin" / "nis_cache"


def _tools(args) -> Tools:
    """--njoy and --prepro-dir, then NJOY_EXE and PREPRO_EXE_DIR, then the programs from 'tools'."""
    default = build_tools.default_tools(args.tools)
    njoy = args.njoy or os.environ.get("NJOY_EXE") or default.njoy
    prepro = args.prepro_dir or os.environ.get("PREPRO_EXE_DIR") or default.prepro_dir
    return Tools(njoy=Path(njoy), prepro_dir=Path(prepro))


def _select(files: list, only: list) -> list:
    if not only:
        return files
    wanted = set(only)
    return [f for f in files if endf.read_material(f).name in wanted]


def _report(results: list, output: Path) -> int:
    failed = [r for r in results if not r.ok]
    done = len(results) - len(failed)
    print(f"{done} of {len(results)} nuclides done, results in {output} (summary.json)", file=sys.stderr)
    if failed:
        print(f"{len(failed)} failed (the logs are in {output / 'logs'}):", file=sys.stderr)
        for r in sorted(failed, key=lambda r: r.name or str(r.source)):
            print(f"    {r.name or r.source}: {r.error}", file=sys.stderr)
    return 1 if failed else 0


def _combine(folder: Path, output: Path, title: str, nis_cache: Path) -> None:
    files = sorted(Path(folder).glob("*.asc"))
    if not files:
        raise SystemExit(f"no processed files (*.asc) in {folder}")
    count = endf.combine(files, output, title)
    print(f"{count} materials in {output}", file=sys.stderr)
    if nis_cache and Path(nis_cache).is_file():
        cache = Path(output).with_suffix(".nisreac")
        subprocess.run([str(nis_cache), "--reaction-data", str(output), "--output", str(cache)], check=True)
        print(f"NIS cache in {cache}", file=sys.stderr)
    elif nis_cache:
        print(f"no {nis_cache}, so no NIS cache made (build NIS, or give --nis-cache)", file=sys.stderr)


def cmd_tools(args) -> int:
    t = build_tools.build(args.tools, args.jobs, args.njoy_ref, args.prepro_ref, args.fortran, args.flags)
    print(f"NJOY: {t.njoy}\nPREPRO: {t.prepro_dir}", file=sys.stderr)
    return 0


def cmd_fetch(args) -> int:
    path = fetch.fetch_tendl(args.year, args.raw, args.url, unpack=not args.no_unpack)
    if path.is_dir():
        print(f"{len(endf.find_endf_files(path))} neutron files in {path}", file=sys.stderr)
    else:
        print(f"the archive is {path} (run can read it, it is not necessary to unpack it)", file=sys.stderr)
    return 0


def _process(sources: list, output: Path, args) -> list:
    """Folders and files of ENDF files, or one tar archive of them."""
    if len(sources) == 1 and is_archive(sources[0]):
        return process_archive(sources[0], output, _tools(args), args.jobs, args.temperature, args.work_dir,
                               args.keep, args.force, only=tuple(args.only))
    files = []
    for source in sources:
        files += endf.find_endf_files(source) if source.is_dir() else [source]
    files = _select(files, args.only)
    if not files:
        raise SystemExit("no ENDF neutron files")
    return process_all(files, output, _tools(args), args.jobs, args.temperature, args.work_dir, args.keep,
                       args.force)


def cmd_run(args) -> int:
    return _report(_process(args.sources, args.output, args), args.output)


def cmd_combine(args) -> int:
    _combine(args.folder, args.output, args.title, args.nis_cache)
    return 0


def cmd_tendl(args) -> int:
    # the files are read from the archive, so the 14 GB of TENDL are never all on the disk
    archive = fetch.fetch_tendl(args.year, args.raw, args.url, unpack=False)
    out = args.processed or DEFAULT_OUT / f"tendl{args.year}"
    status = _report(_process([archive], out, args), out)
    output = args.output or ROOT / "data" / f"tendl{args.year}.reac"
    _combine(out, output, f"TENDL-{args.year} Incident Neutron File (1102 groups)", args.nis_cache)
    return status


def _add_run_options(p) -> None:
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 1,
                   help="nuclides at the same time (default: the number of cores)")
    p.add_argument("--only", nargs="*", default=[], help="only these nuclides, e.g. Fe56 Co59 Ag110m")
    p.add_argument("--temperature", type=float, default=293.6, help="K")
    p.add_argument("--tools", type=Path, default=DEFAULT_TOOLS, help="the folder of 'tools'")
    p.add_argument("--njoy", type=Path, default=None, help="the NJOY program (or NJOY_EXE)")
    p.add_argument("--prepro-dir", type=Path, default=None, help="the folder of the PREPRO programs (or PREPRO_EXE_DIR)")
    p.add_argument("--work-dir", type=Path, default=None,
                   help="for the temporary files, default /dev/shm (memory) if there is one")
    p.add_argument("--keep", action="store_true", help="keep the work folders")
    p.add_argument("--force", action="store_true", help="process again the nuclides that are done")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    p = commands.add_parser("tools", help="get and build NJOY2016 and PREPRO")
    p.add_argument("--tools", type=Path, default=DEFAULT_TOOLS, help="where to put them")
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 1)
    p.add_argument("--njoy-ref", default=None, help="a tag or branch of NJOY2016, default the latest")
    p.add_argument("--prepro-ref", default=None, help="a tag or branch of PREPRO, default the latest")
    p.add_argument("--fortran", default="gfortran")
    p.add_argument("--flags", default="-O2", help="flags for the PREPRO codes")
    p.set_defaults(func=cmd_tools)

    p = commands.add_parser("fetch", help="download the TENDL neutron files")
    p.add_argument("--year", default=fetch.LATEST_TENDL, help=f"the TENDL release, default {fetch.LATEST_TENDL}")
    p.add_argument("--raw", type=Path, default=DEFAULT_RAW, help="where to put them")
    p.add_argument("--url", default=None, help="the address of the archive, if not the usual one")
    p.add_argument("--no-unpack", action="store_true", help="only download the archive (about 14 GB less disk)")
    p.set_defaults(func=cmd_fetch)

    p = commands.add_parser("run", help="process ENDF files (or folders of them)")
    p.add_argument("sources", nargs="+", type=Path,
                   help="ENDF files, folders of them, or one tar archive of them (e.g. TENDL-n.tgz, not unpacked)")
    p.add_argument("--output", type=Path, required=True, help="the folder for the processed files")
    _add_run_options(p)
    p.set_defaults(func=cmd_run)

    p = commands.add_parser("combine", help="join the processed files into one file")
    p.add_argument("folder", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--title", default="Incident Neutron File (1102 groups)", help="the first line of the file")
    p.add_argument("--nis-cache", type=Path, default=DEFAULT_NIS_CACHE, help="the nis_cache program")
    p.set_defaults(func=cmd_combine)

    p = commands.add_parser("tendl", help="fetch, run and combine a TENDL release")
    p.add_argument("--year", default=fetch.LATEST_TENDL)
    p.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    p.add_argument("--url", default=None)
    p.add_argument("--processed", type=Path, default=None, help="default processing/out/tendl<year>")
    p.add_argument("--output", type=Path, default=None, help="default data/tendl<year>.reac")
    p.add_argument("--nis-cache", type=Path, default=DEFAULT_NIS_CACHE)
    _add_run_options(p)
    p.set_defaults(func=cmd_tendl)

    args = parser.parse_args()
    try:
        return args.func(args)
    except (ValueError, FileNotFoundError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
