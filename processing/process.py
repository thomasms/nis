"""
Get the nuclear data for NIS: download the decay data and the fission yields, and make the reaction data from
raw ENDF neutron files with NJOY and PREPRO (see processing/README.md).

Commands:
    bootstrap  get all the data the tests and the examples need (data/jeff311.dec, endfbvii1.dec, jeff311.nfy,
               jeff311.sfy and jeff4.reac). The reaction data takes a long time (--no-reactions: not it)
    libraries  the libraries that can be downloaded, and their files
    decay      download decay data: data/<library>.dec
    yields     download fission yields: data/<library>.nfy and data/<library>.sfy
    reactions  download the neutron files of a library and make data/<library>.reac (fetch, run, combine)
    tools      get and build NJOY2016 and PREPRO
    fetch      download the raw ENDF files of a library (TENDL)
    run        process ENDF files to 1102 group cross sections, one file for each nuclide (in parallel)
    combine    join the processed files into one reaction data file (and make the NIS cache)
    tendl      all of it for a TENDL release: fetch, run, combine

Examples:
    python3 processing/process.py bootstrap
    python3 processing/process.py decay jeff4 endfbviii0
    python3 processing/process.py yields jeff4
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

from nisdata import endf, fetch, libraries, tools as build_tools
from nisdata.pipeline import BAD_VERSIONS, VERSION, Tools
from nisdata.runner import is_archive, process_all, process_archive

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_TOOLS = HERE / "tools"
DEFAULT_RAW = HERE / "raw"
DEFAULT_OUT = HERE / "out"
DEFAULT_NIS_CACHE = ROOT / "build" / "bin" / "nis_cache"
DEFAULT_DATA = ROOT / "data"


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


def _title(text: str) -> str:
    """The first line of a reaction data file: it gives the version of the processing."""
    return f"{text} (1102 groups, nisdata {VERSION})"


def _up_to_date(reaction_file: Path) -> bool:
    """True if the reaction data file is there and was made by this version of the processing."""
    if not libraries.is_there(reaction_file):
        return False
    with open(reaction_file, errors="replace") as f:
        return f"nisdata {VERSION})" in f.readline()


def _finish(out: Path, output: Path, text: str, status: int, nis_cache: Path) -> int:
    """
    Join the processed files. If some nuclides failed, the file has no version mark, so the next run tries the
    failed nuclides again (the others are done, so they are not processed again).
    """
    if status == 0:
        _combine(out, output, _title(text), nis_cache)
        return 0
    _combine(out, output, f"{text} (1102 groups, INCOMPLETE)", nis_cache)
    print(f"WARNING: {output} does not have the nuclides that failed (see above). Run the same command again to "
          f"try them again", file=sys.stderr)
    return status


def _failed(out: Path) -> int:
    """The number of nuclides that failed in the last run of the processed folder (from its summary.json)."""
    summary = Path(out) / "summary.json"
    if not summary.is_file():
        return 0
    try:
        import json
        return int(json.loads(summary.read_text()).get("failed", 0))
    except (ValueError, OSError):
        return 0


def _use_existing(output: Path, args, out: Path = None) -> bool:
    """
    True if the reaction data file can be used as it is: made by this version, or an older one with --allow-old.
    out: the folder of the processed files; if nuclides failed in its last run, the file is not complete.
    """
    if args.force:
        return False
    if out is not None and libraries.is_there(output) and _failed(out):
        print(f"{output} does not have all the nuclides ({_failed(out)} failed in the last run, see "
              f"{Path(out) / 'summary.json'}), so they are tried again", file=sys.stderr)
        return False
    if _up_to_date(output):
        print(f"{output} is there, not made again (--force makes it again)", file=sys.stderr)
        return True
    if not libraries.is_there(output):
        return False
    with open(output, errors="replace") as f:
        first = f.readline()
    if "INCOMPLETE" in first:
        print(f"{output} does not have all the nuclides (some failed before), so they are tried again",
              file=sys.stderr)
        return False
    bad = [v for v in BAD_VERSIONS if f"nisdata {v})" in first]
    if bad:
        print(f"{output} is from version {bad[0]} of the processing, which has an error in the products of MT=5: "
              f"it is made again (also with --allow-old)", file=sys.stderr)
        return False
    if args.allow_old:
        print(f"{output} is from an older version of the processing (not {VERSION}), it is used as it is "
              f"(--allow-old). It has no products of MT=5 (n,anything): for TENDL they start at 30 MeV, so this "
              f"changes nothing for a spectrum below 30 MeV (FNS)", file=sys.stderr)
        return True
    print(f"{output} is from an older version of the processing (not {VERSION}), so it is made again. It has no "
          f"products of MT=5 (n,anything). To use it as it is, stop now (Ctrl-C) and give --allow-old "
          f"(make: ALLOW_OLD=1)", file=sys.stderr)
    return False


def _check_version(output: Path) -> None:
    """
    The processed files of an older version of the processing are moved to <output>/version-<n>/, so that
    they are made again. A stopped run then continues with the files of this version only.
    """
    output = Path(output)
    stamp = output / "VERSION"
    old = list(output.glob("*.asc"))
    done = stamp.read_text().strip() if stamp.is_file() else ("1" if old else None)
    if done is not None and done != str(VERSION) and old:
        keep = output / f"version-{done}"
        keep.mkdir(exist_ok=True)
        for path in old:
            path.replace(keep / path.name)
        print(f"the {len(old)} files in {output} are from version {done} of the processing, this is version "
              f"{VERSION}: all are made again (the old files are in {keep}, you can delete them)", file=sys.stderr)
    output.mkdir(parents=True, exist_ok=True)
    stamp.write_text(f"{VERSION}\n")


def _process(sources: list, output: Path, args) -> list:
    """Folders and files of ENDF files, or one tar archive of them."""
    _check_version(output)
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
    stamp = Path(args.folder) / "VERSION"
    current = stamp.is_file() and stamp.read_text().strip() == str(VERSION)
    title = args.title or (_title("Incident Neutron File") if current else "Incident Neutron File (1102 groups)")
    _combine(args.folder, args.output, title, args.nis_cache)
    return 0


def cmd_tendl(args) -> int:
    output = args.output or ROOT / "data" / f"tendl{args.year}.reac"
    out = args.processed or DEFAULT_OUT / f"tendl{args.year}"
    if _use_existing(output, args, out):
        return 0
    folder = Path(args.raw) / f"TENDL-{args.year}-n"
    archives = sorted(Path(args.raw).glob(f"TENDL-{args.year}-n.t*"))
    archives = [a for a in archives if not a.name.endswith(".part")]
    if not archives and not args.url and folder.is_dir() and endf.find_endf_files(folder):
        # unpacked before (process.py fetch), and the archive is not there: no new download
        print(f"the TENDL-{args.year} files are in {folder}, they are not downloaded again", file=sys.stderr)
        source = folder
    else:
        if not archives:
            print(f"the raw TENDL-{args.year} files are not in {args.raw} (no TENDL-{args.year}-n.tgz and no folder "
                  f"TENDL-{args.year}-n/), so they are downloaded. The processed files in {out} are not enough: "
                  f"to make the data again, the raw files are processed again. To use the reaction data file that you "
                  f"have, stop now (Ctrl-C) and give --allow-old (make: ALLOW_OLD=1), or give --raw if the raw files "
                  f"are in another folder", file=sys.stderr, flush=True)
        # the files are read from the archive, so the 14 GB of TENDL are never all on the disk
        source = fetch.fetch_tendl(args.year, args.raw, args.url, unpack=False)
    status = _report(_process([source], out, args), out)
    return _finish(out, output, f"TENDL-{args.year} Incident Neutron File", status, args.nis_cache)


def cmd_libraries(args) -> int:
    for lib in libraries.LIBRARIES.values():
        print(f"{lib.key:12s} {lib.title:14s} {', '.join(lib.parts()):26s} licence: {lib.licence}")
    return 0


def _make(args, part: str) -> int:
    for key in args.libraries:
        libraries.make(libraries.get(key), part, args.raw, args.data, args.force)
    return 0


def cmd_decay(args) -> int:
    return _make(args, "decay")


def cmd_yields(args) -> int:
    for key in args.libraries:
        lib = libraries.get(key)
        for part in ("nfy", "sfy"):
            if getattr(lib, part):
                libraries.make(lib, part, args.raw, args.data, args.force)
    return 0


def _reactions(key: str, args) -> int:
    lib = libraries.get(key)
    output = args.output or Path(args.data) / f"{lib.key}.reac"
    out = args.processed or DEFAULT_OUT / lib.key
    if _use_existing(output, args, out):
        return 0
    folder = libraries.fetch_neutron(lib, args.raw)
    status = _report(_process([folder], out, args), out)
    return _finish(out, output, f"{lib.title} Incident Neutron File", status, args.nis_cache)


def cmd_reactions(args) -> int:
    return _reactions(args.library, args)


def cmd_bootstrap(args) -> int:
    for key in ("jeff311", "endfbvii1"):
        libraries.make(libraries.get(key), "decay", args.raw, args.data, args.force)
    for part in ("nfy", "sfy"):
        libraries.make(libraries.get("jeff311"), part, args.raw, args.data, args.force)
    if args.no_reactions:
        return 0
    if not (args.njoy or args.prepro_dir or os.environ.get("NJOY_EXE")):
        tools = build_tools.default_tools(args.tools)
        if not Path(tools.njoy).is_file():
            print("no NJOY and PREPRO, so build them first", file=sys.stderr)
            build_tools.build(args.tools, args.jobs, None, None, "gfortran", "-O2")
    return _reactions("jeff4", args)


def _add_data_options(p, force: bool = True) -> None:
    p.add_argument("--data", type=Path, default=DEFAULT_DATA, help="where to put the NIS data files (default data/)")
    p.add_argument("--raw", type=Path, default=DEFAULT_RAW, help="where to put the downloads")
    if force:  # else --force of the run options
        p.add_argument("--force", action="store_true", help="make the files again, also if they are there")


def _add_allow_old(p) -> None:
    p.add_argument("--allow-old", action="store_true",
                   help="use a reaction data file of an older version of the processing as it is (no MT=5 products)")


def _add_reaction_options(p) -> None:
    _add_allow_old(p)
    p.add_argument("--processed", type=Path, default=None, help="default processing/out/<library>")
    p.add_argument("--output", type=Path, default=None, help="default data/<library>.reac")
    p.add_argument("--nis-cache", type=Path, default=DEFAULT_NIS_CACHE)


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

    p = commands.add_parser("libraries", help="the libraries that can be downloaded")
    p.set_defaults(func=cmd_libraries)

    for name, func, text in (("decay", cmd_decay, "download decay data: data/<library>.dec"),
                             ("yields", cmd_yields, "download fission yields: data/<library>.nfy and .sfy")):
        p = commands.add_parser(name, help=text)
        p.add_argument("libraries", nargs="+", choices=list(libraries.LIBRARIES), metavar="library",
                       help=f"one or more of: {', '.join(libraries.LIBRARIES)}")
        _add_data_options(p)
        p.set_defaults(func=func)

    p = commands.add_parser("reactions", help="download the neutron files and make data/<library>.reac")
    p.add_argument("library", choices=[k for k, v in libraries.LIBRARIES.items() if v.neutron])
    _add_data_options(p, force=False)
    _add_reaction_options(p)
    _add_run_options(p)
    p.set_defaults(func=cmd_reactions)

    p = commands.add_parser("bootstrap", help="get all the data that the tests and the examples need")
    p.add_argument("--no-reactions", action="store_true",
                   help="only the decay data and the yields (minutes), not data/jeff4.reac (hours)")
    _add_data_options(p, force=False)
    _add_reaction_options(p)
    _add_run_options(p)
    p.set_defaults(func=cmd_bootstrap)

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
    p.add_argument("--title", default=None, help="the first line of the file")
    p.add_argument("--nis-cache", type=Path, default=DEFAULT_NIS_CACHE, help="the nis_cache program")
    p.set_defaults(func=cmd_combine)

    p = commands.add_parser("tendl", help="fetch, run and combine a TENDL release")
    p.add_argument("--year", default=fetch.LATEST_TENDL)
    p.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    p.add_argument("--url", default=None)
    p.add_argument("--processed", type=Path, default=None, help="default processing/out/tendl<year>")
    p.add_argument("--output", type=Path, default=None, help="default data/tendl<year>.reac")
    p.add_argument("--nis-cache", type=Path, default=DEFAULT_NIS_CACHE)
    _add_allow_old(p)
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
