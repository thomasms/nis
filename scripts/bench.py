"""
NIS benchmarks: one set of tools for all the benchmarks in the standard format (benchmarks/README.md).

Commands:
    list       the benchmarks and their cases
    check      check the case files, and that each deck (<case>.inp) is made from its case file
    convert    make a benchmark in the standard format from its original files
    decks      make the NIS decks from the case files again
    run        run the cases with NIS, write the results (one CSV for all benchmarks)
    stats      C/E statistics of a results file
    plot       plots of a results file: a summary for each benchmark, or one plot for each case

Examples:
    uv run python scripts/bench.py list
    uv run python scripts/bench.py run fns --decay-data data/jeff311.dec --reaction-data data/jeff4.reac \\
        --output results.csv
    uv run python scripts/bench.py stats results.csv
    uv run python scripts/bench.py plot results.csv --summary fns_summary.png
"""

import argparse
import os
import sys
import time
from pathlib import Path

from nisbench import FormatError, load_benchmarks
from nisbench.format import BENCHMARKS_DIR, load_case, write_benchmark, write_case


def _add_selection(parser) -> None:
    parser.add_argument("benchmarks", nargs="*", help="benchmark names, all if none")
    parser.add_argument("--cases", nargs="*", default=[],
                        help="only these cases, or the cases that start with them (e.g. Co or Co/2000exp_5min)")
    parser.add_argument("--root", type=Path, default=BENCHMARKS_DIR, help="the benchmarks folder")


def cmd_list(args) -> int:
    for b in load_benchmarks(args.benchmarks, args.root, args.cases):
        points = sum(len(m.values) for c in b.cases for m in c.measurements)
        quantities = sorted({m.quantity for c in b.cases for m in c.measurements})
        print(f"{b.name:10s} {b.category:8s} {len(b.cases):4d} cases {points:6d} points  "
              f"{', '.join(quantities)}  {b.title}")
        if args.cases or args.verbose:
            for c in b.cases:
                print(f"    {c.name:28s} {c.title}")
    return 0


def cmd_check(args) -> int:
    from nisbench.deck import deck_is_current
    errors = 0
    for b in load_benchmarks(args.benchmarks, args.root, args.cases):
        stale = [c for c in b.cases if not deck_is_current(c)]
        for c in stale:
            print(f"{c.deck_path}: missing or not made from {c.path.name}, run: scripts/bench.py decks {b.name}")
        errors += len(stale)
        print(f"{b.name}: {len(b.cases)} cases, {len(b.cases) - len(stale)} decks up to date", file=sys.stderr)
    return 1 if errors else 0


def cmd_convert(args) -> int:
    from nisbench.convert import converters
    from nisbench.deck import write_deck
    available = converters()
    if args.benchmark not in available:
        raise SystemExit(f"no converter for {args.benchmark!r}, known: {', '.join(available)}")
    converter = available[args.benchmark]
    folder = args.root / converter.NAME
    cases = converter.convert(args.source, args.only)
    if not cases:
        raise SystemExit(f"no cases in {args.source}, it must be {converter.SOURCE_HELP}")
    write_benchmark(converter.METADATA, folder)
    for case in cases:
        path = folder / f"{case.name}.json"
        write_case(case, path)
        write_deck(load_case(path))
    print(f"{len(cases)} cases in {folder}", file=sys.stderr)
    return 0


def cmd_decks(args) -> int:
    from nisbench.deck import write_deck
    for b in load_benchmarks(args.benchmarks, args.root, args.cases):
        for c in b.cases:
            write_deck(c)
        print(f"{b.name}: {len(b.cases)} decks", file=sys.stderr)
    return 0


def cmd_run(args) -> int:
    from nisbench.engines import NuclearData, reaction_cache, run
    from nisbench.results import calculation_names, case_rows, write_results
    benchmarks = load_benchmarks(args.benchmarks, args.root, args.cases)
    cases = [c for b in benchmarks for c in b.cases]
    if not cases:
        raise SystemExit("no cases")
    print(f"{len(cases)} cases", file=sys.stderr)

    reaction = None
    if args.reaction_data:
        nis_cache = args.nis_cache or (args.nis.parent / "nis_cache" if args.nis else None)
        reaction = reaction_cache(args.reaction_data, nis_cache)
    data = NuclearData(decay=args.decay_data, reaction=reaction, fission_yields=args.fission_yields,
                       sf_yields=args.sf_yields)

    start = time.time()
    rows = []
    for case, results in zip(cases, run(cases, data, args.engine, args.nis, args.jobs)):
        case_result = case_rows(case, results)
        rows += case_result
        if args.verbose:
            ce = [float(r["nis"]) / float(r["experiment"]) for r in case_result if float(r["experiment"]) > 0]
            print(f"{case.benchmark} {case.name:28s} NIS C/E {sum(ce) / len(ce):7.3f}", file=sys.stderr)
    write_results(args.output, rows, calculation_names(cases))
    print(f"done in {time.time() - start:.1f} s, results in {args.output}", file=sys.stderr)
    return 0


def cmd_stats(args) -> int:
    from nisbench.results import read_results
    from nisbench.stats import report
    rows, calculations = read_results(args.results)
    print(report(rows, calculations, args.worst))
    return 0


def cmd_plot(args) -> int:
    from nisbench import plot
    from nisbench.results import read_results
    rows, calculations = read_results(args.results)
    pairs = sorted({(r["benchmark"], r["quantity"]) for r in rows})
    if args.benchmark:
        pairs = [p for p in pairs if p[0] == args.benchmark]
    if args.quantity:
        pairs = [p for p in pairs if p[1] == args.quantity]
    if not pairs:
        raise SystemExit("nothing to plot")
    if args.summary:
        if len(pairs) > 1:
            raise SystemExit(f"the results have {len(pairs)} benchmarks and quantities, give --benchmark and --quantity")
        benchmark, quantity = pairs[0]
        title = None
        if (args.root / benchmark / "benchmark.json").exists():
            import json
            name = json.loads((args.root / benchmark / "benchmark.json").read_text()).get("title", benchmark)
            title = f"{name}: C/E for each case"
        plot.summary(rows, calculations, benchmark, quantity, args.summary, columns=args.columns, title=title)
        print(f"plot in {args.summary}", file=sys.stderr)
    if args.cases_dir:
        selected = [r for r in rows if (r["benchmark"], r["quantity"]) in pairs]
        files = plot.cases(selected, calculations, args.cases_dir, args.format)
        print(f"{len(files)} plots in {args.cases_dir}", file=sys.stderr)
    if not args.summary and not args.cases_dir:
        raise SystemExit("give --summary and/or --cases-dir")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    p = commands.add_parser("list", help="the benchmarks and their cases")
    _add_selection(p)
    p.add_argument("-v", "--verbose", action="store_true", help="list the cases too")
    p.set_defaults(func=cmd_list)

    p = commands.add_parser("check", help="check the case files and the decks")
    _add_selection(p)
    p.set_defaults(func=cmd_check)

    p = commands.add_parser("convert", help="make a benchmark from its original files")
    p.add_argument("benchmark", help="the converter, e.g. fns")
    p.add_argument("--source", type=Path, required=True, help="the folder of the original files")
    p.add_argument("--only", nargs="*", default=[], help="only some cases (the converter says how to name them)")
    p.add_argument("--root", type=Path, default=BENCHMARKS_DIR, help="the benchmarks folder")
    p.set_defaults(func=cmd_convert)

    p = commands.add_parser("decks", help="make the NIS decks from the case files again")
    _add_selection(p)
    p.set_defaults(func=cmd_decks)

    p = commands.add_parser("run", help="run the cases with NIS")
    _add_selection(p)
    p.add_argument("--decay-data", type=Path, required=True)
    p.add_argument("--reaction-data", type=Path, default=None,
                   help="ENDF reaction data or a cache of it (the cache is made first if needed)")
    p.add_argument("--fission-yields", type=Path, default=None, help="neutron induced fission yields (ENDF)")
    p.add_argument("--sf-yields", type=Path, default=None, help="spontaneous fission yields (ENDF)")
    p.add_argument("--engine", choices=("pynis", "nis"), default="pynis",
                   help="pynis: the python bindings (fast), nis: the nis program with each deck")
    p.add_argument("--nis", type=Path, default=None, help="the nis program, for --engine nis")
    p.add_argument("--nis-cache", type=Path, default=None, help="the nis_cache program, by default next to --nis")
    p.add_argument("--output", type=Path, default=Path("results.csv"))
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="number of cases to run at the same time")
    p.add_argument("-v", "--verbose", action="store_true", help="the mean C/E of each case")
    p.set_defaults(func=cmd_run)

    p = commands.add_parser("stats", help="C/E statistics of results files")
    p.add_argument("results", nargs="+", type=Path)
    p.add_argument("--worst", type=int, default=25, help="number of cases to list where NIS differs most")
    p.set_defaults(func=cmd_stats)

    p = commands.add_parser("plot", help="plots of results files")
    p.add_argument("results", nargs="+", type=Path)
    p.add_argument("--benchmark", default=None)
    p.add_argument("--quantity", default=None)
    p.add_argument("--summary", type=Path, default=None, help="the file of the summary plot")
    p.add_argument("--columns", type=int, default=2, help="columns of the summary plot")
    p.add_argument("--cases-dir", type=Path, default=None, help="a folder for one plot for each case")
    p.add_argument("--format", default="png", help="file format of the case plots, e.g. png, pdf, svg")
    p.add_argument("--root", type=Path, default=BENCHMARKS_DIR, help="the benchmarks folder (for the titles)")
    p.set_defaults(func=cmd_plot)

    args = parser.parse_args()
    try:
        return args.func(args)
    except FormatError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
