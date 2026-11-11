"""Process many ENDF files at the same time, one nuclide for each worker."""

import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

from .endf import read_material
from .pipeline import Result, Tools, process_file


def process_all(files: list, output_dir: Path, tools: Tools, jobs: int = None, temperature: float = 293.6,
                work_root: Path = None, keep: bool = False, force: bool = False, quiet: bool = False) -> list:
    """
    Process the files with `jobs` workers (default: the number of cores). Each worker runs the programs of
    one nuclide, so the work is in the NJOY and PREPRO processes and threads are enough here. The largest
    files start first, so that no long nuclide is left at the end to run alone.
    """
    tools.check()
    check_unique(files)
    jobs = jobs or os.cpu_count() or 1
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ordered = sorted(files, key=lambda p: Path(p).stat().st_size, reverse=True)

    results = []
    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(process_file, f, output_dir, tools, temperature, work_root, keep, force)
                   for f in ordered]
        for i, future in enumerate(as_completed(futures), 1):
            r = future.result()
            results.append(r)
            if not quiet:
                state = "skip (done)" if r.skipped else ("ok" if r.ok else f"FAILED: {r.error}")
                print(f"[{i}/{len(ordered)}] {r.name or r.source.name:8s} {r.seconds:7.1f} s  {state}",
                      file=sys.stderr, flush=True)
    elapsed = time.perf_counter() - start

    write_summary(output_dir / "summary.json", results, elapsed, jobs)
    return results


def check_unique(files: list) -> None:
    """Each nuclide only one time: the output file is named after the nuclide."""
    by_name = {}
    for f in files:
        try:
            by_name.setdefault(read_material(f).name, []).append(str(f))
        except (ValueError, OSError):
            # process_file reports it
            pass
    duplicates = {name: paths for name, paths in by_name.items() if len(paths) > 1}
    if duplicates:
        text = "\n".join(f"    {name}: {', '.join(paths)}" for name, paths in sorted(duplicates.items()))
        raise ValueError(f"more than one file for the same nuclide, give only one:\n{text}")


def write_summary(path: Path, results: list, elapsed: float, jobs: int) -> None:
    def entry(r: Result) -> dict:
        d = asdict(r)
        d["source"] = str(r.source)
        d["output"] = str(r.output) if r.output else None
        return d
    failed = [r for r in results if not r.ok]
    summary = {"files": len(results), "processed": sum(r.ok and not r.skipped for r in results),
               "skipped": sum(r.skipped for r in results), "failed": len(failed),
               "jobs": jobs, "seconds": round(elapsed, 1),
               "results": sorted((entry(r) for r in results), key=lambda d: d["name"] or d["source"])}
    Path(path).write_text(json.dumps(summary, indent=1) + "\n")
