"""Process many ENDF files at the same time, one nuclide for each worker."""

import json
import os
import shutil
import sys
import tarfile
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

from .endf import NEUTRON_SUBLIBRARY, read_material
from .pipeline import Result, Tools, default_work_root, process_file


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


def is_archive(path: Path) -> bool:
    return Path(path).is_file() and tarfile.is_tarfile(path)


def process_archive(archive: Path, output_dir: Path, tools: Tools, jobs: int = None, temperature: float = 293.6,
                    work_root: Path = None, keep: bool = False, force: bool = False, quiet: bool = False,
                    only: tuple = ()) -> list:
    """
    Process the ENDF neutron files in a tar archive (e.g. TENDL-n.tgz) without unpacking all of it: the archive
    is read one time from start to end, each file is written to the work folder just before it is processed
    and deleted after. Only a few files (2 for each worker) are out of the archive at one time, so a library
    of 14 GB needs no more disk than its archive and the outputs. `only`: nuclide names, all if empty.
    """
    tools.check()
    jobs = jobs or os.cpu_count() or 1
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="nisdata_archive_", dir=work_root or default_work_root()))
    slots = threading.BoundedSemaphore(2 * jobs)
    lock = threading.Lock()
    results, seen = [], {}

    def report(r: Result) -> None:
        with lock:
            results.append(r)
            if not quiet:
                state = "skip (done)" if r.skipped else ("ok" if r.ok else f"FAILED: {r.error}")
                print(f"[{len(results)}] {r.name or r.source.name:8s} {r.seconds:7.1f} s  {state}",
                      file=sys.stderr, flush=True)

    def work(staged: Path, member: str) -> None:
        try:
            r = process_file(staged, output_dir, tools, temperature, work_root, keep, force)
            r.source = Path(f"{archive}:{member}")
            report(r)
        finally:
            staged.unlink(missing_ok=True)
            slots.release()

    start = time.perf_counter()
    try:
        with ThreadPoolExecutor(max_workers=jobs) as pool, tarfile.open(archive, "r|*") as tar:
            futures = []
            for member in tar:
                if not member.isfile():
                    continue
                slots.acquire()
                staged = staging / f"{len(futures)}_{Path(member.name).name}"
                with tar.extractfile(member) as source, open(staged, "wb") as target:
                    shutil.copyfileobj(source, target, length=1 << 20)
                try:
                    material = read_material(staged)
                except (ValueError, OSError):
                    material = None
                if material is None or material.nsub != NEUTRON_SUBLIBRARY or (only and material.name not in only):
                    staged.unlink()
                    slots.release()
                    continue
                if material.name in seen:
                    staged.unlink()
                    slots.release()
                    report(Result(source=Path(f"{archive}:{member.name}"), name=material.name,
                                  error=f"the same nuclide as {seen[material.name]}, not processed"))
                    continue
                seen[material.name] = member.name
                futures.append(pool.submit(work, staged, member.name))
            for future in futures:
                future.result()
    finally:
        shutil.rmtree(staging, ignore_errors=True)
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
