"""
Converters: each one reads the original files of a benchmark and writes the standard format.

A converter is a module here with
    NAME           the benchmark name (the folder in benchmarks/)
    METADATA       the benchmark.json content (without format and version)
    SOURCE_HELP    what the original files are, for the --source help
    convert(source: Path, only: list) -> list of Case (without path)
"""

import importlib
import pkgutil


def converters() -> dict:
    out = {}
    for module in pkgutil.iter_modules(__path__):
        m = importlib.import_module(f"{__name__}.{module.name}")
        out[m.NAME] = m
    return out
