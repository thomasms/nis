"""
nisbench: the standard benchmark format of NIS, and the tools to use it.

Each benchmark (for example fns) is a folder in benchmarks/ with a benchmark.json and one
<case>.json for each case. The case file is the source of truth: the NIS input deck <case>.inp
next to it is made from it. The format is in benchmarks/README.md.

The command line is scripts/bench.py.
"""

from .format import (BENCHMARKS_DIR, Benchmark, Case, FormatError, Measurement, load_benchmark,
                     load_benchmarks, load_case, read_group_structure)

__all__ = ["BENCHMARKS_DIR", "Benchmark", "Case", "FormatError", "Measurement", "load_benchmark",
           "load_benchmarks", "load_case", "read_group_structure"]
