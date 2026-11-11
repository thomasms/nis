"""Paths to the nuclear data, the benchmarks and the programs, for the python tests."""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent.parent
DATA_DIR = ROOT / "data"
FNS_DECKS_DIR = ROOT / "benchmarks" / "fns"
# the NIS results of all the FNS decks, the FNS regression test compares with it
FNS_REFERENCE = FNS_DECKS_DIR / "nis_reference.csv"

JEFF311_DECAY = DATA_DIR / "jeff311.dec"
ENDFBVII_DECAY = DATA_DIR / "endfbvii1.dec"
JEFF4_REACTION = DATA_DIR / "jeff4.reac"
# fission yields, from the NEA data bank (see the README)
JEFF311_NFY = DATA_DIR / "jeff311.nfy"
JEFF311_SFY = DATA_DIR / "jeff311.sfy"

# the nis program, from the CMake build (or give its path in NIS_BIN)
NIS_BIN = Path(os.environ.get("NIS_BIN", ROOT / "build" / "bin" / "nis"))

# the benchmark scripts are not a package, so add their folder to the path
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


def has_data(path: Path, min_size: int = 100_000) -> bool:
    """True if the data file is there. Without git lfs, the file is a small text pointer file."""
    return path.exists() and path.stat().st_size > min_size
