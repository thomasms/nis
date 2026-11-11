"""
Download the JEFF-3.1.1 fission yields from the NEA data bank (data.oecd-nea.org) and check their MD5:

    data/jeff311.nfy   neutron induced fission yields (19 fissile nuclides)
    data/jeff311.sfy   spontaneous fission yields (Cm242, Cm244, Cf252)

Needs curl. Usage:
    python scripts/get_fission_yields.py [--output-dir data]
"""

import argparse
import hashlib
import sys
import subprocess
from pathlib import Path

# the NEA data bank record "JEFF-3.1.1 Evaluated Data: fission product yields data"
_RECORD = "https://data.oecd-nea.org/api/records/yv03w-dxz31/files"
_FILES = {
    "jeff311.nfy": ("JEFF311NFY.ASC", "e0b8eb05e3a8c08a4ca3bda4d39fc262"),
    "jeff311.sfy": ("JEFF311SFY.ASC", "03d496b0c6a52d5290ffb07cf7c47a53"),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path, default=Path("data"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for name, (key, md5) in _FILES.items():
        path = args.output_dir / name
        if path.exists() and hashlib.md5(path.read_bytes()).hexdigest() == md5:
            print(f"{path} is already there", file=sys.stderr)
            continue
        # curl, as the data bank refuses the python downloader
        content = subprocess.run(["curl", "-sSfL", "--max-time", "300", f"{_RECORD}/{key}/content"],
                                 capture_output=True, check=True).stdout
        if hashlib.md5(content).hexdigest() != md5:
            print(f"{key}: the MD5 is not correct, the file is not written", file=sys.stderr)
            return 1
        path.write_bytes(content)
        print(f"{path}: {len(content)} bytes", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
