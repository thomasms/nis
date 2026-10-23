"""
Write the ALL_ZAIS list in src/nis/nuclides.hpp from the current list and ENDF decay files.

Every nuclide (ground state and isomers) with an MF=8 MT=457 section in the decay files
is added, so decays to isomers are not lost to the sink.

Usage:
    python scripts/make_zai_list.py data/jeff311.dec data/endfbvii1.dec
"""

import argparse
import re
from pathlib import Path

_NUCLIDES_FILE = Path(__file__).parent.parent / "src" / "nis" / "nuclides.hpp"
_ARRAY = re.compile(
    r"(constexpr const int MAX_SIZE = )\d+(;.*?constexpr std::array<int, MAX_SIZE> ALL_ZAIS = \{\n)(.*?)(\n\s*\};)",
    re.DOTALL,
)


def endf_float(text: str) -> float:
    # ENDF writes 1.002300+3 for 1.0023e+3
    return float(re.sub(r"(?<=\d)([+-])(\d+)$", r"e\1\2", text.strip()))


def decay_zais(path: Path) -> set:
    zais = set()
    previous = None
    with open(path) as f:
        for line in f:
            key = line[66:75]
            # the first line of each MF=8 MT=457 section: ZA, AWR, LIS, LISO, NST, NSP
            if key[4:] == " 8457" and key != previous:
                za = int(round(endf_float(line[0:11])))
                liso = int(line[33:44])
                zais.add(za * 10 + liso)
            previous = key
    return zais


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("decay_files", nargs="+", type=Path)
    args = parser.parse_args()

    text = _NUCLIDES_FILE.read_text()
    match = _ARRAY.search(text)
    current = {int(v) for v in re.findall(r"\d+", match.group(3))}

    zais = set(current)
    for path in args.decay_files:
        zais |= decay_zais(path)
    # the free neutron (Z=0) has no element name, and is not needed
    zais = sorted(z for z in zais if z >= 10000)

    indent = " " * 12
    lines = []
    for i in range(0, len(zais), 14):
        lines.append(indent + ", ".join(str(z) for z in zais[i:i + 14]) + ",")
    body = "\n".join(lines).rstrip(",")

    text = text[:match.start()] + match.group(1) + str(len(zais)) + match.group(2) + body + match.group(4) + text[match.end():]
    _NUCLIDES_FILE.write_text(text)
    print(f"{len(zais)} nuclides, was {len(current)}, {sum(1 for z in zais if z % 10)} isomers")


if __name__ == "__main__":
    main()
