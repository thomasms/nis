"""
Read the few things the processing needs from an ENDF file, as text (fast, and no strict parser
that stops on small format errors), and join the processed files into one tape.
"""

import re
from dataclasses import dataclass
from pathlib import Path

# NSUB of the incident neutron sublibrary
NEUTRON_SUBLIBRARY = 10

ELEMENTS = (
    "n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr "
    "Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W "
    "Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds "
    "Rg Cn Nh Fl Mc Lv Ts Og"
).split()

_FLOAT = re.compile(r"^([+-]?[0-9.]+)([+-]\d+)$")


def endf_float(text: str) -> float:
    """An ENDF number: 1.803700+4, 1.8037E+4 or 18037.0."""
    text = text.strip()
    if not text:
        return 0.0
    match = _FLOAT.match(text)
    if match:
        return float(f"{match.group(1)}e{match.group(2)}")
    return float(text.replace("D", "E").replace("d", "e"))


def endf_int(text: str) -> int:
    text = text.strip()
    return int(text) if text else 0


def mat_mf_mt(line: str) -> tuple:
    return endf_int(line[66:70]), endf_int(line[70:72]), endf_int(line[72:75])


@dataclass(frozen=True)
class Material:
    """The header of the first material of an ENDF file."""
    mat: int
    za: int
    liso: int
    nsub: int

    @property
    def z(self) -> int:
        return self.za // 1000

    @property
    def a(self) -> int:
        return self.za % 1000

    @property
    def zai(self) -> int:
        return self.za * 10 + self.liso

    @property
    def name(self) -> str:
        """E.g. Fe56, Co58m, Ag110m, Hf178n (the same names as in the NIS inputs)."""
        symbol = ELEMENTS[self.z] if self.z < len(ELEMENTS) else f"Z{self.z}"
        state = "" if self.liso == 0 else "mnopqrstuvwxyz"[self.liso - 1]
        return f"{symbol}{self.a}{state}"


def read_material(path: Path) -> Material:
    """The MAT, ZA, LISO and NSUB of the first material (MF1/MT451). ValueError if it is not ENDF."""
    with open(path, errors="replace") as f:
        lines = [f.readline() for _ in range(4)]
    # the first line is the tape identification (TPID)
    head = lines[1:]
    if any(len(line.rstrip("\n")) < 75 or mat_mf_mt(line)[1:] != (1, 451) for line in head):
        raise ValueError(f"{path} is not an ENDF file (no MF1/MT451 at the start)")
    mat = mat_mf_mt(head[0])[0]
    za = round(endf_float(head[0][0:11]))
    liso = endf_int(head[1][33:44])
    nsub = endf_int(head[2][44:55])
    return Material(mat=mat, za=za, liso=liso, nsub=nsub)


def find_endf_files(root: Path, sublibrary: int = NEUTRON_SUBLIBRARY) -> list:
    """All the ENDF files under root (in all folders) of one sublibrary, by default incident neutrons."""
    files = []
    for path in sorted(p for p in Path(root).rglob("*") if p.is_file()):
        try:
            material = read_material(path)
        except (ValueError, OSError, UnicodeDecodeError):
            continue
        if material.nsub == sublibrary:
            files.append(path)
    return files


def materials_text(path: Path) -> str:
    """The materials of a tape without its first line (TPID) and its end line (TEND, MAT = -1)."""
    lines = Path(path).read_text().splitlines(keepends=True)
    body = lines[1:]
    while body and (not body[-1].strip() or mat_mf_mt(body[-1])[0] == -1):
        body.pop()
    return "".join(body)


def tend_line() -> str:
    return f"{'':66s}{-1:4d}{0:2d}{0:3d}{0:5d}\n"


def tpid_line(title: str) -> str:
    return f"{title[:66]:66s}{1:4d}{0:2d}{0:3d}{0:5d}\n"


def combine(files: list, output: Path, title: str) -> int:
    """Join one material files into one tape, ordered by ZA then LISO. Returns the number of materials."""
    ordered = sorted(files, key=lambda p: read_material(p).zai)
    with open(output, "w") as out:
        out.write(tpid_line(title))
        for path in ordered:
            out.write(materials_text(path))
        out.write(tend_line())
    return len(ordered)
