"""
IRDFF-II spectrum-averaged cross sections: from the IRDFF-II files to the standard format.

In a standard neutron spectrum (Cf-252 spontaneous fission, U-235 thermal fission), the measured spectrum-averaged
cross section (SACS) of a reaction is the production rate of the product over the flux. NIS gets it from its
reaction matrix of the same spectrum, so each case tests the reaction data of one reaction directly (and its isomer
data: In115(n,n')In115m, Nb93(n,2n)Nb92m, ...).

The original files, all in one folder (<source>):

    irdff_sp_<name>.g       (or .endf)  a spectrum of the IRDFF-II spectra page (nds.iaea.org/IRDFF/spectra.html):
                            MF3/MT261, 725 groups (histogram) or pointwise (lin-lin), in 1/eV. The MAT gives the
                            spectrum: 9861 Cf-252(s.f.), 9228 U-235(nth,f). With both forms, the 725 groups are used:
                            the case files then give only the group values (group structure IRDFF-725), and the
                            cross sections are the same as with the pointwise form (to 4 digits).
    irdff_sp_<name>_RR_UNC.pdf (or .lst, .txt)  optional: the SACS that IRDFF-II calculates in the spectrum (the
                            calculation "IRDFF-II"). Read with pdftotext (poppler) for a .pdf.

The measured SACS are in measured.csv: <source>/measured.csv if it is there, else benchmarks/irdff2/measured.csv
(the values of the IRDFF-II paper, in the repository). Its columns:

    spectrum,table,notation,target,product,value_mb,uncertainty_percent,reference

spectrum is Cf252 or U235, target a nuclide (Al27), product a nuclide in all its states (Na24), or only its ground
state (Al26g), or an isomer (Al26m). table, notation (the IRDFF-II name of the reaction, Al27a) and reference are
only to find the value in the paper. Lines that start with # are comments.

Each case is one reaction in one spectrum: 1 atom of the target, one step with the spectrum, and the cross section
of the product at the end of the step. The case name is <spectrum>/<target>-<product>, e.g. Cf252/Al27-Na24.
"""

import csv
import math
import re
import shutil
import subprocess
from pathlib import Path

from ..format import BENCHMARKS_DIR, GROUP_STRUCTURES_DIR, Case, Measurement, read_group_structure

NAME = "irdff2"
METADATA = {
    "name": NAME,
    "title": "IRDFF-II spectrum-averaged cross sections",
    "category": "other",
    "description": "Measured spectrum-averaged cross sections of activation (dosimetry) reactions in the Cf-252 "
                   "spontaneous fission and U-235 thermal fission neutron spectra, with the IRDFF-II spectra. "
                   "Each case is one reaction in one spectrum. It tests the reaction data of the reaction (and its "
                   "isomer data) in a well known fast neutron spectrum.",
    "references": [
        "A. Trkov et al., IRDFF-II: A New Neutron Metrology Library, Nuclear Data Sheets 163 (2020) 1-108",
        "W. Mannhart, Status of the Cf-252 fission neutron spectrum evaluation with regard to recent experiments, "
        "INDC(NDS)-220 (1989) 305",
    ],
    "source": "The IRDFF-II spectra (nds.iaea.org/IRDFF), the measured values of the IRDFF-II paper (measured.csv), "
              "converted with scripts/bench.py convert irdff2. IRDFF-II is the cross section that IRDFF-II calculates "
              "in the spectrum (its RR_UNC listing).",
}
SOURCE_HELP = ("a folder with the IRDFF-II spectra (irdff_sp_*.endf or .g), measured.csv and optionally the "
               "RR_UNC listings (irdff_sp_*_RR_UNC.pdf)")

CALCULATION = "IRDFF-II"
# MAT of the IRDFF-II spectrum -> (name in the case names and measured.csv, title)
SPECTRA = {9861: ("Cf252", "Cf-252 spontaneous fission"), 9228: ("U235", "U-235 thermal fission")}

ELEMENTS = (
    "n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr "
    "Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W "
    "Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf"
).split()

# MT -> (change of Z, change of A) of the residual, for the reactions of the RR_UNC listings that make a nuclide
_MT_RESIDUAL = {16: (0, -1), 17: (0, -2), 22: (-2, -4), 28: (-1, -1), 37: (0, -3), 102: (0, 1), 103: (-1, 0),
                104: (-1, -1), 105: (-1, -2), 106: (-2, -2), 107: (-2, -3), 152: (0, -4)}
# MT -> the light nuclide of a gas production (all the reactions that make it). A light residual of one reaction,
# such as the He4 of Li6(n,t), is not used: it is only a part of the production.
_GAS_PRODUCTION = {203: "H1", 204: "H2", 205: "H3", 206: "He3", 207: "He4"}

MEASURED = BENCHMARKS_DIR / NAME / "measured.csv"

# the time step: the cross section does not depend on the length or the flux, so 1 s and 1 n/m2/s
_STEP = (1.0, 1.0)


# ---------------------------------------------------------------- nuclide names

def nuclide_name(za: int) -> str:
    return f"{ELEMENTS[za // 1000]}{za % 1000}"


_PRODUCT = re.compile(r"^([A-Z][a-z]?\d{1,3})([gmn]?)$")


def product_state(token: str) -> tuple:
    """(product, all_states) of a product of measured.csv: Na24 is all the states, Al26g only the ground state,
    Al26m the isomer."""
    match = _PRODUCT.match(token)
    if not match:
        raise ValueError(f"{token!r} is not a product such as Na24, Al26g or Al26m")
    nuclide, state = match.groups()
    if state == "":
        return nuclide, True
    return (nuclide, False) if state == "g" else (nuclide + state, False)


def product_token(product: str, all_states: bool) -> str:
    """The inverse of product_state."""
    return product if all_states or product[-1] in "mn" else product + "g"


# ---------------------------------------------------------------- the spectra

def _float(text: str) -> float:
    """An ENDF number: 1.234567+5, 1.2345E-12, .001000000 or 0.0."""
    text = text.strip()
    if not text:
        return 0.0
    match = re.match(r"^([+-]?[0-9.]+)([+-]\d+)$", text)
    return float(f"{match.group(1)}e{match.group(2)}") if match else float(text)


def read_spectrum(path: Path) -> tuple:
    """(MAT, energies in eV, values in 1/eV, interpolation laws, breaks) of the MF3/MT261 TAB1 of an IRDFF-II file."""
    lines = [line for line in Path(path).read_text(errors="replace").splitlines()
             if len(line) >= 75 and line[70:75] == " 3261"]
    if not lines:
        raise ValueError(f"{path}: no MF3/MT261 (spectrum)")
    mat = int(lines[0][66:70])
    head = lines[1]
    nr, np_ = int(head[44:55]), int(head[55:66])

    def numbers(start: int, count: int) -> tuple:
        values, i = [], start
        while len(values) < count:
            values += [_float(lines[i][11 * k:11 * k + 11]) for k in range(6)]
            i += 1
        return values[:count], i

    interpolation, i = numbers(2, 2 * nr)
    xy, _ = numbers(i, 2 * np_)
    breaks = [int(v) for v in interpolation[0::2]]
    laws = [int(v) for v in interpolation[1::2]]
    return mat, xy[0::2], xy[1::2], laws, breaks


def group_spectrum(energies: list, values: list, laws: list, breaks: list) -> tuple:
    """(bounds, group values) of a TAB1 spectrum: the integral over each interval between two points, exact for a
    histogram (INT=1) or lin-lin (INT=2). Intervals of zero width, from 0 eV, or at the ends with a zero integral
    are left out (NIS needs bounds above 0)."""
    bounds, groups = [], []
    region = 0
    for i in range(len(energies) - 1):
        while i + 2 > breaks[region]:
            region += 1
        law = laws[region]
        if law not in (1, 2):
            raise ValueError(f"interpolation law {law} is not supported (only 1 and 2)")
        low, high = energies[i], energies[i + 1]
        if high <= low or low <= 0.0:
            continue
        integral = values[i] * (high - low) if law == 1 else 0.5 * (values[i] + values[i + 1]) * (high - low)
        if bounds and bounds[-1] != low:
            raise ValueError(f"the spectrum has a gap at {low} eV")
        if not bounds:
            bounds.append(low)
        bounds.append(high)
        groups.append(integral)
    while groups and groups[-1] == 0.0:
        groups.pop()
        bounds.pop()
    while groups and groups[0] == 0.0:
        groups.pop(0)
        bounds.pop(0)
    return bounds, groups


def _group_structure(bounds: list) -> str:
    """The name of the group structure in benchmarks/group_structures with these bounds, or None."""
    for path in sorted(GROUP_STRUCTURES_DIR.glob("*.txt")):
        known = read_group_structure(path.stem)
        if len(known) == len(bounds) and all(math.isclose(a, b, rel_tol=1e-6) for a, b in zip(known, bounds)):
            return path.stem
    return None


def case_spectrum(energies: list, values: list, laws: list, breaks: list) -> tuple:
    """(bounds, group values, group structure name) of a TAB1 spectrum. A histogram on the bounds of a known group
    structure gives all its groups, with the name. Else the groups of group_spectrum, with no name."""
    if laws == [1] and _group_structure(energies) is not None:
        groups = [v * (high - low) for v, low, high in zip(values, energies, energies[1:])]
        return energies, groups, _group_structure(energies)
    return (*group_spectrum(energies, values, laws, breaks), None)


def _spectrum_files(source: Path) -> dict:
    """MAT -> the file of the spectrum (the groups if there are both)."""
    found = {}
    for path in sorted(Path(source).glob("irdff_sp_*")):
        if path.suffix not in (".endf", ".g"):
            continue
        mat = read_spectrum(path)[0]
        if mat not in found or path.suffix == ".g":
            found[mat] = path
    return found


# ---------------------------------------------------------------- the IRDFF-II calculation (RR_UNC listing)

# "  57   41 13027 11024g 8.66807 1.0167E+00 +/- 1.804E-02 ..." (the first number is the line of the listing)
_RR_LINE = re.compile(r"^\s*\d+\s+\d+\s+(\d+)\s+(\d+|\d+[gm])\s+\S+\s+([0-9.]+E[+-]\d+)\s+\+/-\s+([0-9.]+E[+-]\d+)")


def _listing_text(path: Path) -> str:
    if path.suffix.lower() != ".pdf":
        return path.read_text(errors="replace")
    if shutil.which("pdftotext") is None:
        raise RuntimeError(f"{path}: pdftotext (poppler) is needed to read a .pdf, or give the listing as text")
    return subprocess.run(["pdftotext", "-layout", str(path), "-"], capture_output=True, text=True, check=True).stdout


def read_rr_unc(path: Path) -> tuple:
    """(MAT of the spectrum, {(target, product, all_states): SACS in mb}) of an RR_UNC listing. Only the reactions
    of one nuclide (not of an element, such as Fe) that make one nuclide."""
    text = _listing_text(Path(path))
    match = re.search(r"Spectrum MAT No\.\s*:\s*(\d+)", text)
    if not match:
        raise ValueError(f"{path}: not an RR_UNC listing (no spectrum MAT)")
    values = {}
    for line in text.splitlines():
        m = _RR_LINE.match(line)
        if not m:
            continue
        za, reaction, value = int(m.group(1)), m.group(2), float(m.group(3))
        if za % 1000 == 0:
            continue
        if reaction.isdigit() and int(reaction) in _GAS_PRODUCTION:
            key = (nuclide_name(za), _GAS_PRODUCTION[int(reaction)], True)
        elif reaction.isdigit():
            if int(reaction) not in _MT_RESIDUAL:
                continue
            dz, da = _MT_RESIDUAL[int(reaction)]
            residual = za + 1000 * dz + da
            if residual // 1000 <= 2:
                continue
            key = (nuclide_name(za), nuclide_name(residual), True)
        else:
            product, all_states = product_state(nuclide_name(int(reaction[:-1])) + reaction[-1])
            key = (nuclide_name(za), product, all_states)
        values[key] = value
    return int(match.group(1)), values


# ---------------------------------------------------------------- the measured values

def read_measured(path: Path) -> list:
    """The rows of measured.csv: (spectrum, target, product, all_states, value_mb, uncertainty_mb, reference)."""
    rows = []
    with open(path, newline="") as f:
        lines = [line for line in f if line.strip() and not line.lstrip().startswith("#")]
    for i, row in enumerate(csv.DictReader(lines), 2):
        try:
            product, all_states = product_state(row["product"].strip())
            value = float(row["value_mb"])
            rows.append((row["spectrum"].strip(), row["target"].strip(), product, all_states, value,
                         value * float(row["uncertainty_percent"]) / 100.0, (row.get("reference") or "").strip()))
        except (KeyError, ValueError, TypeError) as e:
            raise ValueError(f"{path}, row {i}: {e}") from e
    return rows


# ---------------------------------------------------------------- the cases

def make_case(spectrum: str, title: str, bounds: list, groups: list, target: str, product: str, all_states: bool,
              value: float, uncertainty: float, calculated: float = None, group_structure: str = None) -> Case:
    token = product_token(product, all_states)
    measurement = Measurement("cross_section", "mb", [0], [value], [uncertainty],
                              {CALCULATION: [calculated]} if calculated is not None else {},
                              target=target, product=product, all_states=all_states)
    return Case(benchmark=NAME, name=f"{spectrum}/{target}-{token}", title=f"{target} to {token} in the {title} spectrum",
                mass_kg=None, atoms={target: 1.0}, steps=[_STEP], measurements=[measurement],
                spectrum=(bounds, groups), group_structure=group_structure)


def convert(source: Path, only: list = None) -> list:
    """only: spectrum names (Cf252, U235) or case names (Cf252/Al27-Na24), all if empty."""
    source = Path(source)
    spectra = {}
    for mat, path in _spectrum_files(source).items():
        if mat in SPECTRA:
            _, energies, values, laws, breaks = read_spectrum(path)
            spectra[SPECTRA[mat][0]] = (SPECTRA[mat][1], *case_spectrum(energies, values, laws, breaks))
    if not spectra:
        return []

    calculated = {}
    for path in sorted(source.glob("irdff_sp_*_RR_UNC.*")):
        mat, values = read_rr_unc(path)
        if mat in SPECTRA:
            calculated[SPECTRA[mat][0]] = values

    measured_path = source / "measured.csv" if (source / "measured.csv").exists() else MEASURED
    cases = []
    for spectrum, target, product, all_states, value, uncertainty, _ in read_measured(measured_path):
        if spectrum not in spectra:
            raise ValueError(f"{measured_path}: no spectrum file for {spectrum}")
        title, bounds, groups, structure = spectra[spectrum]
        calc = calculated.get(spectrum, {}).get((target, product, all_states))
        case = make_case(spectrum, title, bounds, groups, target, product, all_states, value, uncertainty, calc,
                         structure)
        if not only or any(case.name == o or case.name.startswith(o.rstrip("/") + "/") for o in only):
            cases.append(case)
    return cases
