"""
EXFOR 14 MeV cross sections: from the EXFOR tables of the IAEA to the standard format.

The measured cross sections of the reactions of the nuclides of the FNS foils (all the nuclides in the sample of an
FNS case, benchmarks/fns) with 14 MeV (D-T) neutrons, from 13.5 to 15 MeV. NIS gets each one from the reaction matrix
of a narrow group at the energy of the value, so the benchmark tests the reaction data of the reactions that make
the FNS decay heat, one reaction at a time.

The original files (<source>) are a copy of github.com/IAEA-NDS/exfortables_py, the EXFOR tables of the IAEA: one
file for each EXFOR data set (a subentry and its reaction), with the incident energy, the cross section and its
uncertainty. Only the neutron folders of the FNS nuclides are needed:

    git clone --depth 1 --filter=blob:none --sparse https://github.com/IAEA-NDS/exfortables_py.git
    git -C exfortables_py sparse-checkout set n/Fe-56 n/Al-27 ...      (or all of it)

The data is used as EXFOR gives it, nothing is changed or added. A data set is used if:

  1. its reaction code is (<target>(N,<process>)<product>,,SIG): a cross section at an incident energy. Spectrum
     averages (SPA, FIS, MXW, AV), derived (DERIV), evaluated (EVAL) and recommended (RECOM) values, ratios, partial
     cross sections (PAR) and the sum of more than one state (M+, G+M1, ...) are not.
  2. the state of the product is clear: no state, or -T (all the states); -G (the ground state); -M or -M1 (the
     first isomer, m in NIS); -M2 (the second isomer, n). Other codes, such as a level -L, are not used.
  3. NIS tracks the product in that state (it is in the decay data), and it is not the target itself (elastic and
     inelastic scattering).
  4. the process is X (the production of the product by all the reactions of the target), or the product is
     radioactive (it is measured by its activity, so the value is again the production by all the reactions). NIS
     gives this production. A stable product of one reaction, e.g. Ti50 of V51(n,d) from the deuterons, is not used:
     other reactions (V51(n,np)) make it too.

and only its values from 13.5 to 15 MeV with an uncertainty (more than 0) and a value more than 0 are used.

Each case is one reaction: all the data sets of the product from the target, one measurement for each (its source is
the EXFOR subentry, the process and the first author). The case name is <target>/<product>, e.g. Fe56/Mn56 (all
states), Nb93/Nb92m (the isomer) or Se82/Se81g (the ground state).
"""

import re
from pathlib import Path

from ..format import BENCHMARKS_DIR, Case, Measurement, load_benchmark
from .irdff2 import ELEMENTS, product_token

NAME = "exfor14"
# the copy of the EXFOR tables that the case files are made from
EXFOR_TABLES = "github.com/IAEA-NDS/exfortables_py, commit f9c2d73 (2026-09-11)"
METADATA = {
    "name": NAME,
    "title": "EXFOR 14 MeV cross sections",
    "category": "fusion",
    "description": "Measured cross sections of the reactions of the nuclides of the FNS foils with 14 MeV (D-T) "
                   "neutrons, from 13.5 to 15 MeV, from EXFOR. Each case is one reaction (a product from a target), "
                   "with all its EXFOR data sets. It tests the reaction data (and its isomer data) of the reactions "
                   "that make the FNS decay heat, one reaction at a time.",
    "references": [
        "N. Otuka et al., Towards a More Complete and Accurate Experimental Nuclear Reaction Data Library (EXFOR): "
        "International Collaboration Between Nuclear Reaction Data Centres (NRDC), Nuclear Data Sheets 120 (2014) "
        "272-276",
        "The reference of each data set is in its EXFOR entry (nds.iaea.org/EXFOR/<entry>), see the source of each "
        "measurement",
    ],
    "source": f"The EXFOR tables of the IAEA ({EXFOR_TABLES}), converted with scripts/bench.py convert exfor14 "
              f"(the rules to select the data are in scripts/nisbench/convert/exfor14.py). The products are checked "
              f"with the JEFF-3.1.1 decay data (data/jeff311.dec).",
}
SOURCE_HELP = "a copy of github.com/IAEA-NDS/exfortables_py (with the folders n/<target>, e.g. n/Fe-56)"
# the converter needs decay data, to know which products NIS tracks and which are stable
NEEDS_DECAY_DATA = True

ENERGY_MIN_MEV, ENERGY_MAX_MEV = 13.5, 15.0

# EXFOR state code of the product -> (state in NIS, all states)
STATES = {"": (0, True), "-T": (0, True), "-G": (0, False), "-M": (1, False), "-M1": (1, False), "-M2": (2, False)}

_CODE = re.compile(r"^\((\d+)-([A-Z]+)-(\d+)\(N,([^)]+)\)(\d+)-([A-Z]+)-(\d+)((?:-[A-Z0-9+]+)?),,SIG\)$")
_HEADER = re.compile(r"^#\s+([^:]+):\s*(.*?)\s*$")
_NUCLIDE = re.compile(r"^([A-Z][a-z]?)(\d+)$")


def fns_targets() -> list:
    """The nuclides of the samples of the FNS cases (ground states only)."""
    nuclides = set()
    for case in load_benchmark(BENCHMARKS_DIR / "fns").cases:
        nuclides |= {n for n in case.atoms if _NUCLIDE.match(n)}
    return sorted(nuclides)


def exfor_folder(nuclide: str) -> str:
    """Fe56 -> Fe-56, the name of its folder in the EXFOR tables."""
    element, mass = _NUCLIDE.match(nuclide).groups()
    return f"{element}-{mass}"


def zai(z: int, a: int, state: int = 0) -> int:
    return (z * 1000 + a) * 10 + state


def nis_name(z: int, a: int, state: int = 0) -> str:
    return f"{ELEMENTS[z]}{a}{'' if state == 0 else 'mn'[state - 1]}"


def read_table(path: Path) -> tuple:
    """(header, rows) of an EXFOR table: header is {key: text}, rows the (E MeV, dE MeV, value b, uncertainty b) of
    each line that has 4 numbers."""
    header, rows = {}, []
    for line in Path(path).read_text(errors="replace").splitlines():
        if line.startswith("#"):
            match = _HEADER.match(line)
            if match:
                header.setdefault(match.group(1), match.group(2))
            continue
        try:
            values = [float(v) for v in line.split()[:4]]
        except ValueError:
            continue
        if len(values) == 4:
            rows.append(tuple(values))
    return header, rows


def select(header: dict, rows: list, decay: dict) -> tuple:
    """(target, product, all_states, process, points) of a data set, or None if it is not used (see the rules at the
    top). decay: ZAI -> True if stable, for each nuclide of the decay data. points: (E eV, value mb, uncertainty mb)
    from 13.5 to 15 MeV, in the order of the file."""
    match = _CODE.match(header.get("Reaction code", ""))
    if not match:
        return None
    z, _, a, process, pz, _, pa, code = match.groups()
    if code not in STATES:
        return None
    z, a, pz, pa = int(z), int(a), int(pz), int(pa)
    state, all_states = STATES[code]
    product = zai(pz, pa, state)
    if product not in decay or (pz, pa) == (z, a) and state == 0:
        return None
    if process != "X" and decay[product]:
        return None
    points = [(e * 1e6, v * 1e3, dv * 1e3) for e, _, v, dv in rows
              if ENERGY_MIN_MEV <= e <= ENERGY_MAX_MEV and v > 0 and dv > 0]
    if not points:
        return None
    return nis_name(z, a), nis_name(pz, pa, state), all_states, process, points


def source_text(header: dict, process: str) -> str:
    """EXFOR 41240-012-0, (N,P), A.A.Filatenkov 1999"""
    entry = header.get("Entry id", "?").split()[0]
    return f"EXFOR {entry}, (N,{process}), {header.get('First author', '?')} {header.get('Year', '?')}"


def read_decay_data(path: Path) -> dict:
    """ZAI -> True if stable, for each nuclide of the decay data."""
    import pynis as nis
    data = nis.core.decay_data_from_file(str(path))
    return {int(k): bool(data[k].is_stable) for k in data.keys()}


def make_case(target: str, product: str, all_states: bool, data_sets: list) -> Case:
    """data_sets: (source, points) of each EXFOR data set."""
    token = product_token(product, all_states)
    measurements = []
    for source, points in sorted(data_sets):
        points = sorted(points, key=lambda p: p[0])
        measurements.append(Measurement("cross_section", "mb", [], [p[1] for p in points], [p[2] for p in points],
                                        target=target, product=product, all_states=all_states,
                                        energies_ev=[p[0] for p in points], source=source))
    return Case(benchmark=NAME, name=f"{target}/{token}", title=f"{target} to {token}, 13.5 to 15 MeV (EXFOR)",
                mass_kg=None, atoms={target: 1.0}, steps=[], measurements=measurements)


def convert(source: Path, only: list = None, decay_data: Path = None) -> list:
    """only: targets (Fe56) or case names (Fe56/Mn56), all if empty. decay_data: the NIS decay data file."""
    source = Path(source)
    decay = read_decay_data(decay_data or BENCHMARKS_DIR.parent / "data" / "jeff311.dec")
    cases = []
    for target in fns_targets():
        if only and not any(o == target or o.startswith(target + "/") for o in only):
            continue
        reactions = {}
        for path in sorted((source / "n" / exfor_folder(target)).glob("*/xs/*.txt")):
            header, rows = read_table(path)
            selected = select(header, rows, decay)
            if selected is None or selected[0] != target:
                continue
            _, product, all_states, process, points = selected
            reactions.setdefault((product, all_states), []).append((source_text(header, process), points))
        for (product, all_states), data_sets in sorted(reactions.items()):
            case = make_case(target, product, all_states, data_sets)
            if not only or any(case.name == o or o == target for o in only):
                cases.append(case)
    return cases
