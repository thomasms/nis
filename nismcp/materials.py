"""
A material as mass fractions of nuclides ("U235", "Co59", "Am242m") and natural elements ("Fe", "U"), and its
number of atoms of each nuclide.

The atomic masses and the natural abundances come from the periodictable package (NIST and IUPAC values).
"""

import re

AVOGADRO = 6.02214076e23

# where periodictable is wrong for a natural element, the IUPAC values (atom %). The key is the mass number, with
# "m" for a metastable state.
# - U: periodictable has no natural abundance (IUPAC 2013 values).
# - Ta: the natural Ta180 is the metastable state Ta180m (almost stable), not the ground state (8 h half-life).
# - Th: natural Th is all Th232 (IUPAC); the 0.02 % of Th230 in periodictable is from uranium ores only.
_ABUNDANCE_OVERRIDES = {
    "U": {234: 0.0054, 235: 0.7204, 238: 99.2742},
    "Ta": {"180m": 0.01201, 181: 99.98799},
    "Th": {232: 100.0},
}

_ELEMENT = re.compile(r"^([A-Za-z]{1,2})$")
_NUCLIDE = re.compile(r"^([A-Za-z]{1,2})-?(\d{1,3})(m\d?|n)?$")


def _element(symbol: str):
    import periodictable
    symbol = symbol[0].upper() + symbol[1:].lower()
    try:
        return periodictable.elements.symbol(symbol)
    except ValueError:
        raise ValueError(f"{symbol!r} is not an element") from None


def natural_abundances(symbol: str) -> dict:
    """{mass number: atom fraction} of the natural element. The mass number is an int, or a str such as "180m"
    for a metastable state."""
    element = _element(symbol)
    if element.symbol in _ABUNDANCE_OVERRIDES:
        percent = _ABUNDANCE_OVERRIDES[element.symbol]
    else:
        percent = {a: element[a].abundance for a in element.isotopes if element[a].abundance > 0}
    if not percent:
        raise ValueError(f"{element.symbol} has no natural isotopes, give the nuclides (for example "
                         f"{element.symbol}{max(element.isotopes)})")
    total = sum(percent.values())
    return {a: p / total for a, p in percent.items()}


def atomic_mass(symbol: str, a) -> float:
    """The atomic mass in g/mol (u). If the table has no value for the nuclide, the mass number. a is the mass
    number (an int, or a str such as "180m": a metastable state has the mass of the ground state here)."""
    a = int(str(a).rstrip("mn"))
    element = _element(symbol)
    if a in element.isotopes and element[a].mass:
        return float(element[a].mass)
    return float(a)


def canonical_nuclide(name: str, known: dict) -> str:
    """The NIS name of a nuclide ("u235" -> "U235", "Tc-99m" -> "Tc99m", "Am242m1" -> "Am242m").
    known: the NIS names, the nuclides of the decay data."""
    match = _NUCLIDE.match(name.strip())
    if not match:
        raise ValueError(f"cannot read the nuclide {name!r}, give for example U235, Co60, Tc99m")
    symbol, a, state = match.groups()
    symbol = _element(symbol).symbol
    state = (state or "").lower()
    if state in ("m1", "m"):
        state = "m"
    elif state == "m2":
        state = "n"
    nis_name = f"{symbol}{int(a)}{state}"
    if nis_name not in known:
        raise ValueError(f"the decay data has no nuclide {nis_name}")
    return nis_name


def is_element(name: str) -> bool:
    return bool(_ELEMENT.match(name.strip()))


def atoms_from_mass(mass_kg: float, composition: dict, known: dict) -> tuple:
    """The atoms of each nuclide, and notes for the user.

    composition: {name: mass fraction}. A name is a nuclide (U235) or a natural element (Fe). The fractions are
    normalised to a sum of 1, so percent or fractions both work."""
    if not composition:
        raise ValueError("give the material, for example {'U235': 1.0} or {'Fe': 0.7, 'Cr': 0.18, 'Ni': 0.12}")
    if not mass_kg > 0:
        raise ValueError("the mass must be more than 0 kg")
    fractions = {}
    for name, fraction in composition.items():
        fraction = float(fraction)
        if fraction < 0:
            raise ValueError(f"the fraction of {name} must be 0 or more")
        fractions[name] = fractions.get(name, 0.0) + fraction
    total = sum(fractions.values())
    if total <= 0:
        raise ValueError("the sum of the mass fractions must be more than 0")
    notes = []
    if abs(total - 1.0) > 1e-3 and abs(total - 100.0) > 0.1:
        notes.append(f"the mass fractions add to {total:g}, so they were normalised to 1")

    atoms = {}
    grams = mass_kg * 1000.0
    for name, fraction in fractions.items():
        if fraction == 0:
            continue
        element_grams = grams * fraction / total
        if is_element(name):
            element = _element(name)
            abundances = natural_abundances(name)
            # the mean atomic mass of the natural element, from the same abundances
            mean_mass = sum(f * atomic_mass(element.symbol, a) for a, f in abundances.items())
            moles = element_grams / mean_mass
            for a, f in abundances.items():
                nuclide = f"{element.symbol}{a}"
                if nuclide not in known:
                    notes.append(f"{nuclide} is not in the decay data, so it was left out")
                    continue
                atoms[nuclide] = atoms.get(nuclide, 0.0) + moles * f * AVOGADRO
        else:
            nuclide = canonical_nuclide(name, known)
            symbol, a = _NUCLIDE.match(nuclide).groups()[:2]
            moles = element_grams / atomic_mass(symbol, int(a))
            atoms[nuclide] = atoms.get(nuclide, 0.0) + moles * AVOGADRO
    return atoms, notes


def check_atoms(atoms: dict, known: dict) -> dict:
    """The atoms with the NIS names. atoms: {nuclide: number of atoms}."""
    if not atoms:
        raise ValueError("give at least one nuclide")
    result = {}
    for name, value in atoms.items():
        value = float(value)
        if value < 0:
            raise ValueError(f"the atoms of {name} must be 0 or more")
        nuclide = canonical_nuclide(name, known)
        result[nuclide] = result.get(nuclide, 0.0) + value
    return result
