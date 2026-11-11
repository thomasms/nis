"""
The NIS input deck of a case. The deck is made from the case file, so do not edit it: edit the
case file and make the deck again (scripts/bench.py decks).
"""

from pathlib import Path

from .format import Case


def _numbers(values: list) -> str:
    return " ".join(f"{v:.9e}" for v in values)


def deck_text(case: Case) -> str:
    lines = [
        f"# {case.benchmark} benchmark, {case.name}: {case.title}",
        f"# made by scripts/bench.py from {case.path.name if case.path else 'the case file'}, do not edit",
    ]
    if case.mass_kg is not None:
        lines.append(f"# mass_kg {case.mass_kg:.9e}")
    if case.spectrum is not None:
        structure = f" ({case.group_structure})" if case.group_structure else ""
        lines.append(f"# spectrum: {len(case.spectrum[1])} group flux shape{structure}, "
                     "nis normalises it to the flux of each step")
    lines.append("# steps: dt in s, flux in n/m2/s (0 is cooling)")
    lines += [f"dt {dt:.9e} {flux:.9e}" if flux > 0 else f"dt {dt:.9e} 0" for dt, flux in case.steps]
    lines += ["", "nuclides"]
    lines += [f"    {name} {atoms:.9e}" for name, atoms in sorted(case.atoms.items())]
    if case.spectrum is not None:
        bounds, values = case.spectrum
        lines += ["", "spectrum", f"    {_numbers(bounds)}", f"    {_numbers(values)}"]
    return "\n".join(lines) + "\n"


def write_deck(case: Case, path: Path = None) -> Path:
    path = Path(path or case.deck_path)
    path.write_text(deck_text(case))
    return path


def deck_is_current(case: Case) -> bool:
    return case.deck_path.exists() and case.deck_path.read_text() == deck_text(case)
