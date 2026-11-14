"""
Fission decay heat experiments: from the IAEA CoNDERC fission archive to the standard format.

The archive (fission.zip from the CoNDERC web site, www-nds.iaea.org/conderc, fission section) has one CSV file
for each experiment in fission/expDataSystem/<nuclide>_<thermal|fast>/<Author>_<pulse|irradiation s>.csv:

    #Author:,Dickens
    #Irradiation:,Pulse
    #T(s),Beta(MeV/f/s),dBeta(MeV/f/s),Gamma(MeV/f/s),dGamma(MeV/f/s),Total(MeV/f/s),dTotal(MeV/f/s)
    2.7000E+00,6.4312E-01,5.3371E-02,...

The headers do not give the true units. From the reports and a check against NIS for each file:
    - a pulse: t * P(t), MeV per fission (unit MeV/fission)
    - an irradiation of T s at a constant fission rate: P(t) over the fission rate, MeV per fission/s
      (unit MeV/(fission/s)). Baumung (200 s) is t * P(t) over the fissions of the irradiation (MeV/fission).
A "-", an empty cell or a 0 is no value. Some files need more:
    - Jurney measured the gamma heat only; in the U235 file it is in the Total column.
    - Dickens Pu239 is three series of measurements one after the other; each is a case.
    - Fisher (Godiva, 1964) is left out: no one unit makes the beta and gamma values agree with each other.
    - Unat is left out: the fissions are of both U235 and U238.
    - The header of U235_thermal/Yarnell_20000.csv says Akiyama and Pulse; the file name is right.

Each case is a fission source (the independent yields of the nuclide at 0.0253 eV for thermal, 400 keV for fast
fission), so no reaction data is needed: one fission at the start for a pulse, or a fission rate in the
irradiation step. Then one cooling step to each measured time.
"""

import io
import zipfile
from pathlib import Path

from ..format import Case, Measurement

NAME = "fission_decay_heat"
METADATA = {
    "name": NAME,
    "title": "Fission decay heat",
    "category": "fission",
    "description": "Beta, gamma and total decay heat after the fission of U233, U235, U238, Pu239, Pu241 and "
                   "Th232 by thermal or fast neutrons: fission pulses, and irradiations from 10 s to 1e5 s, "
                   "measured from 0.25 s to 300 days after the fissions. The decay heat is per fission (pulse) or "
                   "per fission/s (irradiation).",
    "references": [
        "J.K. Dickens et al., Fission product energy release for times following thermal-neutron fission of U235 "
        "between 2 and 14000 seconds, ORNL/NUREG-14 (1979); J.K. Dickens, ORNL/TM-10094 (1986)",
        "M. Akiyama et al., J. Nucl. Sci. Technol. 24 (1982) 709-722 (gamma) and 803-816 (beta)",
        "A. Tobias, Derivation of decay heat benchmarks for U235 and Pu239 by a least squares fit to measured "
        "data, CEGB (1989)",
        "W.A. Schier, G.P. Couchell, DOE/ER/40723-4, University of Massachusetts Lowell (1997)",
        "K. Baumung, KfK 3262, Kernforschungszentrum Karlsruhe (1981)",
        "J.L. Yarnell, P.J. Bendt, LA-NUREG-6713 (1977) and LA-7452-MS (1978)",
        "M. Lott, G. Lhiaubet, F. Defreche, Journal of Nuclear Energy 27 (1973) 597-605",
        "E. Jurney, P. Bendt, T. England, Fission product gamma spectra, LA-7620-MS (1979)",
        "T.D. MacMahon, R. Wellum, H.W. Wilson, Journal of Nuclear Energy 24 (1970) 493-501",
        "J. Scobie, R.D. Scott, H.W. Wilson, Journal of Nuclear Energy 25 (1971) 1-10",
        "A. McNair et al., Journal of Nuclear Energy 23 (1969) 73-86 (U235) and 697-703 (Pu239)",
        "M. Murphy et al., AEEW-R 1212, UK Atomic Energy Authority (1979)",
        "M. Fleming, J.-Ch. Sublet, Validation of FISPACT-II decay heat and inventory predictions for fission "
        "events, UKAEA-R(18)003 (2018)",
    ],
    "source": "The IAEA CoNDERC fission archive (fission.zip, fission/expDataSystem), converted with "
              "scripts/bench.py convert fission_decay_heat. Each case is a fission source with the independent "
              "fission yields (thermal: 0.0253 eV, fast: 400 keV); run it with --fission-yields.",
}
SOURCE_HELP = "the CoNDERC fission archive (fission.zip) or the folder it was extracted to"

NEUTRON_ENERGY_EV = {"thermal": 0.0253, "fast": 4.0e5}
LEFT_OUT = {"Fisher": "the beta and gamma values do not agree with one unit", "Unat": "fissions of U235 and U238"}
_COLUMNS = (("beta_heat", 1), ("gamma_heat", 3), ("decay_heat", 5))
_TITLES = {"thermal": "thermal", "fast": "fast"}


def _files(source: Path) -> dict:
    """relative path (<nuclide>_<spectrum>/<file>.csv) -> text"""
    source = Path(source)
    out = {}
    if source.is_file():
        with zipfile.ZipFile(source) as z:
            for name in z.namelist():
                parts = name.split("/")
                if name.endswith(".csv") and len(parts) >= 3 and parts[-3] == "expDataSystem":
                    out["/".join(parts[-2:])] = z.read(name).decode("utf-8", "replace")
    else:
        root = next((p for p in source.rglob("expDataSystem") if p.is_dir()), None)
        if root is not None:
            for path in root.glob("*/*.csv"):
                out[f"{path.parent.name}/{path.name}"] = path.read_text(errors="replace")
    return out


def _number(cell: str):
    cell = cell.strip()
    if cell in ("", "-"):
        return None
    value = float(cell)
    return value if value > 0 else None


def _rows(text: str) -> list:
    return [[c for c in line.split(",")] for line in io.StringIO(text).read().splitlines()
            if line.strip() and not line.startswith("#")]


def _series(rows: list) -> list:
    """The rows in series where the time goes up (a file can have more than one series)."""
    series = [[rows[0]]]
    for row in rows[1:]:
        if float(row[0]) <= float(series[-1][-1][0]):
            series.append([])
        series[-1].append(row)
    return series


def _case(relative: str, rows: list, suffix: str = ""):
    folder, filename = relative.split("/")
    nuclide, spectrum = folder.split("_")
    author, irradiation = Path(filename).stem.split("_")
    pulse = irradiation.lower() == "pulse"
    irradiation_s = 0.0 if pulse else float(irradiation)

    times = [float(r[0]) for r in rows]
    values = {}
    for quantity, column in _COLUMNS:
        for i, row in enumerate(rows):
            if author == "Jurney" and quantity != "gamma_heat":
                continue
            # Jurney measured the gamma heat only (in the U235 file it is in a column one to the right)
            source_column = column
            if author == "Jurney":
                source_column = next((c for c in (3, 5, 6) if len(row) > c + 1 and _number(row[c]) is not None), None)
                if source_column is None:
                    continue
            if len(row) <= source_column + 1:
                continue
            value = _number(row[source_column])
            if value is None:
                continue
            uncertainty = _number(row[source_column + 1]) or 0.0
            values.setdefault(quantity, []).append((i, value, uncertainty))
    if not values:
        return None

    # pulse: one fission at the start; irradiation: 1 fission/s (Baumung: 1 fission in all)
    rate = 0.0 if pulse else (1.0 / irradiation_s if author == "Baumung" else 1.0)
    steps, rates = [], []
    if not pulse:
        steps.append((irradiation_s, 0.0))
        rates.append(rate)
    previous = 0.0
    for t in times:
        steps.append((t - previous, 0.0))
        rates.append(0.0)
        previous = t
    first = 0 if pulse else 1
    unit = "MeV/fission" if pulse or author == "Baumung" else "MeV/(fission/s)"
    measurements = [Measurement(quantity, unit, [first + i for i, _, _ in v], [x for _, x, _ in v],
                                [u for _, _, u in v]) for quantity, v in values.items()]

    what = "fission pulse" if pulse else f"{irradiation_s:g} s irradiation"
    return Case(benchmark=NAME, name=f"{folder}/{author}_{irradiation}{suffix}",
                title=f"{nuclide} {_TITLES[spectrum]} fission, {what}, {author}{suffix.replace('_', ' ')}",
                mass_kg=None, atoms={}, steps=steps, measurements=measurements,
                fission_source={"nuclide": nuclide, "neutron_energy_ev": NEUTRON_ENERGY_EV[spectrum],
                                "fissions_at_start": 1.0 if pulse else 0.0},
                fission_rates=None if pulse else rates)


def convert(source: Path, only: list = None) -> list:
    cases = []
    for relative, text in sorted(_files(source).items()):
        folder, filename = relative.split("/")
        nuclide = folder.split("_")[0]
        author = Path(filename).stem.split("_")[0]
        if author in LEFT_OUT or nuclide in LEFT_OUT:
            continue
        rows = _rows(text)
        if not rows:
            continue
        series = _series(rows)
        for n, part in enumerate(series):
            case = _case(relative, part, f"_series{n + 1}" if len(series) > 1 else "")
            if case is not None and (not only or any(case.name.startswith(o) for o in only)):
                cases.append(case)
    return cases
