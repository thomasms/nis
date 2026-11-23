"""
The production of each product of MT=5 (n,anything) in MF10, made from the MF6 data.

Many evaluations (TENDL, and the JEFF-4.0 files that come from TENDL) give the products of MT=5 only in MF6:
a cross section in MF3/MT5, and for each product a yield y(E) in MF6/MT5. PREPRO SIXPAK makes MF9 only for the
isomers of the other reactions, so without this step the MT=5 products are lost. For each product (also the light
particles and the photons, as in the UKAEA files) this writes the cross section sigma_5(E) * y(E) in MF10/MT5,
with LFS = the isomer flag (LIP) of the product. NIS reads the LFS of MT=5 as the isomer number.
"""

import bisect
from pathlib import Path

from .endf import endf_float, endf_int, mat_mf_mt

MT = 5
MF_XS, MF_PRODUCTS, MF_OUT = 3, 6, 10


def _lines_of(text_lines: list, mf: int, mt: int) -> list:
    return [line for line in text_lines if len(line) >= 75 and mat_mf_mt(line)[1:] == (mf, mt)]


class _Reader:
    """Reads the records of one section, line by line (the 6 fields of 11 characters)."""

    def __init__(self, lines: list):
        self.lines, self.i = lines, 0

    def cont(self) -> tuple:
        line = self.lines[self.i]
        self.i += 1
        return (endf_float(line[0:11]), endf_float(line[11:22]), endf_int(line[22:33]), endf_int(line[33:44]),
                endf_int(line[44:55]), endf_int(line[55:66]))

    def numbers(self, n: int) -> list:
        values = []
        while len(values) < n:
            line = self.lines[self.i]
            self.i += 1
            values += [endf_float(line[11 * k:11 * k + 11]) for k in range(6)]
        return values[:n]

    def skip(self, n: int) -> None:
        self.i += (n + 5) // 6

    def tab1(self) -> tuple:
        head = self.cont()
        nr, np_ = head[4], head[5]
        interp = self.numbers(2 * nr)
        xy = self.numbers(2 * np_)
        breaks = [int(interp[2 * k]) for k in range(nr)]
        laws = [int(interp[2 * k + 1]) for k in range(nr)]
        return head, xy[0::2], xy[1::2], breaks, laws

    def skip_tab2(self) -> tuple:
        head = self.cont()
        self.skip(2 * head[4])
        return head

    def skip_list(self) -> tuple:
        head = self.cont()
        self.skip(head[4])
        return head


def _skip_law(reader: _Reader, law: int) -> None:
    """The energy-angle data of one product of MF6 (it is not needed here)."""
    if law in (0, 3, 4):
        return
    if law in (1, 2, 5):
        for _ in range(reader.skip_tab2()[5]):
            reader.skip_list()
    elif law == 6:
        reader.cont()
    elif law == 7:
        for _ in range(reader.skip_tab2()[5]):
            for _ in range(reader.skip_tab2()[5]):
                reader.tab1()
    else:
        raise ValueError(f"MF6/MT5: LAW={law} is not known")


def _value(x: list, y: list, breaks: list, laws: list, e: float) -> float:
    """The TAB1 at e: 0 out of its range, histogram for INT=1, else linear."""
    if not x or e < x[0] or e > x[-1]:
        return 0.0
    j = bisect.bisect_right(x, e) - 1
    if j >= len(x) - 1:
        return y[-1]
    law = laws[bisect.bisect_left(breaks, j + 2)] if laws else 2
    if law == 1 or x[j + 1] == x[j]:
        return y[j]
    return y[j] + (y[j + 1] - y[j]) * (e - x[j]) / (x[j + 1] - x[j])


def _number(v: float) -> str:
    """An ENDF number in 11 characters, 1.234567+5."""
    if v == 0.0:
        return " 0.000000+0"
    mantissa, exponent = f"{v:.6e}".split("e")
    return f"{mantissa}{int(exponent):+d}".rjust(11)


def _records(mat: int, rows: list, seq_start: int) -> list:
    """Lines of 6 fields with MAT, MF10, MT5 and sequence numbers."""
    lines = []
    for k, row in enumerate(rows):
        lines.append(f"{row:66s}{mat:4d}{MF_OUT:2d}{MT:3d}{(seq_start + k) % 100000:5d}\n")
    return lines


def _fields(values: list) -> list:
    rows = []
    for k in range(0, len(values), 6):
        rows.append("".join(values[k:k + 6]))
    return rows


def residual_production(text_lines: list) -> list:
    """
    The lines of MF10/MT5 (the section and its SEND line) from the MF3/MT5 and MF6/MT5 of the material, or []
    if the material has no MF6/MT5 (or no MF3/MT5, or MF10/MT5 already).
    """
    if _lines_of(text_lines, MF_OUT, MT):
        return []
    xs_lines, product_lines = _lines_of(text_lines, MF_XS, MT), _lines_of(text_lines, MF_PRODUCTS, MT)
    if not xs_lines or not product_lines:
        return []
    mat = mat_mf_mt(xs_lines[0])[0]

    xs = _Reader(xs_lines)
    za, awr = xs.cont()[0:2]
    _, energies, sigma, xs_breaks, xs_laws = xs.tab1()

    reader = _Reader(product_lines)
    nk = reader.cont()[4]
    products = []
    for _ in range(nk):
        head, ye, yv, breaks, laws = reader.tab1()
        zap, lip, law = int(round(head[0])), head[2], head[3]
        _skip_law(reader, law)
        # all the energies of the cross section and of the yield, in the range of the yield
        grid = sorted(set(e for e in energies + ye if ye[0] <= e <= ye[-1]))
        values = [_value(energies, sigma, xs_breaks, xs_laws, e) * _value(ye, yv, breaks, laws, e) for e in grid]
        if any(v > 0.0 for v in values):
            products.append((zap, lip, grid, values))

    if not products:
        return []
    rows = ["".join([_number(za), _number(awr), f"{0:11d}", f"{0:11d}", f"{len(products):11d}", f"{0:11d}"])]
    for zap, lip, grid, values in products:
        rows.append("".join([_number(0.0), _number(0.0), f"{zap:11d}", f"{lip:11d}", f"{1:11d}", f"{len(grid):11d}"]))
        rows.append(f"{len(grid):11d}{2:11d}")
        pairs = []
        for e, v in zip(grid, values):
            pairs += [_number(e), _number(v)]
        rows += _fields(pairs)
    lines = _records(mat, rows, 1)
    lines.append(f"{'':66s}{mat:4d}{MF_OUT:2d}{0:3d}{99999:5d}\n")
    return lines


def add_residual_production(source: Path, output: Path, products_from: Path = None) -> int:
    """
    Writes source to output, with MF10/MT5 made from the MF3/MT5 and MF6/MT5 of products_from (default: source).
    MF10/MT5 goes before the other MF10 sections (or in a new MF10, before MF > 10). Returns the number of
    products, 0 if there is no MT5 to add (then output is a copy of source).
    """
    lines = Path(source).read_text().splitlines(keepends=True)
    data = Path(products_from or source).read_text().splitlines(keepends=True)
    section = residual_production(data)
    if section and _lines_of(lines, MF_OUT, MT):
        section = []
    if section:
        mat = mat_mf_mt(section[0])[0]
        keys = [mat_mf_mt(line) if len(line) >= 75 else (None, None, None) for line in lines]
        mf10 = [i for i, k in enumerate(keys) if k[0] == mat and k[1] == MF_OUT]
        if mf10:
            # before the first MF10 section with MT > 5, else after the last MF10 line (before the FEND)
            after = [i for i in mf10 if keys[i][2] > MT]
            at, new_file = (after[0] if after else mf10[-1] + 1), False
        else:
            # a new MF10 before the first file after it, else before the MEND of the material
            later = [i for i, k in enumerate(keys) if k[0] == mat and k[1] > MF_OUT]
            last = max(i for i, k in enumerate(keys) if k[0] == mat)
            at, new_file = (later[0] if later else last + 1), True
        lines[at:at] = section + ([f"{'':66s}{mat:4d}{0:2d}{0:3d}{0:5d}\n"] if new_file else [])
    Path(output).write_text("".join(lines))
    return int(section[0][44:55]) if section else 0
