"""
The nuclear data libraries that NIS can download, and how to make the NIS data files from them.

Each library gives some of these files in data/ (the name of the file is the key of the library):

    <key>.dec    the decay data, all the nuclides in one ENDF tape
    <key>.nfy    the neutron induced fission yields
    <key>.sfy    the spontaneous fission yields
    <key>.reac   the reaction data: the neutron files, processed with NJOY and PREPRO (see process.py)

The files come from the NEA data bank (data.oecd-nea.org) and the NNDC (www.nndc.bnl.gov). The MD5 of each
download is checked. Some files need small fixes before ENDFtk can read them (see fix_tape).
"""

import hashlib
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path

from . import endf, fetch

DECAY_SUBLIBRARY = 4
SF_YIELD_SUBLIBRARY = 5
NF_YIELD_SUBLIBRARY = 11

_NEA = "https://data.oecd-nea.org/api/records/{record}/files/{name}/content"
_NNDC = "https://www.nndc.bnl.gov/{path}"


@dataclass(frozen=True)
class Source:
    """A file to download: its address, its name on the disk and its MD5."""
    url: str
    name: str
    md5: str


def nea(record: str, name: str, md5: str) -> Source:
    return Source(_NEA.format(record=record, name=name), name, md5)


def nndc(path: str, md5: str) -> Source:
    return Source(_NNDC.format(path=path), path.rpartition("/")[2], md5)


@dataclass(frozen=True)
class Library:
    key: str
    title: str
    licence: str
    decay: Source = None
    nfy: Source = None
    sfy: Source = None
    neutron: Source = None

    def parts(self) -> list:
        return [name for name in ("decay", "nfy", "sfy", "neutron") if getattr(self, name)]


LIBRARIES = {lib.key: lib for lib in [
    Library("jeff311", "JEFF-3.1.1", "CC-BY-4.0",
            decay=nea("4xq18-6xa42", "JEFF311-Evaluated-RDD.zip", "fd55945f064be777f407a2452c39665c"),
            nfy=nea("yv03w-dxz31", "JEFF311NFY.ASC", "e0b8eb05e3a8c08a4ca3bda4d39fc262"),
            sfy=nea("yv03w-dxz31", "JEFF311SFY.ASC", "03d496b0c6a52d5290ffb07cf7c47a53")),
    Library("jeff4", "JEFF-4.0", "CC-BY-4.0",
            decay=nea("tw0c6-t1386", "Radioactive_Decay_Data_JEFF-40.txt", "1dbde55c142501b2359c49d5333d6e70"),
            nfy=nea("trtwt-k2828", "nf_Fission_Yields_JEFF-40.txt", "6eee8b1f4728068894577bb4cd8612af"),
            sfy=nea("tccek-e9s26", "sf_Fission_Yields_JEFF-40.txt", "1dd9436f80cd94b1fdbce2a78f0c1cd0"),
            neutron=nea("e9ajn-a3p20", "JEFF40-Evaluations-Neutron-593.zip", "51d00ee7bf1491d428f9b30a9782e41d")),
    Library("endfbvii1", "ENDF/B-VII.1", "public",
            decay=nndc("endf-b7.1/zips/ENDF-B-VII.1-decay.zip", "ad346888d55427c3f25d39f7fed3e659")),
    Library("endfbviii0", "ENDF/B-VIII.0", "public",
            decay=nndc("endf-b8.0/zips/ENDF-B-VIII.0_decay.zip", "aa80cd0a880d9d7e0905940b868370c3")),
]}


def get(key: str) -> Library:
    if key not in LIBRARIES:
        raise ValueError(f"no library '{key}', the libraries are: {', '.join(LIBRARIES)}")
    return LIBRARIES[key]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download(source: Source, folder: Path) -> Path:
    """Download the file to folder (if it is not there), and check its MD5."""
    path = fetch.download(source.url, Path(folder) / source.name)
    found = md5(path)
    if found != source.md5:
        bad = path.with_suffix(path.suffix + ".bad")
        path.replace(bad)
        raise ValueError(f"{source.name}: the MD5 is {found}, not {source.md5}. The file is now {bad}. "
                         f"If the file at {source.url} changed, the MD5 in nisdata/libraries.py must change too.")
    return path


def fix_tape(text: str, title: str) -> str:
    """
    Small fixes so that ENDFtk can read the tape. They do not change the data that NIS uses.
    - No first line (TPID): add one. The JEFF-4.0 fission yield files start with the first material.
    - A discrete spectrum in MF8/MT457 with NT=8 (the total internal conversion coefficient but not the K and L
      shell ones): ENDFtk reads only NT=6 or 12, so add RICK, dRICK, RICL and dRICL as 0 (NT=12).
      JEFF-4.0 decay has one (MAT 128).
    """
    lines = text.splitlines(keepends=True)
    if lines and endf.mat_mf_mt(lines[0])[1:] == (1, 451) and endf.mat_mf_mt(lines[0])[0] > 0:
        lines.insert(0, endf.tpid_line(title))
    for i, line in enumerate(lines):
        # the LIST head of a discrete spectrum: [ER, dER, 0, 0, NT, 0]
        if (line[70:75] == " 8457" and line[44:55].strip() == "8" and line[22:33].strip() == "0"
                and line[33:44].strip() == "0" and line[55:66].strip() == "0"):
            second = lines[i + 2] if i + 2 < len(lines) else ""
            if not second or second[22:66].strip():
                raise ValueError(f"line {i + 3}: a list of 8 items must have 2 items on its second line")
            lines[i] = line[:44] + f"{12:11d}" + line[55:]
            lines[i + 2] = second[:22] + " 0.000000+0" * 4 + second[66:]
    return "".join(lines)


def tape_from(path: Path, sublibrary: int, title: str) -> str:
    """
    The text of one ENDF tape from a download: a zip of one file for each nuclide (only the files of the
    sublibrary are used, in the order of ZA then LISO), or one tape.
    """
    if zipfile.is_zipfile(path):
        texts = []
        with zipfile.ZipFile(path) as z:
            for member in z.infolist():
                if member.is_dir():
                    continue
                text = z.read(member).decode("ascii", errors="replace")
                head = text.splitlines(keepends=True)[:4]
                try:
                    material = endf.material_from_lines(head, member.filename)
                except ValueError:
                    continue
                if material.nsub == sublibrary:
                    texts.append(text)
        if not texts:
            raise ValueError(f"no ENDF files of sublibrary {sublibrary} in {path}")
        return endf.combine_texts([fix_tape(t, title) for t in texts], title)
    return fix_tape(Path(path).read_text(errors="replace"), title)


def _write(text: str, target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    partial.write_text(text)
    partial.replace(target)
    count = len({endf.mat_mf_mt(line)[0] for line in text.splitlines()[1:] if line[70:75] == " 1451"})
    print(f"{target}: {count} materials", file=sys.stderr)
    return target


_PARTS = {
    "decay": ("dec", DECAY_SUBLIBRARY, "Radioactive Decay Data"),
    "nfy": ("nfy", NF_YIELD_SUBLIBRARY, "Neutron Induced Fission Product Yields"),
    "sfy": ("sfy", SF_YIELD_SUBLIBRARY, "Spontaneous Fission Product Yields"),
}


def is_there(path: Path) -> bool:
    """True if the file is there and is not a git lfs pointer file (a small text file, from an old checkout)."""
    path = Path(path)
    if not path.is_file():
        return False
    if path.stat().st_size < 1000:
        with open(path, "rb") as f:
            return not f.read(40).startswith(b"version https://git-lfs")
    return True


def make(lib: Library, part: str, raw: Path, data: Path, force: bool = False) -> Path:
    """Make data/<key>.<dec|nfy|sfy> of a library. A file that is there is not made again, unless force."""
    suffix, sublibrary, what = _PARTS[part]
    source = getattr(lib, part)
    if source is None:
        raise ValueError(f"{lib.title} has no {what.lower()} here")
    target = Path(data) / f"{lib.key}.{suffix}"
    if is_there(target) and not force:
        print(f"{target} is there, not made again (--force makes it again)", file=sys.stderr)
        return target
    path = download(source, Path(raw) / lib.key)
    return _write(tape_from(path, sublibrary, f"{lib.title} {what}"), target)


def fetch_neutron(lib: Library, raw: Path) -> Path:
    """Download and unpack the neutron files of a library. Returns the folder of the files."""
    if lib.neutron is None:
        raise ValueError(f"{lib.title} has no neutron files here")
    folder = Path(raw) / f"{lib.key}-n"
    archive = download(lib.neutron, Path(raw) / lib.key)
    if not (folder.is_dir() and any(folder.iterdir())):
        print(f"unpack {archive} to {folder}", file=sys.stderr)
        fetch.unpack_zip(archive, folder) if zipfile.is_zipfile(archive) else fetch.unpack_archive(archive, folder)
    return folder
