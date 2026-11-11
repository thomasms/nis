"""Download and unpack the ENDF files of a library (for now: the TENDL neutron files)."""

import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

# the neutron files of a TENDL release. The site: https://tendl.imperial.ac.uk (it was at
# https://tendl.web.psi.ch before). If the site moves again, give the address with --url.
TENDL_URL = "https://tendl.imperial.ac.uk/tendl_{year}/tar_files/TENDL-n.tgz"
LATEST_TENDL = "2025"


def download(url: str, target: Path) -> Path:
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_file():
        print(f"{target} is there, not downloaded again", file=sys.stderr)
        return target
    partial = target.with_suffix(target.suffix + ".part")
    print(f"download {url}", file=sys.stderr, flush=True)
    with urllib.request.urlopen(url) as response, open(partial, "wb") as out:
        shutil.copyfileobj(response, out, length=1 << 20)
    partial.replace(target)
    return target


def unpack(archive: Path, folder: Path) -> Path:
    """Unpack a tar archive safely (no paths out of the folder, no links or devices)."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tar:
        if hasattr(tarfile, "data_filter"):
            tar.extractall(folder, filter="data")
        else:
            root = folder.resolve()
            for member in tar.getmembers():
                path = (folder / member.name).resolve()
                if not (member.isfile() or member.isdir()) or root not in (path, *path.parents):
                    raise ValueError(f"unsafe member in {archive}: {member.name}")
            tar.extractall(folder)
    return folder


def fetch_tendl(year: str, folder: Path, url: str = None) -> Path:
    """Download TENDL-<year> neutron files into <folder>/TENDL-<year>-n/ and return that folder."""
    folder = Path(folder)
    url = url or TENDL_URL.format(year=year)
    archive = download(url, folder / f"TENDL-{year}-n{''.join(Path(url).suffixes) or '.tgz'}")
    return unpack(archive, folder / f"TENDL-{year}-n")
