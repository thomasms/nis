"""Download and unpack the ENDF files of a library (for now: the TENDL neutron files)."""

import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

# the neutron files of a TENDL release. The site: https://tendl.imperial.ac.uk (it was at
# https://tendl.web.psi.ch before). If the site moves again, give the address with --url.
TENDL_URL = "https://tendl.imperial.ac.uk/tendl_{year}/tar_files/TENDL-n.tgz"
LATEST_TENDL = "2025"


def _size(n: float) -> str:
    for unit in ("B", "kB", "MB", "GB"):
        if n < 1000 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.2f} {unit}"
        n /= 1000


def _duration(seconds: float) -> str:
    seconds = int(seconds)
    if seconds >= 3600:
        return f"{seconds // 3600}h{seconds % 3600 // 60:02d}m"
    return f"{seconds // 60}m{seconds % 60:02d}s"


class Progress:
    """The progress of a download on stderr: one line that changes in a terminal, else a line every 30 s."""

    def __init__(self, total: int = None, done: int = 0, stream=sys.stderr):
        self.total, self.done, self.stream = total, done, stream
        self.start_done, self.start = done, time.monotonic()
        self.terminal = hasattr(stream, "isatty") and stream.isatty()
        self.interval = 0.5 if self.terminal else 30.0
        self.shown = 0.0

    def text(self) -> str:
        elapsed = time.monotonic() - self.start
        speed = (self.done - self.start_done) / elapsed if elapsed > 0 else 0.0
        if self.total:
            fraction = min(1.0, self.done / self.total)
            width = 30
            bar = "#" * int(width * fraction) + "-" * (width - int(width * fraction))
            left = f", {_duration((self.total - self.done) / speed)} left" if speed > 0 else ""
            return (f"[{bar}] {fraction:6.1%}  {_size(self.done)} of {_size(self.total)}"
                    f"  {_size(speed)}/s{left}")
        return f"{_size(self.done)}  {_size(speed)}/s"

    def add(self, n: int) -> None:
        self.done += n
        now = time.monotonic()
        if now - self.shown >= self.interval:
            self.shown = now
            self.show()

    def show(self) -> None:
        if self.terminal:
            self.stream.write("\r" + self.text() + "\033[K")
        else:
            self.stream.write(self.text() + "\n")
        self.stream.flush()

    def end(self) -> None:
        self.show()
        if self.terminal:
            self.stream.write("\n")
            self.stream.flush()


def _copy(response, out, progress: Progress, chunk: int = 1 << 20) -> None:
    while True:
        data = response.read(chunk)
        if not data:
            break
        out.write(data)
        progress.add(len(data))


def download(url: str, target: Path, tries: int = 50) -> Path:
    """Download to <target>.part, then rename. A stopped download continues from where it stopped."""
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_file():
        print(f"{target} is there, not downloaded again", file=sys.stderr)
        return target
    partial = target.with_suffix(target.suffix + ".part")
    print(f"download {url}", file=sys.stderr, flush=True)
    for attempt in range(1, tries + 1):
        done = partial.stat().st_size if partial.is_file() else 0
        request = urllib.request.Request(url, headers={"Range": f"bytes={done}-"} if done else {})
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                # 206: the server sends the rest. 200: it sends all of it again
                mode = "ab" if done and response.status == 206 else "wb"
                total = response.headers.get("Content-Range", "").rpartition("/")[2] or \
                    response.headers.get("Content-Length")
                start = done if mode == "ab" else 0
                if start:
                    print(f"continue from {_size(start)}", file=sys.stderr, flush=True)
                progress = Progress(int(total) if total and total.isdigit() else None, start)
                try:
                    with open(partial, mode) as out:
                        _copy(response, out, progress)
                finally:
                    progress.end()
            if total and total.isdigit() and partial.stat().st_size != int(total):
                raise OSError(f"got {partial.stat().st_size} of {total} bytes")
            break
        except urllib.error.HTTPError as e:
            # 416: the part file is already complete
            if e.code == 416 and done:
                break
            raise
        except OSError as e:
            if attempt == tries:
                raise
            print(f"download stopped ({e}), continue ({attempt}/{tries})", file=sys.stderr, flush=True)
            time.sleep(min(60, 5 * attempt))
    partial.replace(target)
    return target


def unpack_archive(archive: Path, folder: Path) -> Path:
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


def fetch_tendl(year: str, folder: Path, url: str = None, unpack: bool = True) -> Path:
    """
    Download the TENDL-<year> neutron files to <folder>/TENDL-<year>-n.tgz. With unpack, unpack them into
    <folder>/TENDL-<year>-n/ and return that folder, else return the archive.
    """
    folder = Path(folder)
    url = url or TENDL_URL.format(year=year)
    archive = download(url, folder / f"TENDL-{year}-n{''.join(Path(url).suffixes) or '.tgz'}")
    return unpack_archive(archive, folder / f"TENDL-{year}-n") if unpack else archive
