"""Times as text ("1 d", "30 min", "2.5h") or as a number of seconds, and the format of numbers for the output."""

import math
import re

SECOND = 1.0
MINUTE = 60.0
HOUR = 3600.0
DAY = 86400.0
YEAR = 365.25 * DAY

_UNITS = {
    "ns": 1e-9, "us": 1e-6, "ms": 1e-3,
    "s": SECOND, "sec": SECOND, "secs": SECOND, "second": SECOND, "seconds": SECOND,
    "m": MINUTE, "min": MINUTE, "mins": MINUTE, "minute": MINUTE, "minutes": MINUTE,
    "h": HOUR, "hr": HOUR, "hrs": HOUR, "hour": HOUR, "hours": HOUR,
    "d": DAY, "day": DAY, "days": DAY,
    "w": 7 * DAY, "wk": 7 * DAY, "week": 7 * DAY, "weeks": 7 * DAY,
    "y": YEAR, "yr": YEAR, "yrs": YEAR, "year": YEAR, "years": YEAR, "a": YEAR,
}
_TIME = re.compile(r"^\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*([A-Za-z]*)\s*$")


def parse_time(value, what: str = "time") -> float:
    """Seconds from a number (seconds) or a text with a unit, for example "1 d", "30min", "1e3 s", "10 y"."""
    if isinstance(value, bool):
        raise ValueError(f"{what}: give a number of seconds or a text like '1 d'")
    if isinstance(value, (int, float)):
        seconds = float(value)
    else:
        match = _TIME.match(str(value))
        if not match:
            raise ValueError(f"{what}: cannot read {value!r}, give for example '1 d', '30 min' or 3600")
        number, unit = match.groups()
        unit = unit.lower()
        if unit and unit not in _UNITS:
            raise ValueError(f"{what}: unknown unit {unit!r} in {value!r}, use s, min, h, d, w or y")
        seconds = float(number) * _UNITS.get(unit, SECOND)
    if not math.isfinite(seconds) or seconds < 0:
        raise ValueError(f"{what}: the time must be 0 or more")
    return seconds


def format_time(seconds: float) -> str:
    """A short text for a time: '1 d', '5 min', '3.2 y'."""
    if seconds == 0:
        return "0 s"
    for unit, size in (("y", YEAR), ("d", DAY), ("h", HOUR), ("min", MINUTE), ("s", SECOND)):
        if seconds >= size or unit == "s":
            value = seconds / size
            return f"{value:.4g} {unit}"
    return f"{seconds:.4g} s"


def sig(value: float, digits: int = 4) -> float:
    """The value with a few significant digits, so the output to the model is short."""
    if value == 0 or not math.isfinite(value):
        return value
    return float(f"{value:.{digits}g}")
