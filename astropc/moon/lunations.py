"""Lunations
"""

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"

from dataclasses import dataclass
from enum import Enum
from math import cos, radians, sin

from astropc.mathutils import reduce_deg
from astropc.timeutils import day_of_year, is_leapyear


@dataclass
class QuarterDataMixin:
    title: str
    coeff: float


class Quarter(QuarterDataMixin, Enum):
    """Lunar Quarter."""

    NEW_MOON = "New Moon", 0.0
    FIRST_QUARTER = "First Quarter", 0.25
    FULL_MOON = "Full Moon", 0.5
    LAST_QUARTER = "Last Quarter", 0.75


def _calculate_delta(
    quarter: Quarter, t: float, ms: float, mm: float, f: float
) -> float:
    tms = ms + ms
    tmm = mm + mm
    tf = f + f

    if quarter in (Quarter.NEW_MOON, Quarter.FULL_MOON):
        delta = (
            (1.734e-1 - 3.93e-4 * t) * sin(ms)
            + 2.1e-3 * sin(tms)
            - 4.068e-1 * sin(mm)
            + 1.61e-2 * sin(tmm)
            - 4e-4 * sin(mm + tmm)
            + 1.04e-2 * sin(tf)
            - 5.1e-3 * sin(ms + mm)
            - 7.4e-3 * sin(ms - mm)
            + 4e-4 * sin(tf + ms)
            - 4e-4 * sin(tf - ms)
            - 6e-4 * sin(tf + mm)
            + 1e-3 * sin(tf - mm)
            + 5e-4 * sin(ms + tmm)
        )
    else:
        delta = (
            (0.1721 - 0.0004 * t) * sin(ms)
            + 0.0021 * sin(tms)
            - 0.6280 * sin(mm)
            + 0.0089 * sin(tmm)
            - 0.0004 * sin(tmm + mm)
            + 0.0079 * sin(tf)
            - 0.0119 * sin(ms + mm)
            - 0.0047 * sin(ms - mm)
            + 0.0003 * sin(tf + ms)
            - 0.0004 * sin(tf - ms)
            - 0.0006 * sin(tf + mm)
            + 0.0021 * sin(tf - mm)
            + 0.0003 * sin(ms + tmm)
            + 0.0004 * sin(ms - tmm)
            - 0.0003 * sin(tms + mm)
        )
        w = 0.0028 - 0.0004 * cos(ms) + 0.0003 * cos(ms)
        if quarter == Quarter.LAST_QUARTER:
            w = -w
        delta += w

    return delta


def find_closest_phase(quarter: Quarter, year: int, month: int, day: int) -> float:
    """Find DJD of a quarter, closest to the given date.

    Args:
        quarter: which quarter to search
        year: year
        month month
        day: Julian day (DJD) of the event

    """
    n = 366 if is_leapyear(year) else 365
    y = year + day_of_year(year, month, day) / n
    k = (
        round((y - 1900) * 12.3685) + quarter.coeff
    )  # TODO: find a better way to access fields
    t = k / 1236.85
    t2 = t * t
    t3 = t2 * t

    def assemble(a: float, b: float, c: float, d: float) -> float:
        return radians(reduce_deg(a + b * k + c * t2 + d * t3))

    c = radians(166.56 + (132.87 - 9.173e-3 * t) * t)
    j = (
        0.75933 + 29.53058868 * k + 0.0001178 * t2 - 1.55e-07 * t3 + 3.3e-4 * sin(c)
    )  # mean lunar phase

    ms = assemble(359.2242, 29.105356080, -0.0000333, -0.00000347)
    mm = assemble(306.0253, 385.81691806, 0.0107306, 0.00001236)
    f = assemble(21.2964, 390.67050646, -0.0016528, -0.00000239)

    return j + _calculate_delta(quarter, t, ms, mm, f)
