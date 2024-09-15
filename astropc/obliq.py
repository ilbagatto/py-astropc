"""Calculates obliquity of the ecliptic.

Source: P.Duffett-Smith, "Astronomy with Your PC", 2 edition.
"""

from .timeutils import DAYS_PER_CENT

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


def calc_obliquity(djd: float, deps: float = 0.0) -> float:
    """Calculate obliquity of ecliptic.

    Args:
        djd (float): number of Julian days since 1900 Jan. 0.5.
        deps (float, optional): nutation in ecliptic obliquity in arc-degrees. Defaults to 0.0.

    Returns:
        float: obliquity of ecliptic in degrees.

    With the second argument, calculates *true obliquity*. Otherwise, *true obliquity*.
    """
    t = djd / DAYS_PER_CENT
    c = (((-0.00181 * t) + 0.0059) * t + 46.845) * t
    return 23.45229444 - (c / 3600) + deps
