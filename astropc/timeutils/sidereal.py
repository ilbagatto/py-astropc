"""
 Converts a given civil time into the local sidereal time and vice-versa.

 Sidereal and Civil time

 Sidereal time is reckoned by the daily transit of a fixed point in space
 (fixed with respect to the distant stars), 24 hours of sidereal time elapsing
 between an successive transits. The sidereal day is thus shorter than the
 solar day by nearely 4 minutes, and although the solar and sidereal time
 agree once a year, the difference between them grows systematically as the
 months pass in the sense that sidereal time runs faster than solar time.
 *Sidereal time* (ST) is used extensively by astronomers since it is the time
 kept by the star.

Caveats

 Times may be converted quite easily from UT to Greenwich mean sidereal time
 (SG) since there is a small range of sidereal times which occurs twice on
 the same calendar date. The ambiguous period is from 23h 56m 04s UT to
 Oh 03m 56s UT, i.e. about 4 minutes either side of midnight. The routine
 given here correctly converts SG to UT in the period before midnight, but
 not in the period after midnight when the ambiguity must be resolved by other
 means.

 -- Peter Duffett-Smith, "Astronomy with your PC"
"""

from astropc.mathutils import to_range

from .julian import cal_date, djd_midnight, jul_day

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"

_SID_RATE = 0.9972695677
_AMBIG_DELTA = 6.552e-2


def _tnaught(djd: float) -> float:
    (year, _, _) = cal_date(djd)
    dj0 = jul_day(year, 1, 0.0)
    t = dj0 / 36525
    x = 6.57098e-2 * (djd - dj0) - (
        24
        - (6.6460656 + (5.1262e-2 + (t * 2.581e-5)) * t)
        - (2400 * (t - ((year - 1900) / 100)))
    )
    return x


def djd_to_sidereal(djd: float, lng: float = 0.0) -> float:
    """Converts civil to Local Sidereal time.

    Args:
        djd: number of Julian days elapsed since 1900, Jan 0.5.
        lng: optional Geographic longitude, negative eastwards, 0.0 by default

    Returns:
        local sidereal time
    """
    djm = djd_midnight(djd)
    utc = (djd - djm) * 24
    t0 = _tnaught(djm)
    gst = (1.0 / _SID_RATE) * utc + t0
    lst = gst - lng / 15
    return to_range(lst, 24.0)


def sidereal_to_utc(lst: float, djd: float, lng: float = 0.0) -> tuple[float, bool]:
    """Converts Local Sidereal time, lst, to civil time.

    Args:
        djd: a number of Julian days elapsed since 1900, Jan 0.5.
        lng: optional Geographic longitude, negative eastwards, 0.0 by default

    Returns:
        * Universal Time
        * flag, which is True if the result is ambiguous (see Caveats)
    """
    djm = djd_midnight(djd)
    t0 = to_range(_tnaught(djm), 24.0)
    gst = lst + lng / 15
    utc = to_range(gst - t0, 24.0) * _SID_RATE
    return utc, utc < _AMBIG_DELTA
