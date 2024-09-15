"""
Transforming between various types of celestial coordinates.
"""

import enum
import functools
from math import acos, asin, atan2, cos, degrees, radians, sin, tan

from .mathutils import PI2, reduce_rad

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


def _results_to_deg(func):  # type: ignore
    @functools.wraps(func)
    def wrapper(*args, **kwargs):  # type: ignore
        a, b = func(*args, **kwargs)
        return degrees(a), degrees(b)

    return wrapper


class ConversionType(enum.IntEnum):
    """Direction of conversion between the Ecliptic and the Equator."""

    """equatorial -> ecliptical"""
    EQU_TO_ECL = 1
    """ecliptical -> equatorial"""
    ECL_TO_EQU = -1


def _equecl(x: float, y: float, e: float, k: ConversionType) -> tuple[float, float]:
    """Base convertion routine.

    Args:
        x: longitude or right ascension in radians
        y: latitude or declination in radians
        k: direction of the convertion

    Returns:
        the pair of result coordinates in radians.
    """
    sinE = sin(e)
    cosE = cos(e)
    sinX = sin(x)
    a = atan2(sinX * cosE + k.value * (tan(y) * sinE), cos(x))
    b = asin(sin(y) * cosE - k.value * (cos(y) * sinE * sinX))
    return (reduce_rad(a), b)


def _equhor(x: float, y: float, phi: float) -> tuple[float, float]:
    """Converts between azimuth/altitude and hour-angle/declination.

    The equations are symmetrical in the two pairs of coordinates so that
    exactly the same code may be used to convert in either direction, there
    is no need to specify direction with a swich (see Dufett-Smith, page 35).

    Returns The pair of result coordinates.
    All angular values are in radians.
    """
    sx = sin(x)
    sy = sin(y)
    sphi = sin(phi)
    cx = cos(x)
    cy = cos(y)
    cphi = cos(phi)
    sq = (sy * sphi) + (cy * cphi * cx)
    q = asin(sq)
    cp = (sy - (sphi * sq)) / (cphi * cos(q))
    p = acos(cp)
    if sx > 0:
        p = PI2 - p

    return (p, q)


@_results_to_deg  # type: ignore
def ecl2equ(lmbda: float, beta: float, eps: float) -> tuple[float, float]:
    """Converts ecliptical to equatorial coordinates.

    Args:
        lmbda: longitude, degrees
        beta: latitude, degrees
        eps: obliquity of the ecliptic, degrees

    Returns:
        the pair of equatorial coordinates, (alpha, delta), in degrees.
    """
    return _equecl(
        radians(lmbda), radians(beta), radians(eps), ConversionType.ECL_TO_EQU
    )


@_results_to_deg  # type: ignore
def equ2ecl(alpha: float, delta: float, eps: float) -> tuple[float, float]:
    """Converts equatorial to ecliptical coordinates.

    Args:
        alpha: right ascension, degrees
        delta: declination, degrees
        eps: obliquity of the ecliptic, degrees

    Returns:
        the pair of ecliptic coordinates, (lambda, beta), in degrees.
    """
    return _equecl(
        radians(alpha), radians(delta), radians(eps), ConversionType.EQU_TO_ECL
    )


@_results_to_deg  # type: ignore
def equ2hor(h: float, delta: float, phi: float) -> tuple[float, float]:
    """Converts equatorial to horizontal coordinates.

    Args:
        h: local hour angle, in degrees, measured westwards from the South.
           h = LST - RA (RA = Right Ascension)
        delta: declination, in degrees
        phi: latitude, in degrees, positive in the Nothern hemisphere, negative in the Southern.

    Returns:
        * azimuth, in degrees, measured westward from the South
        * altitude, in degrees, positive above the horizon
    """

    return _equhor(radians(h), radians(delta), radians(phi))


@_results_to_deg  # type: ignore
def hor2equ(az: float, alt: float, phi: float) -> tuple[float, float]:
    """Converts horizontal to equatorial coordinates.

    Args:
        az: azimuth, in degrees, measured westwards from the South.

            `h = LST - RA` (RA = Right Ascension)

        alt: altitude, in degrees, positive above the horizon
        phi: latitude, in degrees, positive in the Nothern hemisphere, negative in the Southern.

    Returns:
        * hour angle, in degrees
        * declination*, in degrees
    """
    return _equhor(radians(az), radians(alt), radians(phi))
