"""Core mathematical routines for astronomic calculations.
"""

from collections import namedtuple
from functools import partial, reduce
from math import fabs, fmod, modf, pi
from operator import mul

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"

PI2 = pi * 2
"""360 degrees in radians"""


PI_HALF = pi / 2.0
"""`pi / 2`"""


def polynome(t: float, a: float, *args: float) -> float:
    """Calculates polynome:

        a1 + a2*t + a3*t*t + a4*t*t*t...

    Args:
        t: number of Julian centuries elapsed since 1900, Jan 0.5.
        a: the first, required coefficient.
        terms: optionally, any number of additional coefficients.

    Returns:
        result of the expression


    >>> polynome(10.0, 1.0, 2.0, 3.0)
    321.0
    """
    coeffs = (a,) + args
    tmul = partial(mul, t)
    x = reduce(lambda a, b: tmul(a) + b, reversed(coeffs))
    return x


def to_range(x: float, r: float) -> float:
    """Reduces x to given range

        0 >= x < r

    Args:
        x: number to normalize.
        r: range

    Returns:
        normalized number

    >>> to_range(-700, 360)
    20.0
    """
    a = fmod(x, r)
    return a + r if a < 0 else a


def reduce_deg(x: float) -> float:
    """Reduces a number to range 0 .. 360

    Args:
        x: arc-degrees to normalize

    Returns:
        arc-degrees in range 0 .. 360

    >>> reduce_deg(-700)
    20.0
    """
    return to_range(x, 360)


def reduce_rad(x: float) -> float:
    """Reduces x to 0 >= x < `PI2`

    Args:
        x: radians to normalize

    Returns:
        radians in range 0 .. `PI2`

    >>> round(reduce_rad(12.89), 4)
    0.3236
    """
    return to_range(x, PI2)


def frac(x: float) -> float:
    """Fractional part of a number.
    The result always keeps sign of the argument.

    Args:
        x: a number

    Returns:
        fractional part of x

    >>> frac(-5.5)
    -0.5
    """
    return modf(x)[0]


def frac360(x: float) -> float:
    """
    Used with polinomial function for better accuracy.

    Args:
        x: arc-degrees to normalize

    Returns:
        arc-degrees in range 0 .. 360

    >>> round(frac360(862.7609301843507), 4)
    273.9349
    """
    return frac(x) * 360.0


def ddd(d: int, m: int = 0, s: float = 0.0) -> float:
    """Converts sexagesimal values to decimal arc-degrees or hours.

    Args:
        d: degrees or hours
        m: minutes
        s: seconds

    Returns:
        decimal degrees or hours

    In case of negative input value only the first non-zero element will be negative.

    >>> ddd(-55, 45, 0)
    -55.75
    >>> ddd(0, -45)
    -0.75
    """

    def sum_frac(a: float, b: float) -> float:
        return fabs(a) / 60 + fabs(b)

    x = reduce(sum_frac, (s, m, d))
    # if there is at least one negative argument, the result should be negative
    try:
        next(a for a in (d, m, s) if a < 0)
    except StopIteration:
        return x
    else:
        return -x


def dms(x: float) -> tuple[int, int, float]:
    """Converts decimal hours (or degrees) to sexagesimal values.

    Args:
        x: arc-degrees or hours

    Returns:
        a tuple of (degrees, minutes, seconds) or (hours, minutes, seconds)

    >>> dms(55.75)
    (55, 45, 0.0)
    """
    f, i = modf(abs(x))
    d = int(i)
    f, i = modf(f * 60)
    m = int(i)
    s = f * 60

    if x < 0:
        if d != 0:
            d = -d
        elif m != 0:
            m = -m
        else:
            s = -s

    return (d, m, s)


def zdms(x: float) -> tuple[int, int, int, float]:
    """Converts decimal degrees to zodiac sign number, zodiac degrees, minutes and seconds.

    Args:
        x: arc-degrees

    Returns:

    * zodiac sign (0 for Aries, 11 for Pisces)
    * degrees (0-29)
    * minutes (0-59)
    * seconds (0-59.9)

    >>> zdms(320.25)  # Aquarius, 20:15
    (10, 20, 15, 0.0)
    """
    d, m, s = dms(x)
    z = d // 30
    d = d % 30
    return z, d, m, s


def shortest_arc_deg(a: float, b: float) -> float:
    """Calculate shortest arc in dergees between a and b.

    Args:
        a: the first point, in arc-degrees
        b: the second point, in arc-degrees

    Returns:
        the shortest distance, in arc-degrees

    >>> shortest_arc_deg(10, 270)
    100

    """

    x = fabs(a - b)
    return 360.0 - x if x > 180 else x


def shortest_arc_rad(a: float, b: float) -> float:
    """Calculate shortest arc in radians between a and b.

    Args:
        a: the first point, in radians
        b: the second point, in radians

    Returns:
        the shortest distance, in radians

    """
    x = fabs(a - b)
    return PI2 - x if x > pi else x


def diff_angle(a: float, b: float) -> float:
    """Angle `b - a`, accounting for circular values.

    This allows us to directly compare angles which cross through 0:
    `359 degrees... 0 degrees... 1 degree...` etc.

    Parameters a and b should be in the range `0..360`. The
    result will be in the range `-180..180`.

    Args:
        a (float): the first angle, in arc-degrees
        b (float): the second angle, in arc-degrees

    Returns:
        float: angle in arc-degrees
    """
    x = b + 360 - a if b < a else b - a
    return x - 360 if x > 180 else x


Polar = namedtuple("Polar", "phi rho")
Polar.__doc__ = """Polar coordinates.

Args:
    rho: radial coordinate.
    phi: angular coordinate.
"""
