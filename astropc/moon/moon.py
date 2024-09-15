"""The Moon and Lunar Node.
"""

from collections import namedtuple
from functools import partial
from math import cos, radians, sin

from astropc.nutation import calc_nutation
from astropc.mathutils import frac, polynome, reduce_deg
from astropc.timeutils.julian import DAYS_PER_CENT

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


MoonPosition = namedtuple("MoonPosition", "lmbda beta delta parallax motion")
MoonPosition.__doc__ = (
    """Moon ecliptic coordinates (lambda, beta), parallax and daily motion."""
)

_M = (27.32158213, 365.2596407, 27.55455094, 29.53058868, 27.21222039, 6798.363307)

_MOON_ORBIT = {
    # Mean longitude
    "L": (218.3164477, 481267.88123421, -0.0015786, 1.0 / 538841, -(1.0 / 65194000)),
    # Mean elongation
    "D": (297.8501921, 445267.1114034, -0.0018819, 1.0 / 545868, -(1.0 / 113065000)),
    # Mean anomaly
    "M": (134.9633964, 477198.8675055, 0.0087414, 1.0 / 69699, -(1.0 / 14712000)),
    # Argument of latitude (mean distance of the Moon from its ascending node)
    "F": (93.272095, 483202.0175233, -0.0036539, -(1.0 / 3526000), 1.0 / 863310000),
}

_SUN_ORBIT = {
    # Mean anomaly
    "M": (357.5291092, 35999.0502909, -0.0001536, 1.0 / 24490000)
}


def mean_node(t: float) -> float:
    """Mean Lunar Node.

    Args:
        t: number of Julian centuries elapsed since 1900, Jan 0.5.

    Returns:
        longitude in degrees
    """
    return reduce_deg(
        polynome(t, 125.0445479, -1934.1362891, 0.0020754, 1.0 / 467441, 1.0 / 60616000)
    )


def true_position(djd: float) -> MoonPosition:
    """True position of the Moon.

    Args:
        djd: number of Julian days since 1900 Jan. 0.5.

    Returns:
        `MoonPosition` record.
    """

    t = djd / 36525
    t2 = t * t

    m = [360 * frac(djd / x) for x in _M]

    ld = 270.434164 + m[0] - (1.133e-3 - 1.9e-6 * t) * t2  # Moon's mean longitude
    ms = 358.475833 + m[1] - (1.5e-4 + 3.3e-6 * t) * t2  # mean anomaly of the Sun
    md = 296.104608 + m[2] + (9.192e-3 + 1.44e-5 * t) * t2  # mean anomaly
    de = 350.737486 + m[3] - (1.436e-3 - 1.9e-6 * t) * t2  # mean elongation
    f = (
        11.250889 + m[4] - (3.211e-3 + 3e-7 * t) * t2
    )  # mean distance of Moon from its ascending node
    n = (
        259.183275 - m[5] + (2.078e-3 + 2.2e-5 * t) * t2
    )  # longitude of Moon's asc. node
    a = radians(51.2 + 20.2 * t)
    sa = sin(a)
    sn = sin(radians(n))
    b = 346.56 + (132.87 - 9.1731e-3 * t) * t
    sb = 3.964e-3 * sin(radians(b))
    c = radians(n + 275.05 - 2.3 * t)
    sc = sin(c)
    ld += 2.33e-4 * sa + sb + 1.964e-3 * sn
    ms -= 1.778e-3 * sa
    md += 8.17e-4 * sa + sb + 2.541e-3 * sn
    f += sb - 2.4691e-2 * sn - 4.328e-3 * sc
    de += 2.011e-3 * sa + sb + 1.964e-3 * sn
    e = 1 - (2.495e-3 + 7.52e-6 * t) * t
    e2 = e * e
    ms = radians(ms)
    n = radians(n)
    de = radians(de)
    f = radians(f)
    md = radians(md)

    de2 = de + de
    de3 = de2 + de
    de4 = de2 + de2
    md2 = md + md
    md3 = md2 + md
    ms2 = ms + ms
    f2 = f + f
    f3 = f2 + f
    #
    # ecliptic longitude
    #
    l = (  # noqa: E741
        6.28875 * sin(md)
        + 1.274018 * sin(de2 - md)
        + 6.58309e-1 * sin(de2)
        + 2.13616e-1 * sin(md2)
        - e * 1.85596e-1 * sin(ms)
        - 1.14336e-1 * sin(f2)
        + 5.8793e-2 * sin(2 * (de - md))
        + 5.7212e-2 * e * sin(de2 - ms - md)
        + 5.332e-2 * sin(de2 + md)
        + 4.5874e-2 * e * sin(de2 - ms)
        + 4.1024e-2 * e * sin(md - ms)
        - 3.4718e-2 * sin(de)
        - e * 3.0465e-2 * sin(ms + md)
        + 1.5326e-2 * sin(2 * (de - f))
        - 1.2528e-2 * sin(f2 + md)
        - 1.098e-2 * sin(f2 - md)
        + 1.0674e-2 * sin(de4 - md)
        + 1.0034e-2 * sin(md3)
        + 8.548e-3 * sin(de4 - md2)
        - e * 7.91e-3 * sin(ms - md + de2)
        - e * 6.783e-3 * sin(de2 + ms)
        + 5.162e-3 * sin(md - de)
        + e * 5e-3 * sin(ms + de)
        + 3.862e-3 * sin(de4)
        + e * 4.049e-3 * sin(md - ms + de2)
        + 3.996e-3 * sin(2 * (md + de))
        + 3.665e-3 * sin(de2 - md3)
        + e * 2.695e-3 * sin(md2 - ms)
        + 2.602e-3 * sin(md - 2 * (f + de))
        + e * 2.396e-3 * sin(2 * (de - md) - ms)
        - 2.349e-3 * sin(md + de)
        + e2 * 2.249e-3 * sin(2 * (de - ms))
        - e * 2.125e-3 * sin(md2 + ms)
        - e2 * 2.079e-3 * sin(ms2)
        + e2 * 2.059e-3 * sin(2 * (de - ms) - md)
        - 1.773e-3 * sin(md + 2 * (de - f))
        - 1.595e-3 * sin(2 * (f + de))
        + e * 1.22e-3 * sin(de4 - ms - md)
        - 1.11e-3 * sin(2 * (md + f))
        + 8.92e-4 * sin(md - de3)
        - e * 8.11e-4 * sin(ms + md + de2)
        + e * 7.61e-4 * sin(de4 - ms - md2)
        + e2 * 7.04e-4 * sin(md - 2 * (ms + de))
        + e * 6.93e-4 * sin(ms - 2 * (md - de))
        + e * 5.98e-4 * sin(2 * (de - f) - ms)
        + 5.5e-4 * sin(md + de4)
        + 5.38e-4 * sin(4 * md)
        + e * 5.21e-4 * sin(de4 - ms)
        + 4.86e-4 * sin(md2 - de)
        + e2 * 7.17e-4 * sin(md - ms2)
    )

    lmbda = reduce_deg(ld + l)

    #
    # ecliptic latitude
    #
    g = (
        5.128189 * sin(f)
        + 0.280606 * sin(md + f)
        + 0.277693 * sin(md - f)
        + 0.173238 * sin(de2 - f)
        + 0.055413 * sin(de2 + f - md)
        + 0.046272 * sin(de2 - f - md)
        + 0.032573 * sin(de2 + f)
        + 0.017198 * sin(md2 + f)
        + 0.009267 * sin(de2 + md - f)
        + 0.008823 * sin(md2 - f)
        + e * 0.008247 * sin(de2 - ms - f)
        + 0.004323 * sin(2 * (de - md) - f)
        + 0.0042 * sin(de2 + f + md)
        + e * 0.003372 * sin(f - ms - de2)
        + e * 0.002472 * sin(de2 + f - ms - md)
        + e * 0.002222 * sin(de2 + f - ms)
        + e * 0.002072 * sin(de2 - f - ms - md)
        + e * 0.001877 * sin(f - ms + md)
        + 0.001828 * sin(de4 - f - md)
        - e * 0.001803 * sin(f + ms)
        - 0.00175 * sin(f3)
        + e * 0.00157 * sin(md - ms - f)
        - 0.001487 * sin(f + de)
        - e * 0.001481 * sin(f + ms + md)
        + e * 0.001417 * sin(f - ms - md)
        + e * 0.00135 * sin(f - ms)
        + 0.00133 * sin(f - de)
        + 0.001106 * sin(f + md3)
        + 0.00102 * sin(de4 - f)
        + 0.000833 * sin(f + de4 - md)
        + 0.000781 * sin(md - f3)
        + 0.00067 * sin(f + de4 - md2)
        + 0.000606 * sin(de2 - f3)
        + 0.000597 * sin(2 * (de + md) - f)
        + e * 0.000492 * sin(de2 + md - ms - f)
        + 0.00045 * sin(2 * (md - de) - f)
        + 0.000439 * sin(md3 - f)
        + 0.000423 * sin(f + 2 * (de + md))
        + 0.000422 * sin(de2 - f - md3)
        - e * 0.000367 * sin(ms + f + de2 - md)
        - e * 0.000353 * sin(ms + f + de2)
        + 0.000331 * sin(f + de4)
        + e * 0.000317 * sin(de2 + f - ms + md)
        + e2 * 0.000306 * sin(2 * (de - ms) - f)
        - 0.000283 * sin(md + f3)
    )
    w1 = 0.0004664 * cos(n)
    w2 = 0.0000754 * cos(c)
    beta = g * (1 - w1 - w2)

    # horizontal parallax
    hp = (
        0.950724
        + 0.051818 * cos(md)
        + 0.009531 * cos(de2 - md)
        + 0.007843 * cos(de2)
        + 0.002824 * cos(md2)
        + 0.000857 * cos(de2 + md)
        + e * 0.000533 * cos(de2 - ms)
        + e * 0.000401 * cos(de2 - md - ms)
        + e * 0.00032 * cos(md - ms)
        - 0.000271 * cos(de)
        - e * 0.000264 * cos(ms + md)
        - 0.000198 * cos(f2 - md)
        + 0.000173 * cos(md3)
        + 0.000167 * cos(de4 - md)
        - e * 0.000111 * cos(ms)
        + 0.000103 * cos(de4 - md2)
        - 0.000084 * cos(md2 - de2)
        - e * 0.000083 * cos(de2 + ms)
        + 0.000079 * cos(de2 + md2)
        + 0.000072 * cos(de4)
        + e * 0.000064 * cos(de2 - ms + md)
        - e * 0.000063 * cos(de2 + ms - md)
        + e * 0.000041 * cos(ms + de)
        + e * 0.000035 * cos(md2 - ms)
        - 0.000033 * cos(md3 - de2)
        - 0.00003 * cos(md + de)
        - 0.000029 * cos(2 * (f - de))
        - e * 0.000029 * cos(md2 + ms)
        + e2 * 0.000026 * cos(2 * (de - ms))
        - 0.000023 * cos(2 * (f - de) + md)
        + e * 0.000019 * cos(de4 - ms - md)
    )

    # distance from Earth in A.U.
    delta = 8.794 / (hp * 3600)

    # angular speed
    dm = (
        13.176397
        + 1.434006 * cos(md)
        + 0.280135 * cos(de2)
        + 0.251632 * cos(de2 - md)
        + 0.097420 * cos(md2)
        - 0.052799 * cos(f2)
        + 0.034848 * cos(de2 + md)
        + 0.018732 * cos(de2 - ms)
        + 0.010316 * cos(de2 - ms - md)
        + 0.008649 * cos(ms - md)
        - 0.008642 * cos(f2 + md)
        - 0.007471 * cos(ms + md)
        - 0.007387 * cos(de)
        + 0.006864 * cos(md2 + md)
        + 0.006650 * cos(de4 - md)
        + 0.003523 * cos(de2 + md2)
        + 0.003377 * cos(de4 - md2)
        + 0.003287 * cos(de4)
        - 0.003193 * cos(ms)
        - 0.003003 * cos(de2 + ms)
        + 0.002577 * cos(md - ms + de2)
        - 0.002567 * cos(f2 - md)
        - 0.001794 * cos(de2 - md2)
        - 0.001716 * cos(md - f2 - de2)
        - 0.001698 * cos(de2 + ms - md)
        - 0.001415 * cos(de2 + f2)
        + 0.001183 * cos(md2 - ms)
        + 0.001150 * cos(de + ms)
        - 0.001035 * cos(de + md)
        - 0.001019 * cos(f2 + md2)
        - 0.001006 * cos(ms + md2)
    )

    return MoonPosition(lmbda, beta, delta, hp, dm)


def apparent(
    djd: float, dpsi: float | None = None, ignore_light_travel: bool = True
) -> MoonPosition:
    """Apparent position of the Moon,
    with respect of nutation, aberration and optionally, light-time travel.

    Args:
        djd (float): number of Julian days since 1900 Jan. 0.5.
        dpsi (float | None, optional): nutation in longitude, arc-degrees. Defaults to None.
        ignore_light_travel (bool, optional): ignore light travel time? Defaults to True.

    Returns:
        `MoonPosition` record.
    """
    true_pos = true_position(djd)
    if dpsi is None:
        dpsi = calc_nutation(djd / DAYS_PER_CENT).dpsi

    return MoonPosition(
        true_pos.lmbda + dpsi,
        true_pos.beta,
        true_pos.delta,
        true_pos.parallax,
        true_pos.motion,
    )


def _assemble(t: float, terms: tuple[float, ...]) -> float:
    return radians(reduce_deg(polynome(t, *terms)))


def lunar_node(djd: float, true_node: bool = True) -> float:
    """Longitude of Lunar Node.

    Args:
        djd (float): number of Julian days since 1900 Jan. 0.5.
        true_node: If True, then the result will refer to the *true equinox of the date*.
        Defaults to True.

    Returns:
        float: longitude in arc-degrees.
    """
    t = (djd - 36525) / 36525  # convert DJD to centuries since epoch 2000.0

    mn = polynome(
        t, 125.0445479, -1934.1362891, 0.0020754, 1.0 / 467441, 1.0 / 60616000
    )

    if true_node:
        assemble_with_t = partial(_assemble, t)
        d = assemble_with_t(_MOON_ORBIT["D"])
        m = assemble_with_t(_MOON_ORBIT["M"])
        f = assemble_with_t(_MOON_ORBIT["F"])
        ms = assemble_with_t(_SUN_ORBIT["M"])
        nd = (
            mn
            - 1.4979 * sin(2 * (d - f))
            - 0.1500 * sin(ms)
            - 0.1226 * sin(2 * d)
            + 0.1176 * sin(2 * f)
            - 0.0801 * sin(2 * (m - f))
        )
    else:
        nd = mn

    return reduce_deg(nd)
