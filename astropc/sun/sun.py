"""Position of the Sun.
"""

from functools import partial
from math import cos, degrees, floor, radians, sin
from astropc import nutation
from astropc.mathutils import PI2, Polar, frac360, polynome, reduce_deg
from astropc.kepler import eccentric_anomaly, true_anomaly
from astropc.timeutils.julian import DAYS_PER_CENT

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


_ABERRATION = 5.69e-3
"""aberration in degrees"""


def mean_longitude(t: float) -> float:
    """Mean longitude of the Sun.

    Args:
        t (float): number of Julian days elapsed since 1900, Jan 0.5.

    Returns:
        float: arc-degrees
    """
    return reduce_deg(2.7969668e2 + 3.025e-4 * t * t + frac360(1.000021359e2 * t))


def mean_anomaly(t: float) -> float:
    """Mean anomaly of the Sun

    Args:
        t (float): number of Julian days elapsed since 1900, Jan 0.5.

    Returns:
        float: arc-degrees
    """
    return reduce_deg(
        3.5847583e2 - (1.5e-4 + 3.3e-6 * t) * t * t + frac360(9.999736042e1 * t)
    )


def _calc_perturbations(t: float, a: float, b: float) -> float:
    return radians(a + frac360(b * t))


def true_geocentric(t: float, ms: float | None = None) -> Polar:
    """Calculates true geocentric position of the Sun for the mean equinox of date.

    Args:
        t (float): number of Julian days elapsed since 1900, Jan 0.5.
        ms (float | None, optional): mean anomaly of the Sun. Defaults to None.

    If ms is not provided by the caller, it will be calculated automatically.

    Returns:
        Polar: true geocentric longitude in degrees and the Sun-Earth distance in AU.
    """
    if ms is None:
        ms = mean_anomaly(t)

    ls = mean_longitude(t)
    ma = radians(ms)
    s = polynome(t, 1.675104e-2, -4.18e-5, -1.26e-7)  # eccentricity
    ea = eccentric_anomaly(s, ma - PI2 * floor(ma / PI2))  # eccentric anomaly
    nu = true_anomaly(s, ea)  # true anomaly
    t2 = t * t

    calc_pert = partial(_calc_perturbations, t)

    a = calc_pert(153.23, 6.255209472e1)  # Venus
    b = calc_pert(216.57, 1.251041894e2)  # ?
    c = calc_pert(312.69, 9.156766028e1)  # ?
    d = calc_pert(350.74 - 1.44e-3 * t2, 1.236853095e3)  # Moon
    h = calc_pert(353.4, 1.831353208e2)  # ?
    e = radians(231.19 + 20.2 * t)  # inequality of long period

    # correction in orbital longitude
    dl = (
        1.34e-3 * cos(a)
        + 1.54e-3 * cos(b)
        + 2e-3 * cos(c)
        + 1.79e-3 * sin(d)
        + 1.78e-3 * sin(e)
    )
    # correction in radius-vector
    dr = (
        5.43e-6 * sin(a)
        + 1.575e-5 * sin(b)
        + 1.627e-5 * sin(c)
        + 3.076e-5 * cos(d)
        + 9.27e-6 * sin(h)
    )
    lsn = reduce_deg(degrees(nu) + ls - ms + dl)
    rsn = 1.0000002 * (1 - s * cos(ea)) + dr
    return Polar(lsn, rsn)


def apparent(
    djd: float, dpsi: float | None = None, ignore_light_travel: bool = True
) -> Polar:
    """Apparent position of the Sun with respect of nutation,
    aberration and optionally, light-time travel.

    Args:
        djd (float): number of Julian days since 1900 Jan. 0.5.
        dpsi (float | None, optional): nutation in longitude, degrees.
        ignore_light_travel (bool, optional): ignore light travel time? Defaults to True.

    Returns:
        Polar: apparent longitude of the Sun (degrees) and its distance from Earth (AU).
    """

    t = djd / DAYS_PER_CENT

    if dpsi is None:
        dpsi = nutation.calc_nutation(t).dpsi

    tgeo = true_geocentric(t)
    lmbda = tgeo.phi + dpsi  # nutation
    lmbda -= _ABERRATION  # aberration

    if not ignore_light_travel:
        dt = 1.365 * tgeo.rho  # seconds
        lmbda -= dt * 15 / 3600  # convert to degrees and substract

    return Polar(lmbda, tgeo.rho)
