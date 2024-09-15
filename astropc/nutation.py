"""Effects of nutation on the ecliptic longitude, dpsi, and on
the obliquity of the ecliptic, deps, with accuracy of about 1 arcsecond.

Source: P.Duffett-Smith, "Astronomy with Your PC", 2 edition.
"""

from collections import namedtuple
from math import cos, radians, sin

from astropc.mathutils import frac360

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"

Nutation = namedtuple("Nutation", "dpsi deps")
Nutation.__doc__ = """Nutation on the ecliptic longitude, dpsi, and on
the obliquity of the ecliptic, deps.
"""


def calc_nutation(t: float) -> Nutation:
    """Calculates effects of nutation.

    Args:
        t (float): number of Julian days elapsed since 1900, Jan 0.5.

    Returns:
        Nutation record.
    """
    t2 = t * t

    ls = radians(2.796967e2 + 3.030e-4 * t2 + frac360(1.000021358e2 * t))
    ms = radians(3.584758e2 - 1.500e-4 * t2 + frac360(9.999736056e1 * t))
    ld = radians(2.704342e2 - 1.133e-3 * t2 + frac360(1.336855231e3 * t))
    md = radians(2.961046e2 + 9.192e-3 * t2 + frac360(1.325552359e3 * t))
    nm = radians(2.591833e2 + 2.078e-3 * t2 - frac360(5.372616667 * t))
    tls = ls + ls
    tld = ld + ld
    tnm = nm + nm

    dpsi = (
        (-17.2327 - 1.737e-2 * t) * sin(nm)
        + (-1.2729 - 1.3e-4 * t) * sin(tls)
        + 2.088e-1 * sin(tnm)
        - 2.037e-1 * sin(tld)
        + (1.261e-1 - 3.1e-4 * t) * sin(ms)
        + 6.75e-2 * sin(md)
        - (4.97e-2 - 1.2e-4 * t) * sin(tls + ms)
        - 3.42e-2 * sin(tld - nm)
        - 2.61e-2 * sin(tld + md)
        + 2.14e-2 * sin(tls - ms)
        - 1.49e-2 * sin(tls - tld + md)
        + 1.24e-2 * sin(tls - nm)
        + 1.14e-2 * sin(tld - md)
    )

    deps = (
        (9.21 + 9.1e-4 * t) * cos(nm)
        + (5.522e-1 - 2.9e-4 * t) * cos(tls)
        - 9.04e-2 * cos(tnm)
        + 8.84e-2 * cos(tld)
        + 2.16e-2 * cos(tls + ms)
        + 1.83e-2 * cos(tld - nm)
        + 1.13e-2 * cos(tld + md)
        - 9.3e-3 * cos(tls - ms)
        - 6.6e-3 * cos(tls - nm)
    )

    return Nutation(
        dpsi / 3600, deps / 3600
    )  # 1965-2-1 11:46 dpsi = -0.0042774118548615766
