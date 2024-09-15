"""Kepler equation
"""

from math import atan, cos, fabs, sin, sqrt, tan

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"

_DELTA = 1e-7  # precision for Kepler equation


def eccentric_anomaly(s: float, m: float, _ea: float | None = None) -> float:
    """Solve Kepler equation.

    The function calls itself recursively.

    Args:
        s: s (< 1), the eccentricity
        m: mean anomaly, radians

    Returns:
        the eccentric anomaly

    """
    if _ea is None:
        _ea = m
    dla = _ea - (s * sin(_ea)) - m
    if fabs(dla) < _DELTA:
        return _ea

    return eccentric_anomaly(s, m, _ea - dla / (1 - (s * cos(_ea))))


def true_anomaly(s: float, ea: float) -> float:
    """Find true anomaly.

    Args:
        s: eccentricity
        ea: eccentric anomaly in radians

    Returns:
        true anomaly, radians

    """
    return 2 * atan(sqrt((1 + s) / (1 - s)) * tan(ea / 2))
