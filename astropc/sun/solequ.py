from collections import namedtuple
from enum import Enum, auto
from math import radians, sin

from astropc.mathutils import shortest_arc_deg

from .sun import apparent

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


class SolEquType(Enum):
    """Event types"""

    MARCH_EQUINOX = auto()
    JUNE_SOLSTICE = auto()
    SEPTEMBER_EQUINOX = auto()
    DECEMBER_SOLSTICE = auto()


_DELTA = 1e-6


SolEquEvent = namedtuple("SolEquEvent", "djd sun")
SolEquEvent.__doc__ = """Solstice/quinox circumstances.
djd is time of the event, number of Julian days since 1900 Jan. 0.5.
sun is apparent longitude of the Sun, arc-degrees
"""


def sol_equ(year: int, event_type: SolEquType) -> SolEquEvent:
    """Calculate circumstances of Solstice/Equinox event.

    Args:
        year (int): year
        event_type (SolEquType): event type

    Returns:
        SolEquEvent: time of the event and longitude of the Sun.
    """
    match event_type:
        case SolEquType.MARCH_EQUINOX:
            k = 0
        case SolEquType.JUNE_SOLSTICE:
            k = 1
        case SolEquType.SEPTEMBER_EQUINOX:
            k = 2
        case SolEquType.DECEMBER_SOLSTICE:
            k = 3
    k90 = k * 90
    dj = (year + k / 4.0) * 365.2422 - 693878.7  # shorter but less exact way
    while True:
        sg = apparent(dj, ignore_light_travel=True)
        x = sg.phi
        dj += 58.0 * sin(radians(k90 - x))
        if shortest_arc_deg(k90, x) < _DELTA:
            break
    return SolEquEvent(dj, x)
