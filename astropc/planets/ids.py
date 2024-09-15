from enum import StrEnum

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


class PlanetId(StrEnum):
    """Identifiers of the 8 planets."""

    MERCURY = "Mercury"
    VENUS = "Venus"
    MARS = "Mars"
    JUPITER = "Jupiter"
    SATURN = "Saturn"
    URANUS = "Uranus"
    NEPTUNE = "Neptune"
    PLUTO = "Pluto"
