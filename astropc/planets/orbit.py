"""Osculating orbital elements"""

from dataclasses import dataclass
from math import radians

from astropc.mathutils import frac360, polynome, reduce_deg
from astropc.timeutils.julian import DAYS_PER_CENT

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


class Terms:
    """Orbital terms"""

    def __init__(self, *args: float) -> None:
        self._data = args

    @property
    def terms(self) -> tuple[float, ...]:
        return self._data

    def assemble(self, t: float) -> float:
        """Instaniate osculating elements of an orbit.

        Args:
            t (float): time in centuries from epoch 1900.0

        Returns:
            float: arc-degrees
        """
        return reduce_deg(polynome(t, *self._data))


class MLTerms(Terms):
    """Mean Lonitude, a special case of Terms."""

    def __init__(self, a: float, b: float = 0, c: float = 0, d: float = 0) -> None:
        super().__init__(a, b, c, d)

    def assemble(self, t: float) -> float:
        """The mean longitude increases by 360 deg. for every rotation of the PlanetId
        about the Sun.

        In order to preserve accuracy, it is is expressed in such a manner that integer
        rotations are subtracted from the second term of the expression  before adding
        the other terms.

        Args:
            t (float): time in centuries from epoch 1900.0

        Returns:
            float: arc-degrees
        """
        b = frac360(self._data[1] * t)
        return reduce_deg(
            self._data[0] + b + (self._data[3] * t + self._data[2]) * t * t
        )


@dataclass
class OrbitInstance:
    """A record holding an orbit instantiated for a given moment of time.

    All angular values are in radians.
    """

    perihelion: float
    """argument of perihelion"""

    eccentricity: float
    """eccentricity"""

    mean_node: float
    """mean ascending node"""

    inclination: float
    """inclination"""

    major_semiaxis: float
    """major semi-axis"""

    mean_anomaly: float
    """mean anomaly"""

    daily_motion: float
    """mean daily motion"""


class OElements:
    """Osculating elements of an orbit"""

    def __init__(
        self,
        mean_longitude: MLTerms,
        perihelion: Terms,
        eccentricity: Terms,
        inclination: Terms,
        mean_node: Terms,
        major_semiaxis: float,
    ) -> None:
        self._ml = mean_longitude
        self._ph = perihelion
        self._ec = eccentricity
        self._in = inclination
        self._nd = mean_node
        self._sa = major_semiaxis
        self._dm = (
            self._ml.terms[1] * 9.856263e-3
            + (self._ml.terms[2] + self._ml.terms[3]) / DAYS_PER_CENT
        )

    @property
    def mean_longitude(self) -> MLTerms:
        """mean longitude terms"""
        return self._ml

    @property
    def perihelion(self) -> Terms:
        """argument of perihelion"""
        return self._ph

    @property
    def eccentricity(self) -> Terms:
        """eccentricity"""
        return self._ec

    @property
    def inclination(self) -> Terms:
        """inclination"""
        return self._in

    @property
    def mean_node(self) -> Terms:
        """mean ascending node"""
        return self._nd

    @property
    def major_semiaxis(self) -> float:
        return self._sa

    @property
    def daily_motion(self) -> float:
        """Mean Daily motion"""
        return self._dm

    def assemble_mean_anomaly(self, t: float) -> float:
        """Assemble Mean Anomaly terms.

        Args:
            t (float): time in centuries from epoch 1900.0

        Returns:
            float: Mean Anomaly in arc-degrees
        """
        return reduce_deg(self.mean_longitude.assemble(t) - self.perihelion.assemble(t))

    def instantiate(self, t: float) -> OrbitInstance:
        """Instantiate orbit for a given moment.

        Args:
            t: centuries passed since the epoch 1900,0

        Returns:
            OrbitInstance object
        """
        ph = self._ph.assemble(t)
        return OrbitInstance(
            perihelion=radians(ph),
            eccentricity=self.eccentricity.assemble(t),
            mean_node=radians(self.mean_node.assemble(t)),
            inclination=radians(self.inclination.assemble(t)),
            major_semiaxis=self.major_semiaxis,
            mean_anomaly=radians(self.assemble_mean_anomaly(t)),
            daily_motion=radians(self.daily_motion),
        )
