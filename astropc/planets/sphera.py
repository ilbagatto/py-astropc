"""Provides information that is outside of a planet's
class scope, yet required to calculate its position.

For instance, to calculate Mercury perturbations we need
to know mean anomalies of Venus and Jupiter.
"""

from functools import cache
from math import radians

from astropc import sun
from astropc.mathutils import Polar, reduce_rad
from astropc.nutation import Nutation, calc_nutation
from astropc.obliq import calc_obliquity
from .ids import PlanetId
from .planet import Planet
from astropc.timeutils.julian import DAYS_PER_CENT

from .orbit import OrbitInstance

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


class CelestialSphera:
    """Provides contextual information that is outside of a planet's
    class scope, yet required to calculate its position.

    For instance, to calculate Mercury perturbations, we need
    to know mean anomalies of Venus and Jupiter.
    """

    def __init__(
        self,
        t: float,
        manom_sun: float,
        sun_geo: Polar,
        obliquity: float,
        nut: Nutation,
        apparent: bool = True,
    ) -> None:
        """Initializer.

        Args:
            t (float): Number of Julian days elapsed since 1900, Jan 0.5.
            manom_sun (float): Sun mean anomaly in radians
            sun_geo (Polar): true geocentric coordinates of the Sun
            obliquity (float): obliquity of the ecliptic, degrees.
            nut (Nutation): nutation in longitude and obliquity
            apparent (bool, optional): calculate apparent positions? Defaults to True.
        """
        self._t = t
        self._manom_sun = manom_sun
        self._sun_geo = sun_geo
        self._obliquity = obliquity
        self._nut = nut
        self._apparent = apparent

        # Auxiliraly Sun-related elements needed for calculating perturbations
        aux = [0.0] * 6
        aux[0] = t / 5 + 0.1
        aux[1] = reduce_rad(4.14473 + 5.29691e1 * t)
        aux[2] = reduce_rad(4.641118 + 2.132991e1 * t)
        aux[3] = reduce_rad(4.250177 + 7.478172 * t)
        aux[4] = 5 * aux[2] - 2 * aux[1]
        aux[5] = 2 * aux[1] - 6 * aux[2] + 3 * aux[3]
        self._aux_sun = tuple(aux)

    @property
    def t(self) -> float:
        """Number of Julian days elapsed since 1900, Jan 0.5.

        Actually, this is moment of time for which all calculations are done.
        """
        return self._t

    @property
    def manom_sun(self) -> float:
        """Sun mean anomaly in radians."""
        return self._manom_sun

    @property
    def sun_geo(self) -> Polar:
        """Polar coordinates of the Sun."""
        return self._sun_geo

    @property
    def obliquity(self) -> float:
        """obliquity of the ecliptic, arc-degrees."""
        return self._obliquity

    @property
    def nutation(self) -> Nutation:
        """nutation in longitude and obliquity."""
        return self._nut

    @property
    def apparent(self) -> bool:
        """True for apparent positions."""
        return self._apparent

    def _instantiate_orbit(self, id: PlanetId) -> OrbitInstance:
        pla = Planet.for_id(id)
        return pla.orbit.instantiate(self.t)

    @cache
    def get_orbit_instance(self, id: PlanetId) -> OrbitInstance:
        """Instantiated orbit.

        Once calculated, the result for given **t** is cached.

        Args:
            id (PlanetId): planet identifier.

        Returns:
            OrbitInstance
        """
        return self._instantiate_orbit(id)

    def get_mean_anomaly(self, id: PlanetId, dt: float | None = None) -> float:
        """Mean anomaly of a planet.

        Args:
            id (PlanetId): planet identifier
            dt (float | None, optional): _time correction necessary when calculating true
               (light-time corrected) planetary positions. Defaults to None.

        Returns:
            float: mean anomaly in radians.
        """
        ma = self.get_orbit_instance(id).mean_anomaly
        if dt is not None:
            ma -= radians(dt * Planet.for_id(id).orbit.daily_motion)
        return ma

    @staticmethod
    def create(djd: float, apparent: bool = True) -> "CelestialSphera":
        """Factory that simplifies creation of instances.

        Args:
            djd (float): number of Julian days since 1900 Jan. 0.5.
            apparent (bool, optional): calculate apparent position? Defaults to True.

        Returns:
            CelestialSphera instance.
        """
        t = djd / DAYS_PER_CENT
        ms = sun.mean_anomaly(t)
        sg = sun.true_geocentric(t, ms)
        nu = calc_nutation(t)
        ob = calc_obliquity(djd, deps=nu.deps)

        return CelestialSphera(
            t=t,
            manom_sun=radians(ms),
            sun_geo=sg,
            nut=nu,
            obliquity=ob,
            apparent=apparent,
        )

    @property
    def aux_sun(self) -> tuple[float, ...]:
        return self._aux_sun
