"""Base type and factory of the planets.
"""

from collections import namedtuple
from dataclasses import dataclass
from functools import cache
from math import asin, atan, atan2, cos, degrees, pi, radians, sin, sqrt
from typing import TYPE_CHECKING

from astropc.kepler import eccentric_anomaly, true_anomaly
from astropc.mathutils import reduce_rad

from .ids import PlanetId
from .pert import (
    PertCalculator,
    PertJupiter,
    PertMars,
    PertMercury,
    PertNeptune,
    PertPluto,
    PertSaturn,
    PertUranus,
    PertVenus,
)

if TYPE_CHECKING:
    from .sphera import CelestialSphera

from .orbit import MLTerms, OElements, OrbitInstance, Terms
from .pert import PertRecord

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


HelioRecord = namedtuple("HelioRecord", "ll rpd lpd spsi cpsi rho")
HelioRecord.__doc__ = """Params of calculated planetary heliocentric orbit."""


@dataclass
class EclipticPosition:
    """Ecliptic posiion of a celestial body."""

    lmbda: float = 0.0
    """geocentric ecliptic longitude, arc-degrees"""

    beta: float = 0.0
    """geocentric ecliptic latitude, arc-degrees"""

    delta: float = 0.0
    """distance from the Earth. A.U."""


class Planet:
    """Base class for planets."""

    def __init__(
        self,
        id: PlanetId,
        name: str,
        is_inner: bool,
        orbit: OElements,
        pert_calculator: PertCalculator,
    ) -> None:
        self._id = id
        self._name = name
        self._is_inner = is_inner
        self._orbit = orbit
        self._pert_calculator = pert_calculator

    @property
    def id(self) -> PlanetId:
        """Planet identifier."""
        return self._id

    @property
    def name(self) -> str:
        """Planet name."""
        return self._name

    @property
    def is_inner(self) -> bool:
        """True if the planet is inner."""
        return self._is_inner

    @property
    def orbit(self) -> OElements:
        """Osculating elements of the orbit."""
        return self._orbit

    @property
    def pert_calculator(self) -> PertCalculator:
        """Object responsible for calculating perturbations."""
        return self._pert_calculator

    def _calculate_perturbations(
        self, ctx: "CelestialSphera", dt: float = 0
    ) -> PertRecord:
        """Calculate perturbations."""
        return self.pert_calculator.get_data(ctx, dt)

    @staticmethod
    def _calculate_heliocentric(
        oi: OrbitInstance, ma: float, re: float, lg: float, pert: PertRecord
    ) -> HelioRecord:
        """Core part of heliocentric position calculation.

        Args:
            oi: osculating elements of the obit instantiated for the moment.
            ma: mean anomaly of the planet
            re: Sun-Earth distance
            lg: lonitude of the Earth, radians
            pert: PertRecord instance

        Returns:
            HelioRecord: parameters of the heliocentric orbit.
        """
        s = oi.eccentricity + pert.ds  # eccentricity corrected
        ma = reduce_rad(ma + pert.dm)  # mean anomaly corrected
        ea = eccentric_anomaly(s, ma)  # eccentric anomaly
        nu = true_anomaly(s, ea)  # true anomaly
        # radius-vector
        rp = (oi.major_semiaxis + pert.da) * (1 - s * s) / (1 + s * cos(nu)) + pert.dr
        lp = nu + oi.perihelion + (pert.dml - pert.dm)  # planet's orbital longitude
        lo = lp - oi.mean_node
        sin_lo = sin(lo)
        spsi = sin_lo * sin(oi.inclination)
        psi = sin_lo * sin(oi.inclination)
        y = sin_lo * cos(oi.inclination)
        psi = asin(spsi) + pert.dhl  # heliocentric latitude
        lpd = atan2(y, cos(lo)) + oi.mean_node + radians(pert.dl)
        cpsi = cos(psi)
        ll = lpd - lg

        # distance from the Earth
        rho = sqrt(re * re + rp * rp - 2 * re * rp * cpsi * cos(ll))

        return HelioRecord(
            ll,
            rp * cpsi,
            lpd,
            sin(psi),  # not the same as spsi, for now psi is corrected
            cpsi,
            rho,
        )

    def _get_corrected_helio(
        self,
        ctx: "CelestialSphera",
        oi: OrbitInstance,
        lg: float,
        rg: float,
        dt: float = 0,
        rho: float = 0,
    ) -> HelioRecord:
        """Calculate heliocentric position taking account of the finit light-travel
        time between the Earth and the planet.

        This method is recursive.

            When we view a planet now, we see it in the position it occupied t
            hours ago, given by *t = 0.1386 x RH*, where RH is the distance in AU
            between the Earth and the planet. In this routine, an approximate position
            for the planet is first calculated, neglecting the light-travel time.
            Then a second pass is made through the program using the light-travel
            time based on the approximate position found on the first pass.

            -- Peter Duffett-Smith, p.137-138
        """

        ma = ctx.get_mean_anomaly(self.id, dt)
        pert = self._calculate_perturbations(ctx, dt)
        h = self._calculate_heliocentric(oi, ma, rg, lg, pert)
        if dt == 0:
            # take account of the finit light-travel time between the Earth and the planet.
            # h.rho is the Earth-planet distance
            return self._get_corrected_helio(
                ctx, oi, lg, rg, h.rho * 5.775518e-3, h.rho
            )

        return HelioRecord(
            h.ll,
            h.rpd,
            h.lpd,
            h.spsi,  # not the same as spsi, for now psi is corrected
            h.cpsi,
            rho,
        )

    def geocentric_position(self, ctx: "CelestialSphera") -> "EclipticPosition":
        """Calculate geocentric position of a planet.

        Args:
            ctx: context, containing some key attributes of the moment.

        Returns:
            EclipticPosition: record containing the ecliptic position parameters.
        """
        sg = ctx.sun_geo
        # convert logitude of the Sun to Earth's position
        lg = radians(sg.phi) + pi
        rsn = sg.rho  # Sun-Earth distance
        oi = self.orbit.instantiate(ctx.t)
        # heliocentric position corrected for light-time travel
        h = self._get_corrected_helio(ctx, oi, lg, rsn)

        # Convert to geocentric
        sll = sin(h.ll)
        cll = cos(h.ll)
        # geocentric ecliptic longitude
        lam = (
            atan2(-1 * h.rpd * sll, rsn - h.rpd * cll) + lg + pi
            if self.is_inner
            else atan2(rsn * sll, h.rpd - rsn * cll) + h.lpd
        )
        lam = reduce_rad(lam)
        # geocentric latitude
        bet = atan(h.rpd * h.spsi * sin(lam - h.lpd) / (h.cpsi * rsn * sll))

        if ctx.apparent:
            # nutation
            lam += radians(ctx.nutation.dpsi)
            # aberration
            a = lg - lam
            lam -= 9.9387e-5 * cos(a) / cos(bet)
            lam = reduce_rad(lam)
            bet -= 9.9387e-5 * sin(a) * sin(bet)

        return EclipticPosition(lmbda=degrees(lam), beta=degrees(bet), delta=h.rho)

    @staticmethod
    @cache
    def for_id(id: PlanetId) -> "Planet":
        """Factory method for creating a planet instance.

        Once created, the object is cached.

        Args:
            id (PlanetId): identifier

        Returns:
            Planet instance.
        """
        match (id):
            case PlanetId.MERCURY:
                return Planet(
                    id,
                    "Mercury",
                    True,
                    OElements(
                        mean_longitude=MLTerms(178.179078, 415.2057519, 3.011e-4),
                        perihelion=Terms(75.899697, 1.5554889, 2.947e-4),
                        eccentricity=Terms(2.0561421e-1, 2.046e-5, -3e-8),
                        inclination=Terms(7.002881, 1.8608e-3, -1.83e-5),
                        mean_node=Terms(47.145944, 1.1852083, 1.739e-4),
                        major_semiaxis=3.870986e-1,
                    ),
                    PertMercury(),
                )
            case PlanetId.VENUS:
                return Planet(
                    id,
                    "Venus",
                    True,
                    OElements(
                        mean_longitude=MLTerms(342.767053, 162.5533664, 3.097e-4),
                        perihelion=Terms(130.163833, 1.4080361, -9.764e-4),
                        eccentricity=Terms(6.82069e-3, -4.774e-5, 9.1e-8),
                        inclination=Terms(3.393631, 1.0058e-3, -1e-6),
                        mean_node=Terms(75.779647, 8.9985e-1, 4.1e-4),
                        major_semiaxis=7.233316e-1,
                    ),
                    PertVenus(),
                )
            case PlanetId.MARS:
                return Planet(
                    id,
                    "Mars",
                    False,
                    OElements(
                        mean_longitude=MLTerms(293.737334, 53.17137642, 3.107e-4),
                        perihelion=Terms(3.34218203e2, 1.8407584, 1.299e-4, -1.19e-6),
                        eccentricity=Terms(9.33129e-2, 9.2064e-5, -7.7e-8),
                        inclination=Terms(1.850333, -6.75e-4, 1.26e-5),
                        mean_node=Terms(48.786442, 7.709917e-1, -1.4e-6, -5.33e-6),
                        major_semiaxis=1.5236883,
                    ),
                    PertMars(),
                )
            case PlanetId.JUPITER:
                return Planet(
                    id,
                    "Jupiter",
                    False,
                    OElements(
                        mean_longitude=MLTerms(
                            238.049257, 8.434172183, 3.347e-4, -1.65e-6
                        ),
                        perihelion=Terms(1.2720972e1, 1.6099617, 1.05627e-3, -3.43e-6),
                        eccentricity=Terms(4.833475e-2, 1.6418e-4, -4.676e-7, -1.7e-9),
                        inclination=Terms(1.308736, -5.6961e-3, 3.9e-6),
                        mean_node=Terms(99.443414, 1.01053, 3.5222e-4, -8.51e-6),
                        major_semiaxis=5.202561,
                    ),
                    PertJupiter(),
                )
            case PlanetId.SATURN:
                return Planet(
                    id,
                    "Saturn",
                    False,
                    OElements(
                        mean_longitude=MLTerms(
                            266.564377, 3.398638567, 3.245e-4, -5.8e-6
                        ),
                        perihelion=Terms(9.1098214e1, 1.9584158, 8.2636e-4, 4.61e-6),
                        eccentricity=Terms(5.589232e-2, -3.455e-4, -7.28e-7, 7.4e-10),
                        inclination=Terms(2.492519, -3.9189e-3, -1.549e-5, 4e-8),
                        mean_node=Terms(112.790414, 8.731951e-1, -1.5218e-4, -5.31e-6),
                        major_semiaxis=9.554747,
                    ),
                    PertSaturn(),
                )
            case PlanetId.URANUS:
                return Planet(
                    id,
                    "Uranus",
                    False,
                    OElements(
                        mean_longitude=MLTerms(244.19747, 1.194065406, 3.16e-4, -6e-7),
                        perihelion=Terms(1.71548692e2, 1.4844328, 2.372e-4, -6.1e-7),
                        eccentricity=Terms(4.63444e-2, -2.658e-5, 7.7e-8),
                        inclination=Terms(7.72464e-1, 6.253e-4, 3.95e-5),
                        mean_node=Terms(73.477111, 4.986678e-1, 1.3117e-3),
                        major_semiaxis=19.21814,
                    ),
                    PertUranus(),
                )
            case PlanetId.NEPTUNE:
                return Planet(
                    id,
                    "Neptune",
                    False,
                    OElements(
                        mean_longitude=MLTerms(
                            84.457994, 6.107942056e-1, 3.205e-4, -6e-7
                        ),
                        perihelion=Terms(4.6727364e1, 1.4245744, 3.9082e-4, -6.05e-7),
                        eccentricity=Terms(8.99704e-3, 6.33e-6, -2e-9),
                        inclination=Terms(1.779242, -9.5436e-3, -9.1e-6),
                        mean_node=Terms(130.681389, 1.098935, 2.4987e-4, -4.718e-6),
                        major_semiaxis=30.10957,
                    ),
                    PertNeptune(),
                )
            case PlanetId.PLUTO:
                return Planet(
                    id,
                    "Pluto",
                    False,
                    OElements(
                        mean_longitude=MLTerms(95.3113544, 3.980332167e-1),
                        perihelion=Terms(224.017),
                        eccentricity=Terms(2.5515e-1),
                        inclination=Terms(17.1329),
                        mean_node=Terms(110.191),
                        major_semiaxis=39.8151,
                    ),
                    PertPluto(),
                )
