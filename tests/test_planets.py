from pytest import approx, fixture

from astropc.planets import CelestialSphera, EclipticPosition, Planet, PlanetId
from astropc.planets.orbit import OrbitInstance
from astropc.planets.pert import PertRecord


class TestHeliocentric:
    delta = 1e-6

    @fixture()
    def helio(self):
        oi = OrbitInstance(
            eccentricity=0.20563138612828713,
            major_semiaxis=0.3870986,
            perihelion=1.34752240012577,
            mean_node=0.8402412010285969,
            inclination=0.12225040301524157,
            mean_anomaly=1.7277480419370512,
            daily_motion=0,  # not used
        )
        re = 0.9839698373786032
        lg = 8.379862816965847
        ma = 1.7277480419370512
        pert = PertRecord(
            dl=-0.001369031774972449,
            dr=-0.000013447032242762032,
            dml=0,
            ds=0,
            dm=0,
            da=0,
            dhl=0,
        )
        return Planet._calculate_heliocentric(oi, ma, re, lg, pert)

    def test_ll(self, helio):
        assert approx(helio.ll, abs=self.delta) == -4.920247322912694

    def test_rpd(self, helio):
        assert approx(helio.rpd, abs=self.delta) == 0.4136118768849629

    def test_lpd(self, helio):
        assert approx(helio.lpd, abs=self.delta) == 3.459615494053153

    def test_spsi(self, helio):
        assert approx(helio.spsi, abs=self.delta) == 0.06116731001819705

    def test_cpsi(self, helio):
        assert approx(helio.cpsi, abs=self.delta) == 0.9981275270150292

    def test_rho(self, helio):
        assert approx(helio.rho, abs=self.delta) == 0.9858704400566043


class TestGeocentric:
    delta = 1e-4

    @fixture()
    def sphera(self):
        return CelestialSphera.create(30700.5, apparent=False)

    def _compare(self, id, exp, sphera):
        # EclipticPosition
        pla = Planet.for_id(id)
        pos = pla.geocentric_position(sphera)

        assert approx(pos.lmbda, abs=self.delta) == exp.lmbda
        assert approx(pos.beta, abs=self.delta) == exp.beta
        assert approx(pos.delta, abs=self.delta) == exp.delta

    def test_mercury(self, sphera):
        self._compare(
            PlanetId.MERCURY,
            EclipticPosition(lmbda=275.88530, beta=1.47425, delta=0.98587),
            sphera,
        )

    def test_venus(self, sphera):
        self._compare(
            PlanetId.VENUS,
            EclipticPosition(lmbda=264.15699, beta=1.42582, delta=1.22905),
            sphera,
        )

    def test_mars(self, sphera):
        self._compare(
            PlanetId.MARS,
            EclipticPosition(lmbda=214.98173, beta=1.67762, delta=1.41366),
            sphera,
        )

    def test_jupiter(self, sphera):
        self._compare(
            PlanetId.JUPITER,
            EclipticPosition(lmbda=270.30024, beta=0.29758, delta=6.10966),
            sphera,
        )

    def test_saturn(self, sphera):
        self._compare(
            PlanetId.SATURN,
            EclipticPosition(lmbda=225.37862, beta=2.33550, delta=10.04942),
            sphera,
        )

    def test_uranus(self, sphera):
        self._compare(
            PlanetId.URANUS,
            EclipticPosition(lmbda=252.17354, beta=0.05160, delta=19.63393),
            sphera,
        )

    def test_neptune(self, sphera):
        self._compare(
            PlanetId.NEPTUNE,
            EclipticPosition(lmbda=270.07638, beta=1.16314, delta=31.11160),
            sphera,
        )

    def test_pluto(self, sphera):
        self._compare(
            PlanetId.PLUTO,
            EclipticPosition(lmbda=212.07989, beta=16.88244, delta=29.86118),
            sphera,
        )
