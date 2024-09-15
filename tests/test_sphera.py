from math import radians

from pytest import approx, fixture

from astropc.planets.ids import PlanetId
from astropc.planets.sphera import CelestialSphera

DELTA = 1e-6


@fixture()
def sphera():
    return CelestialSphera.create(30700.5)


class TestSpheraProps:

    def test_t(self, sphera):
        assert approx(sphera.t, abs=DELTA) == 0.8405338809034908

    def test_sungeo_phi(self, sphera):
        assert approx(sphera.sun_geo.phi, abs=DELTA) == 300.1307723107521

    def test_sungeo_rho(self, sphera):
        assert approx(sphera.sun_geo.rho, abs=DELTA) == 0.9839698373786032

    def test_sun_manom(self, sphera):
        assert approx(sphera.manom_sun, abs=DELTA) == radians(16.89671827974547)

    def test_dpsi(self, sphera):
        assert approx(sphera.nutation.dpsi, abs=DELTA) == -0.004176852920102668

    def test_deps(self, sphera):
        assert approx(sphera.nutation.deps, abs=DELTA) == 0.0006849657311651972

    def test_obliquity(self, sphera):
        assert approx(sphera.obliquity, abs=DELTA) == 23.442041099302447


class TestMeanAnomalies:
    def test_manom_mercury(self, sphera):
        ma = sphera.get_mean_anomaly(PlanetId.MERCURY)
        assert approx(ma, abs=DELTA) == 1.7277480419370512

    def test_manom_venus(self, sphera):
        ma = sphera.get_mean_anomaly(PlanetId.VENUS)
        assert approx(ma, abs=DELTA) == 1.3753354318768864

    def test_manom_mars(self, sphera):
        ma = sphera.get_mean_anomaly(PlanetId.MARS)
        assert approx(ma, abs=DELTA) == 3.616595436914378

    def test_manom_jupiter(self, sphera):
        ma = sphera.get_mean_anomaly(PlanetId.JUPITER)
        assert approx(ma, abs=DELTA) == 4.469600429159891

    def test_manom_saturn(self, sphera):
        ma = sphera.get_mean_anomaly(PlanetId.SATURN)
        assert approx(ma, abs=DELTA) == 2.133162332104278

    def test_manom_uranus(self, sphera):
        ma = sphera.get_mean_anomaly(PlanetId.URANUS)
        assert approx(ma, abs=DELTA) == 1.2691334849854374

    def test_manom_neptune(self, sphera):
        ma = sphera.get_mean_anomaly(PlanetId.NEPTUNE)
        assert approx(ma, abs=DELTA) == 3.863368991888066

    def test_manom_pluto(self, sphera):
        ma = sphera.get_mean_anomaly(PlanetId.PLUTO)
        assert approx(ma, abs=DELTA) == 6.138953042601936


class TestCache:
    def test_get_same_orbit_instance(self, sphera, mocker):
        spy = mocker.spy(sphera, "_instantiate_orbit")
        sphera.get_orbit_instance(PlanetId.MERCURY)
        sphera.get_orbit_instance(PlanetId.MERCURY)
        spy.assert_called_once_with(PlanetId.MERCURY)

    def test_get_different_orbit_instance(self, sphera, mocker):
        spy = mocker.spy(sphera, "_instantiate_orbit")
        sphera.get_orbit_instance(PlanetId.MERCURY)
        sphera.get_orbit_instance(PlanetId.JUPITER)
        assert spy.call_count == 2
