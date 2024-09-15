from pytest import approx, fixture
from astropc.planets.orbit import MLTerms, OElements, Terms


DELTA = 1e-6
T = 0.8405338809034908


class TestTerms:

    def test_assemble_standard_terms(self):
        terms = Terms(75.899697, 1.5554889, 0.0002947)
        assert approx(terms.assemble(T), abs=DELTA) == 77.2073463265456

    def test_assemble_ml_terms(self):
        terms = MLTerms(178.179078, 415.2057519, 0.0003011)
        assert approx(terms.assemble(T), abs=DELTA) == 176.2000171915306


class TestOElements:

    @fixture()
    def orbit(self):
        oe = OElements(
            mean_longitude=MLTerms(178.179078, 415.2057519, 3.011e-4),
            perihelion=Terms(75.899697, 1.5554889, 2.947e-4),
            eccentricity=Terms(2.0561421e-1, 2.046e-5, -3e-8),
            inclination=Terms(7.002881, 1.8608e-3, -1.83e-5),
            mean_node=Terms(47.145944, 1.1852083, 1.739e-4),
            major_semiaxis=3.870986e-1,
        )
        return oe.instantiate(T)

    def test_mean_anomaly(self, orbit):
        assert approx(orbit.mean_anomaly, abs=DELTA) == 1.7277480419370512

    def test_daily_motion(self, orbit):
        assert approx(orbit.daily_motion, abs=DELTA) == 0.07142545459475612

    def test_perihlion(self, orbit):
        assert approx(orbit.perihelion, abs=DELTA) == 1.34752240012577

    def test_eccentricity(self, orbit):
        assert approx(orbit.eccentricity, abs=DELTA) == 0.20563138612828713

    def test_mean_node(self, orbit):
        assert approx(orbit.mean_node, abs=DELTA) == 0.8402412010285969

    def test_inclination(self, orbit):
        assert approx(orbit.inclination, abs=DELTA) == 0.12225040301524157

    def test_semiaxis(self, orbit):
        assert approx(orbit.major_semiaxis, abs=DELTA) == 0.3870986
