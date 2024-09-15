from pytest import approx, mark

from astropc.coco import ecl2equ, equ2ecl, equ2hor, hor2equ

DELTA = 1e-4


class TestEquEcl:

    @mark.parametrize(
        "alpha, delta, eps, lmbda",
        [
            (216.7375, 32.351389, 23.445745, 200.318517),
            (0.022917, -87.203333, 23.438719, 277.001739),
        ],
    )
    def test_lambda(self, alpha: float, delta: float, eps: float, lmbda: float):
        got_lmbda, _ = equ2ecl(alpha, delta, eps)
        assert approx(got_lmbda, abs=DELTA) == lmbda

    @mark.parametrize(
        "alpha, delta, eps, beta",
        [
            (216.7375, 32.351389, 23.445745, 43.787175),
            (0.022917, -87.203333, 23.438719, -66.405486),
        ],
    )
    def test_beta(self, alpha: float, delta: float, eps: float, beta: float):
        _, got_beta = equ2ecl(alpha, delta, eps)
        assert approx(got_beta, abs=DELTA) == beta


class TestEclEqu:
    @mark.parametrize(
        "lmbda, beta, eps, alpha",
        [
            (200.318517, 43.787175, 23.445745, 216.7375),
            (277.001739, -66.405486, 23.438719, 0.022917),
        ],
    )
    def test_alpha(self, lmbda, beta, eps, alpha):
        got_alpha, _ = ecl2equ(lmbda, beta, eps)
        assert approx(got_alpha, abs=DELTA) == alpha

    @mark.parametrize(
        "lmbda, beta, eps, delta",
        [
            (200.318517, 43.787175, 23.445745, 32.351389),
            (277.0017389, -66.405486, 23.438719, -87.203333),
        ],
    )
    def test_delta(self, lmbda, beta, eps, delta):
        _, got_delta = ecl2equ(lmbda, beta, eps)
        assert approx(got_delta, abs=DELTA) == delta


class TestEquHor:
    @mark.parametrize(
        "ha, delta, theta, azimuth",
        [
            (8.622222, 14.398611, 51.25, 310.259333),
            (23.316667, -43.0, -20.520278, 161.388611),
        ],
    )
    def test_azimuth(self, ha, delta, theta, azimuth):
        got_azimuth, _ = equ2hor(ha * 15, delta, theta)
        assert approx(got_azimuth, abs=DELTA) == azimuth

    @mark.parametrize(
        "ha, delta, theta, altitude",
        [
            (8.622222, 14.398611, 51.25, -10.972444),
            (23.316667, -43.0, -20.520278, 65.935028),
        ],
    )
    def test_altitude(self, ha, delta, theta, altitude):
        _, got_altitude = equ2hor(ha * 15, delta, theta)
        assert approx(got_altitude, abs=DELTA) == altitude


class TestHorEqu:
    @mark.parametrize(
        "azimuth, altitude, theta, ha",
        [
            (310.259333, -10.972444, 51.25, 8.622222),
            (161.388611, 65.935028, -20.520278, 23.316667),
        ],
    )
    def test_hour_angle(self, azimuth, altitude, theta, ha):
        got_ha, _ = hor2equ(azimuth, altitude, theta)
        assert approx(got_ha / 15, abs=DELTA) == ha

    @mark.parametrize(
        "azimuth, altitude, theta, delta",
        [
            (310.259333, -10.972444, 51.25, 14.398611),
            (161.388611, 65.935028, -20.520278, -43.0),
        ],
    )
    def test_delta(self, azimuth, altitude, theta, delta):
        _, got_delta = hor2equ(azimuth, altitude, theta)
        assert approx(got_delta, abs=DELTA) == delta
