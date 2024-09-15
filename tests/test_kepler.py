from pytest import approx, mark
from astropc.kepler import eccentric_anomaly, true_anomaly

DELTA = 1e-4  # result precision


@mark.parametrize(
    "m, s, e",
    [
        (3.5208387374141448, 0.016718, 3.5147440476661806),
        (0.763009079752865, 0.965, 1.7176273861066755),
    ],
)
def test_eccentric_anomaly(m, s, e):
    assert approx(eccentric_anomaly(s, m), DELTA) == e


@mark.parametrize(
    "s, e, ta",
    [
        (0.016718, 3.5147440476661806, -2.774497552017826),
        (0.965, 1.7176273861066755, 2.9122563898777387),
    ],
)
def test_true_anomaly(s, e, ta):
    assert approx(true_anomaly(s, e), DELTA) == ta
