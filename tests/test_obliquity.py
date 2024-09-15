from pytest import approx
from astropc.obliq import calc_obliquity

DELTA = 1e-4  # result precision


def test_mean_obliquity():
    # P.Duffett-Smith, "Astronomy With Your Personal Computer", p.54
    cases = (
        {"djd": 29120.5, "eps": 23.441917},  # 1979-09-24.0
        {"djd": 36524.5, "eps": 23.439278},  # 2000-01-01.0
    )

    for c in cases:
        got = calc_obliquity(c["djd"])
        assert approx(got, abs=DELTA) == c["eps"]


def test_true_obliquity():
    deps = 9.443 / 3600
    got = calc_obliquity(31875.5, deps)  # 1987-04-10.0
    assert approx(got, abs=DELTA) == 23.4435694
