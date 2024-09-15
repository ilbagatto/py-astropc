from pytest import approx, mark
from astropc.moon.lunations import Quarter, find_closest_phase


@mark.parametrize(
    "year, month, day, djd",
    [
        (1984, 9, 1, 30919.3097),  # 1984, 8, 26, 19, 26
        (1968, 12, 12, 25190.263194),  # 1968, 12, 19, 18, 19
        (2019, 8, 21, 43705.94287),
    ],
)
def test_new_moon(year, month, day, djd):
    assert approx(djd) == find_closest_phase(Quarter.NEW_MOON, year, month, day)


@mark.parametrize(
    "year, month, day, djd",
    [
        (2019, 8, 21, 43712.63302),  # 1984, 8, 26, 19, 26
    ],
)
def test_first_quarter(year, month, day, djd):
    assert approx(djd) == find_closest_phase(Quarter.FIRST_QUARTER, year, month, day)


@mark.parametrize(
    "year, month, day, djd",
    [
        (1984, 9, 1, 30933.79236),  # 1984, 9, 10, 7, 1
        (1965, 2, 1, 23787.52007),
        (2019, 8, 21, 43720.69049),
    ],
)
def test_full_moon(year, month, day, djd):
    assert approx(djd) == find_closest_phase(Quarter.FULL_MOON, year, month, day)


@mark.parametrize(
    "year, month, day, djd",
    [(2044, 1, 1, 52616.49186), (2019, 8, 21, 43728.61252)],  # 1984, 9, 10, 7, 1
)
def test_last_quarter(year, month, day, djd):
    assert approx(djd) == find_closest_phase(Quarter.LAST_QUARTER, year, month, day)
