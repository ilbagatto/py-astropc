from pytest import approx, fixture, mark, raises
from astropc.mathutils import ddd
from astropc.timeutils import CalendarException
from astropc.timeutils.julian import (
    cal_date,
    day_of_year,
    djd_midnight,
    djd_to_datetime,
    djd_zero,
    is_leapyear,
    jul_day,
    weekday,
)


class TestJulianDay:
    delta = 1e-6

    @mark.parametrize(
        "year, month, day, djd",
        [
            (1984, 8, 29.0, 30921.5),
            (1899, 12, 31.5, 0.0),
            (1938, 8, 17.0, 14107.5),
            (1, 1, 1.0, -693596.5),
            (-4713, 7, 12.0, -2414827.5),
            (-4713, 1, 1.5, -2415020.0),
        ],
    )
    def test_civil_to_julian(self, year, month, day, djd):
        got = jul_day(year, month, day)
        assert approx(got, abs=self.delta) == djd

    @mark.parametrize(
        "year, month, day, djd",
        [
            (1984, 8, 29.0, 30921.5),
            (1899, 12, 31.5, 0.0),
            (1938, 8, 17.0, 14107.5),
            (1, 1, 1.0, -693596.5),
            (-4713, 7, 12.0, -2414827.5),
            (-4713, 1, 1.5, -2415020.0),
        ],
    )
    def test_julian_to_civil(self, year, month, day, djd):
        (got_year, got_month, got_day) = cal_date(djd)
        assert got_year == year
        assert got_month == month
        assert approx(got_day, abs=self.delta) == day

    def test_zero_day(self):
        got = jul_day(1900, 1, 0.5)
        assert approx(got, abs=self.delta) == 0.0

    def test_zero_year(self):
        with raises(CalendarException):
            jul_day(0, 12, 1)

    def test_djd_zero(self):
        got = djd_zero(2010)
        assert approx(got, abs=self.delta) == 40176.5


class TestJDMidnight:
    delta = 1e-6

    def before_noon(self):
        assert djd_midnight(23772.99) == approx(23772.5, abs=self.delta)

    def after_noon(self):
        assert djd_midnight(23773.3) == approx(23772.5, abs=self.delta)

    def prev_day_before_midnight(self):
        assert djd_midnight(23772.4) == approx(23771.5, abs=self.delta)

    def prev_day_before_noon(self):
        assert djd_midnight(23771.9) == approx(23771.5, abs=self.delta)

    def next_day_after_midnight(self):
        assert djd_midnight(23773.6) == approx(23773.5, abs=self.delta)


@mark.parametrize(
    "djd, wd",
    [
        (30921.5, 3),
        (0.0, 0),
        (14107.5, 3),
        (-693596.5, 6),
        # Not sure about weekDays of the next two dates; there are controverses;
        # Perl  DateTime module gives weekDays 5 and 4 respectively
        (-2414827.5, 5),
        (-2415020.0, 1),
        (23772.99, 1),
    ],
)
def test_weekdays(djd, wd):
    assert weekday(djd) == wd


@mark.parametrize(
    "year",
    [
        2000,
        2004,
        2008,
        2012,
        2016,
        2020,
        2024,
        2028,
        2032,
        2036,
        2040,
        2044,
        2048,
    ],
)
def test_leap_years(year):
    assert is_leapyear(year)


@mark.parametrize(
    "year",
    [
        2001,
        2003,
        2010,
        2014,
        2017,
        2019,
        2025,
        2026,
        2035,
        2038,
        2045,
        2047,
        2049,
    ],
)
def test_non_leap_years(year):
    assert not is_leapyear(year)


def test_non_leap_day_of_year():
    assert day_of_year(1990, 4, 1) == 91


def test_leap_year_of_year():
    assert day_of_year(2000, 4, 1) == 92


def test_first_day_of_year():
    assert day_of_year(2000, 1, 1) == 1


class TestDateTime:

    @fixture()
    def dt(self):
        djd = jul_day(1965, 2, 1 + ddd(11, 46) / 24)
        return djd_to_datetime(djd)

    def test_year(self, dt):
        assert dt.year == 1965

    def test_month(self, dt):
        assert dt.month == 2

    def test_day(self, dt):
        assert dt.day == 1

    def test_hour(self, dt):
        assert dt.hour == 11

    def test_minute(self, dt):
        assert dt.minute == 46

    def test_second(self, dt):
        assert dt.second == 0
