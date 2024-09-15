"""Converting between civil and Julian dates.

*Julian date* (JD) is the number of days elapsed since mean **UT noon of
January 1st 4713 BC**. This system of time measurement is widely adopted by
the astronomers.

For better precision around the XX century, we use the epoch
**1900 January 0.5** (1989 December 31.5) as the starting point.
(See *"Astronomy With Your Personal Computer"*, p.14.) This kind of Julian
date is referred as 'DJD'. To convert DJD to JD and vise versa, use
`DJD_TO_JD` constant:

    jd = djd + DJD_TO_JD.


The module contains some other usefull calendar-related functions, such as
`weekDay`, `dayOfYear`, `isLeapYear`.

"""

from collections import namedtuple
from datetime import datetime, timezone
from math import floor, trunc, modf, fabs

from astropc.mathutils import dms

__author__ = "ilbagatto"
__license__ = "MIT"
__version__ = "0.0.1"


YearMonthDay = namedtuple("YearMonthDay", "year month day")
YearMonthDay.__doc__ = """Calendar date is expressed as year, month (1 - 12) and day
with hours as fractional part.
"""

DEFAULT_GREG_START = YearMonthDay(1582, 10, 15)
"""Year, month and day when Gregorian calendar was introduced."""


DJD_TO_JD = 2415020
"""Difference in days between DJD and the Standard Juian Date."""

DAYS_PER_CENT = 36525
"""Days per century (36525)."""


class CalendarException(Exception):
    """Base class for calendar related exceptions."""

    def __init__(self, message: str):
        self.message = message

    def __str__(self) -> str:
        return self.message


def after_gregorian(
    year: int,
    month: int,
    day: float | int,
    greg_start: YearMonthDay = DEFAULT_GREG_START,
) -> bool:
    """Does a given date falls to period after introducion of Gregorian calendar?

    Args:
        year: civil year
        month: month (1 - 12)
        day: a day with hours as fractional part
        greg_start: a date when Gregorian calendar was introduced, Oct. 15, 1582 by default.

    Returns:
        True for Gregorian dates.
    """
    if year < greg_start.year:
        return False
    if year > greg_start.year:
        return True
    if month < greg_start.month:
        return False
    if month > greg_start.month:
        return True
    return bool(day >= greg_start.day)


def jul_day(
    year: int,
    month: int,
    day: float,
    gregorian_start: YearMonthDay = DEFAULT_GREG_START,
) -> float:
    """Converts calendar date into Julian days elapsed since **1900, Jan 0.5** (1899 Dec 31.5).

    Args:

        year: civil year, negative for BC era.
        month: month (1 - 12).
        day: a day with hours as fractional part.
        greg_start: a date when Gregorian calendar was introduced, Oct. 15, 1582 by default.

    Raises:
        `CalendarException` if year is zero (astronomical years are not allowed).


    >>> jul_day(1984, 8, 29.0) # August 29, 1984, 0h UTC
    30921.5

    """
    if year == 0:
        raise CalendarException("Zero year not allowed!")

    # convert civil year to astronomic year
    y = year + 1 if year < 0 else year

    m = month
    if month < 3:
        m += 12
        y -= 1

    if after_gregorian(year, month, day, gregorian_start):
        # after Gregorian calendar
        a = trunc(y / 100)
        b = 2 - a + trunc(a / 4)
    else:
        b = 0

    f = 365.25 * y
    c = trunc(f - 0.75 if y < 0 else f) - 694025
    e = trunc(30.6001 * (m + 1))
    return b + c + e + day - 0.5


def cal_date(djd: float) -> tuple[int, int, float]:
    """Converts Julian day (DJD) into the calendar date.

    Args:
        djd: number of Julian days since 1900 Jan. 0.5

    Returns:
        year: civil year
        month: month (1 - 12)
        day: a day with hours as fractional part

    >>> cal_date(30921.5)
    (1984, 8, 29.0)
    """
    d = djd + 0.5
    (f, i) = modf(d)

    if i > -115860:
        a = floor(i / 36524.25 + 9.9835726e-1) + 14
        i += 1 + a - floor(a / 4)

    b = floor(i / 365.25 + 8.02601e-1)
    c = i - floor(365.25 * b + 7.50001e-1) + 416
    g = floor(c / 30.6001)
    da = c - floor(30.6001 * g) + f
    mo = g - (13 if g > 13.5 else 1)
    ye = b + (1900 if mo < 2.5 else 1899)
    # convert astronomical, zero-based year to civil
    if ye < 1:
        ye -= 1

    return ye, mo, da


def djd_midnight(djd: float) -> float:
    """Juian day at Greenwich midnight.

    Args:
        djd: number of Julian days since 1900 Jan. 0.5

    Returns:
        Julian day (DJD) at Greenwich midnight.

    >>> djd_midnight(23772.99)
    23772.5
    """
    f = floor(djd)
    return f + (0.5 if fabs(djd - f) >= 0.5 else -0.5)


def weekday(djd: float) -> int:
    """Day of week.

    Args:
        djd: number of Julian days elapsed since 1900, Jan 0.5.

    Returns:
        number in range (0..6) corresponding to weekDay: 0 for Sunday, 1 for Monday and so on.

    >>> weekday(30921.5)
    3 # Wednesday

    """
    d0 = djd_midnight(djd)
    j0 = d0 + DJD_TO_JD
    return int((j0 + 1.5) % 7)


def is_leapyear(ye: int) -> bool:
    """Is given year a leap-year?

    Args:
        year: civil year

    Returns:
        True for leap-year, else False.

    >>> leap_year(2000)
    True

    >>> leap_year(2001)
    False
    """
    return (ye % 4 == 0) and ((ye % 100 != 0) or (ye % 400 == 0))


def day_of_year(year: int, month: int, day: float) -> int:
    """Number of days in the year up to a particular Gregorian date.

    Args:
        year: civil year
        month: month (1 - 12)
        day: a day with hours as fractional part

    Returns:
        number of days, including the given date (1 .. 366)

    >>> day_of_year(2000, 4, 1)
    92

    >>> day_of_year(2000, 1, 1)
    1

    """
    k = 1 if is_leapyear(year) else 2
    a = floor(275 * month / 9)
    b = floor(k * ((month + 9) / 12.0))
    c = floor(day)
    return a - b + c - 30


def djd_zero(year: int) -> float:
    """
    DJD corresponding to January 0.0 of a given year.

    Zero day is a special case of date: it indicates 12h UT of previous calendar
    date. For instance, 1900 January 0.5 is often used instead ofn899 December 31.5
    to designate start of the astronomical epoch.

    Args:
        year: civil year

    Returns:
        Julian day (DJD)

    >>> djd_zero(2010)
    40176.5
    """
    y = year - 1
    a = trunc(y / 100)
    return trunc(365.25 * y) - a + trunc(a / 4) - 693595.5


def djd_to_datetime(djd: float) -> datetime:
    """Convert djd, a number of Julian days elapsed since 1900, Jan 0.5,
    to Python datetime.datetime object in UTC.

    Args:
        djd: number of Julian days elapsed since 1900, Jan 0.5

    Returns:
        datetime.datetime object.

    """
    (year, month, day) = cal_date(djd)
    f, d = modf(day)
    (ho, mi, se) = dms(f * 24)
    f, s = modf(se)
    ms = f * 1000000  # microsecond
    return datetime(year, month, int(d), ho, mi, int(s), int(ms), tzinfo=timezone.utc)
