"""Time related routines

## Civil vs. Astronomical year

There is disagreement between astronomers and historians about how to count
the years preceding the year 1. Astronomers generally use zero-based system.
The year before the year +1, is the *year zero*, and the year preceding the
latter is the *year -1*. The year which the historians call 585 B.C. is
actually the year -584.

In this module all subroutines accepting year ([isLeapYear], [cal2djd] etc.)
assume that **there is no year zero**. Conversion from the civil to the
astronomical time scale is done internally. Thus, the sequence of years is:
`BC -3, -2, -1, 1, 2, 3, AD`.

## Time

Time is represented by fractional part of a day. For example, 7h30m UT
is `(7 + 30 / 60) / 24 = 0.3125`.

## Zero day

Zero day is a special case of date: it indicates 12h UT of previous calendar
date. For instance, *1900 January 0.5* is often used instead of
*1899 December 31.5* to designate start of the astronomical epoch.

##  Gregorian calendar

_Civil calendar_ in most cases means _proleptic Gregorian calendar_. it is
assumed that Gregorian calendar started at *Oct. 4, 1582*, when it was first
adopted in several European countries. Many other countries still used the
older Julian calendar. In Soviet Russia, for instance, Gregorian system was
accepted on **Jan 26, 1918**.

See
[Wiki article](https://en.wikipedia.org/wiki/Gregorian_calendar#Adoption_of_the_Gregorian_Calendar)
"""

from .solequ import *

# flake8: noqa F401, F403
from .sun import *
