# AstroPC

Library of core routines for practical astronomy.

- [AstroPC](#astropc)
  - [Features](#features)
  - [Installation](#installation)
    - [Virtual environment](#virtual-environment)
    - [Install the package](#install-the-package)
  - [Usage](#usage)
  - [Unit tests](#unit-tests)
  - [Additional information](#additional-information)
    - [Civil vs. Astronomical year](#civil-vs-astronomical-year)
    - [Zero day](#zero-day)
    - [Gregorian calendar](#gregorian-calendar)
  - [Sources](#sources)
  - [How to contribute](#how-to-contribute)
  - [Footnotes](#footnotes)


## Features

* Miscallenous mathematical routines usefull for practical astronomy.
* Convertions between *Civil* and *Julian* dates.
* Calculates difference between *Universial Coordinated* and *Terrestrial Dynamic* time (*Delta-T*).
* Calculates *Sidereal* (*Stellar*) time.
* Calculates *Nutation* and *Obliquity of the ecliptic*.
* Transforming between various types of celestial coordinates.
* Accurate positions of Sun, Moon and the planets, including _Pluto_.
* Time of solstices, equinoxes and lunations.


## Installation

### Virtual environment

Create virtual environment for Python3.12 or later and activate it.

On Linux:

```console
$ python3.12 -m venv .venv
$ . ./venv/bin/activate
```

For details see [Create and Use Virtual Environments](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#create-and-use-virtual-environments).


### Install the package

```console
$ pip install .
```

Or, for development mode:

```console
$ pip install -e '.[dev]'
```

## Usage

See tests/ and examples/ for usage examples.

## Unit tests

From the project root directory:

```console
$ pytest tests
```

## Additional information

### Civil vs. Astronomical year

There is disagreement between astronomers and historians about how to count
the years preceding the year 1. Astronomers generally use zero-based system.
The year before the year +1, is the *year zero*, and the year preceding the
latter is the *year -1*. The year which the historians call *585 B.C.* is
actually the year *-584*.

### Zero day

Zero day is a special case of date: it indicates 12h UT of previous calendar
date. For instance, *1900 January 0.5* is often used instead of
*1899 December 31.5* to designate start of the astronomical epoch.

###  Gregorian calendar

_Civil calendar_ in most cases means _proleptic Gregorian calendar_. it is
assumed that Gregorian calendar started at *Oct. 4, 1582*, when it was first
adopted in several European countries. Many other countries still used the
older Julian calendar. In Soviet Russia, for instance, Gregorian system was
accepted on **Jan 26, 1918**. See
[Wiki article](https://en.wikipedia.org/wiki/Gregorian_calendar#Adoption_of_the_Gregorian_Calendar).

## Sources

The formulae were adopted from the following sources:

* _Peter Duffett-Smith, "Astronomy With Your Personal Computer", Cambridge University Press, 1997_
* _Jean Meeus, "Astronomical Algorithms", 2d edition, Willmann-Bell, 1998_
* _J.L.Lawrence, "Celestial Calculations", The MIT Press, 2018_


## How to contribute

You may contribute to the project by many different ways, starting from refining and correcting its documentation,
especially if you are a native English speaker, and ending with improving the code base. Any kind of testing and
suggestions are welcome.

You may follow the standard Github procedures or, in case you are not comfortable with them, just send your suggestions
to the author by other means.



## Footnotes

[^1]: *Julian Period* 7980 years of numbered days to be used in determining time elapsed between various historical events.
It is the product of three calendar cycles: `28 (solar cycle) × 19 (lunar cycle) × 15 (indiction cycle) = 7980 years`

~~~~
