from pytest import approx

from astropc.sun.solequ import SolEquType, sol_equ
from astropc.sun.sun import apparent, true_geocentric

delta = 1e-4  # result precision


class TestSunPosition:

    cases = (
        {
            "djd": 30916.5,  # 24 Aug 1984 00:00
            "l": 151.01309547440778,
            "r": 1.010993800005251,
            "ap": 151.0035132296576,
        },
        {
            "djd": 30819.10833333333,  # 18 May 1984 14:36
            "l": 57.83143688493146,
            "r": 1.011718488789592,
            "ap": 57.82109236581925,
        },
        {
            "djd": 28804.5,  # 12 Nov 1978 00:00
            "l": 229.2517039627867,
            "r": 0.9898375,
            "ap": 229.2450957063683,
        },
        {
            "djd": 33888.5,  # 1992, Oct. 13 0h
            "l": 199.90600618015975,
            "r": 0.9975999344847888,
            "ap": 199.9047664927989,  # Meeus: 199.90734722222223
        },
    )

    def test_geometric(self):
        for c in self.cases:
            t = c["djd"] / 36525
            geo = true_geocentric(t)
            assert approx(geo.phi, delta) == c["l"]
            assert approx(geo.rho, delta) == c["r"]

    def test_apparent(self):
        for c in self.cases:
            geo = apparent(c["djd"], ignore_light_travel=True)
            assert approx(geo.phi, delta) == c["ap"]


class TestSolEqu:
    cases = (
        # from 'Astronomical algorithms' by Jean Meeus, p.168
        {
            "djd": 22817.39,
            "year": 1962,
            "event": SolEquType.JUNE_SOLSTICE,
            "angle": 90.0,
            "title": 'June solstice (Meeus, "Astronomical Algorithms")',
        },
        # from http://www.usno.navy.mil/USNO/astronomical-applications/data-services/earth-seasons
        {
            "djd": 36603.815972,
            "year": 2000,
            "event": SolEquType.MARCH_EQUINOX,
            "angle": 0.0,
            "title": "March equinox",
        },
        {
            "djd": 36696.575000,
            "year": 2000,
            "event": SolEquType.JUNE_SOLSTICE,
            "angle": 90.0,
            "title": "June solstice",
        },
        {
            "djd": 36790.227778,
            "year": 2000,
            "event": SolEquType.SEPTEMBER_EQUINOX,
            "angle": 180.0,
            "title": "September equinox",
        },
        {
            "djd": 36880.067361,
            "year": 2000,
            "event": SolEquType.DECEMBER_SOLSTICE,
            "angle": 270.0,
            "title": "December solstice",
        },
    )

    def test_djd(self):
        for c in self.cases:
            evt = sol_equ(c["year"], c["event"])
            assert approx(evt.djd, abs=1e-2) == c["djd"]

    def test_sun(self):
        for c in self.cases:
            evt = sol_equ(c["year"], c["event"])
            assert approx(evt.sun, abs=delta) == c["angle"]
