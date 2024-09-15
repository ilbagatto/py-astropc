import math
from typing import Literal

from pytest import approx, mark
from astropc import mathutils as mu


@mark.parametrize(
    "a, b",
    [(-23456789.9, -0.9), (-10.7, -0.7), (0.0, 0.0), (10.7, 0.7), (23456789.9, 0.9)],
)
def test_frac(a: float, b: float):
    assert approx(mu.frac(a)) == b


@mark.parametrize(
    "a, b",
    [
        (65.08828950998063, 31.7842235930254),
        (65.08518145432647, 30.6653235575305),
        (870.1176191104912, 42.3428797768338),
        (862.7609301843507, 273.934866366267),
        (3.496869604337422, 178.873057561472),
    ],
)
def test_frac360(a: float, b: float):
    assert approx(mu.frac360(a)) == b


class TestPolynome:
    delta = 1e-6

    def test_short_polynome(self):
        x = mu.polynome(10.0, 1.0, 2.0, 3.0)
        assert x == approx(321.0, self.delta)

    def test_long_polynome(self):
        x = mu.polynome(
            -0.127296372347707,
            0.409092804222329,
            -0.0226937890431606,
            -7.51461205719781e-06,
            0.0096926375195824,
            -0.00024909726935408,
            -0.00121043431762618,
            -0.000189319742473274,
            3.4518734094999e-05,
            0.000135117572925228,
            2.80707121362421e-05,
            1.18779351871836e-05,
        )
        assert x == approx(0.411961500152426, self.delta)


class TestSexagesimal:
    delta = 1e-2

    @mark.parametrize(
        "ddd, deg, min, sec",
        [
            (-37.583333, -37, 34, 59.9999),
            (55.75, 55, 45, 0.0),
            (-0.166667, 0, -10, 0.0),
            (-0.002778, 0, 0, -10.0),
            (0.0029355, 0, 0, 10.5678),
            (0, 0, 0, 0.0),
        ],
    )
    def test_dms(
        self,
        ddd: float | Literal[0],
        deg: Literal[-37] | Literal[55] | Literal[0],
        min: Literal[34] | Literal[45] | Literal[-10] | Literal[0],
        sec: float,
    ):
        got = mu.dms(ddd)
        assert got[0] == deg
        assert got[1] == min
        assert approx(got[2], abs=self.delta) == sec

    @mark.parametrize(
        "ddd, deg, min, sec",
        [
            (-37.583333, -37, 35, 0.0),
            (55.75, 55, 45, 0.0),
            (-0.166667, 0, -10, 0.0),
            (-0.002778, 0, 0, -10.0),
            (0.0029355, 0, 0, 10.5678),
            (0, 0, 0, 0.0),
        ],
    )
    def test_ddd(
        self,
        ddd: float | Literal[0],
        deg: Literal[-37] | Literal[55] | Literal[0],
        min: Literal[35] | Literal[45] | Literal[-10] | Literal[0],
        sec: float,
    ):
        assert mu.ddd(deg, min, sec) == approx(ddd, abs=self.delta)

    def test_zdms(self):
        x = 320.25
        exp = (10, 20, 15, 0)
        got = mu.zdms(x)
        assert got == exp


class TestRanges:
    delta = 1e-6

    def test_reduce_deg(self):
        arg = (-700, 0, 345, 700, 360, 324070.45)
        exp = (20, 0, 345, 340, 0, 70.45)
        for a, b in zip(arg, exp):
            assert approx(mu.reduce_deg(a), abs=self.delta) == b

    def test_reduce_rad(self):
        arg = (12.89, -12.89, 0.0, 10.0, math.pi, mu.PI2)
        exp = (0.323629385640829, 5.95955592153876, 0.0, 3.71681469282041, math.pi, 0.0)
        for a, b in zip(arg, exp):
            assert approx(mu.reduce_rad(a)) == b


class TestArcs:
    delta = 1e-6

    @mark.parametrize(
        "a, b, x", [(10.0, 270.0, 100.0), (350.0, 20.0, 30.0), (10.0, 20.0, 10.0)]
    )
    def test_degrees(self, a: float, b: float, x: float):
        assert approx(mu.shortest_arc_deg(a, b), abs=self.delta) == x

    @mark.parametrize(
        "a, b, x",
        [
            (0.17453292519943295, 4.71238898038469, 1.7453292519943295),
            (0.17453292519943295, 4.71238898038469, 1.7453292519943295),
            (0.17453292519943295, 4.71238898038469, 1.7453292519943295),
        ],
    )
    def test_radians(self, a: float, b: float, x: float):
        assert approx(mu.shortest_arc_rad(a, b), abs=self.delta) == x

    @mark.parametrize(
        "a, b, expected",
        [
            (75, 10, -65),
            (10, 75, 65),
            (280, 30, 110),
            (30, 280, -110),
        ],
    )
    def test_diff_angle(self, a, b, expected):
        got = mu.diff_angle(a, b)
        assert approx(got, abs=self.delta) == expected
