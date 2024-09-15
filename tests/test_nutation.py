from pytest import approx
from astropc.nutation import calc_nutation
from astropc.timeutils.julian import jul_day


delta = 1e-4
# result precision

cases = (
    {
        "djd": -15804.5,  # 1856 Sept. 23
        "dpsi": -0.00127601021242336,
        "deps": 0.00256293723137559,
    },
    {
        "djd": 36524.5,  # 2000 Jan. 1
        "dpsi": -0.00387728730373955,
        "deps": -0.00159919822661103,
    },
    {
        "djd": 28805.69,  # 1978 Nov 17
        "dpsi": -9.195562346652888e-04,
        "deps": -2.635113483663831e-03,
    },
    {
        "djd": jul_day(1989, 2, 4),
        "dpsi": 0.0023055555555555555,
        "deps": 0.0022944444444444444,
    },
    {
        "djd": jul_day(2000, 1, 1.5),
        "dpsi": -0.003877777777777778,
        "deps": -0.0016,
    },
    {
        "djd": jul_day(1995, 4, 23),
        "dpsi": 0.0026472222222222223,
        "deps": -0.002013888888888889,
    },
    {
        "djd": jul_day(1965, 2, 1),
        "dpsi": -0.0042774118548615766,
        "deps": 0.000425,
    },
)


def test_dpsi():
    for c in cases:
        t = c["djd"] / 36525
        nut = calc_nutation(t)
        assert approx(nut.dpsi, abs=delta) == c["dpsi"]


def test_deps():
    for c in cases:
        t = c["djd"] / 36525
        nut = calc_nutation(t)
        assert approx(nut.deps, abs=delta) == c["deps"]
