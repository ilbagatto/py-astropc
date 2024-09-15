import pytest

from astropc.timeutils.sidereal import djd_to_sidereal, sidereal_to_utc

DELTA = 1e-4

cases = (
    {"djd": 30923.851053, "lst": 7.072111, "utc": 8.425278, "ok": True},  # 1984-08-31.4
    {"djd": 683.498611, "lst": 3.525306, "utc": 23.966667, "ok": False},  # 1901-11-15.0
    {"djd": 682.501389, "lst": 3.526444, "utc": 0.033333, "ok": False},  # 1901-11-14.0
    {
        "djd": 29332.108931,
        "lst": 4.668119,
        "utc": 14.614353,
        "ok": True,
    },  # 1980-04-22.6
)


def test_utc_to_gst():
    for c in cases:
        got = djd_to_sidereal(c["djd"])
        assert pytest.approx(got, abs=DELTA) == c["lst"]


def test_gst_to_utc_non_ambiguous():
    for c in [c for c in cases if c["ok"]]:
        (utc, not_ok) = sidereal_to_utc(c["lst"], c["djd"])
        assert pytest.approx(utc, abs=DELTA) == c["utc"]
        assert not not_ok


def test_gst_to_utc_ambiguous():
    for c in [c for c in cases if not c["ok"]]:
        (utc, not_ok) = sidereal_to_utc(c["lst"], c["djd"])
        # assert pytest.approx(utc, abs=DELTA) == c["utc"]
        assert not_ok
