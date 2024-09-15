import pytest

from astropc.timeutils.deltat import delta_t


@pytest.mark.parametrize(
    "djd, expected, msg",
    [
        (-102146.5, 119.51, "historical start"),  # 1620-05-01
        (-346701.5, 1820.325, "after 948"),  # 950-10-01
        (44020.5, 93.81, "after 2010"),  # 2020-07-10
        (109582.5, 407.2, "after 2100"),  # ?
    ],
)
def test_deltat(djd, expected, msg):
    assert pytest.approx(delta_t(djd), abs=0.1) == expected
