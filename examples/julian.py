# flake8: noqa E402
from datetime import datetime
from pathlib import Path
import sys

libdir = Path(__file__).parent.parent.as_posix()
sys.path.insert(0, libdir)

from astropc.mathutils import ddd
from astropc.timeutils.julian import DJD_TO_JD, jul_day

USAGE = """Usage: python3 julian.py [-h, --help | DATETIME]

DATETIME is a civil date and time in ISO-8601 format,
YYYY-MM-DDT[HH[:MM[:SS[.mmm[uuu]]]]][+HH:MM] like:
"%s".
Current date/time by default.
"""

def main():
    if len(sys.argv) > 1:
        if sys.argv[1] in ("-h", "--help"):
            print(USAGE % datetime.now())
            sys.exit(0)

        try:
            dt = datetime.fromisoformat(sys.argv[1])
        except ValueError:
            print(f"Unexpected date/time: '{sys.argv[1]}'")
            sys.exit(1)
    else:
        dt = datetime.now()

    hm = ddd(dt.hour, dt.minute, dt.second)  # hour, minute and seconds as decimal hours
    djd = jul_day(dt.year, dt.month, dt.day + hm / 24)
    jd = djd + DJD_TO_JD
    print("%12.6f" % jd)
    
if __name__ == '__main__':
    main()
