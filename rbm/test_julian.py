import leapyear
import pytest
import datetime
import calendar
import numpy as np
import pandas as pd

def test_leap_year():
    for year in np.arange(1900, 2100):
        actual = calendar.isleap(year)
        if actual == False:
            actual = 0
        else:
            actual = 1
        assert not bool(vic_lib.leap_year(year, calendars['noleap']))
        assert actual == leapyear.leapyear(year)
