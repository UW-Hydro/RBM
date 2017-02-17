import pytest
import datetime
import calendar
import numpy as np
import pandas as pd

def test_leap_year():
    from rbm.test_function import leapyear
    print('1')
    for year in np.arange(1900, 2100):
        actual = calendar.isleap(year)
        if actual == False:
            actual = 0
        else:
            actual = 1
            print(year)
        assert actual == leapyear.leap_year(year)
