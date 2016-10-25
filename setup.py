#!/usr/bin/env python
import os
import re
import sys
import warnings

from setuptools import setup, find_packages
from setuptools import Command

MAJOR = 2
MINOR = 1
ISRELEASED = False
VERSION = '%d.%d' % (MAJOR, MINOR)
QUALIFIER = ''


DISTNAME = 'rbm'
LICENSE = 'Apache'
AUTHOR = 'RBM Developers'
AUTHOR_EMAIL = 'yifanc6@uw.edu'
URL = 'https://github.com/UW-Hydro/RBM'


INSTALL_REQUIRES = ['Fortran 90']
TESTS_REQUIRE = ['pytest >= 2.7.1']

DESCRIPTION = "River Basin Model"
LONG_DESCRIPTION = """
RBM is a one-dimensional stream-temperature model that solves
the time-dependent equation for the transfer of thermal energy
across the air-water interface (Yearsley, 2009, 2012).
The model uses a semi-Lagrangian, particle tracking numerical
scheme that is highly scalable in space and time. The model
software was written initially for purposes of developing the
Total Maximum Daily Load (TMDL) for water temperature in the
Columbia River system (Yearsley, 2003). The model had been modified
to make use of the output from the large-scale hydrologic model VIC
and regional scale model DHSVM.
"""

# code to extract and write the version copied from pandas
FULLVERSION = VERSION
write_version = True

setup(name=DISTNAME,
      version=FULLVERSION,
      license=LICENSE,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      install_requires=INSTALL_REQUIRES,
      tests_require=TESTS_REQUIRE,
      url=URL,
      packages=find_packages()})
