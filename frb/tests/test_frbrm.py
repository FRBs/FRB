# Module to run tests on surveys
#  Most of these are *not* done with Travis yet
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units

from frb import rm


def test_galacticrm():
    repeater_coord = SkyCoord('05h31m58.698s +33d8m52.59s', frame='icrs')
    RM, RM_err = rm.galactic_rm(repeater_coord)
    # Test
    assert RM.unit == units.rad/units.m**2
    assert np.isclose(RM.value, -17.727994918823242)

