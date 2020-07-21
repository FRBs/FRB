# Module to run tests on galaxy modules
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.table import Table
from astropy import units
from astropy.coordinates import SkyCoord

from frb.galaxies import photom

def test_dust_correct():

    correct = photom.extinction_correction('GMOS_S_r', 0.138)
    assert np.isclose(correct, 1.3818590723917497)

