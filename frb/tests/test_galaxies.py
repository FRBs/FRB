# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.table import Table
from astropy import units
from astropy.coordinates import SkyCoord

from frb.galaxies import nebular


def test_ebv():

    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    # Set width
    ebv = nebular.get_ebv(coord)
    # Test
    assert isinstance(ebv, dict)
    for key in ['meanValue', 'std', 'minValue']:
        assert key in ebv.keys()
    assert np.isclose(float(ebv['meanValue']), 0.0172)
