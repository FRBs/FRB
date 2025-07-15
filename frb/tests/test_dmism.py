# Module to run tests on FRB calculations using DLAs

# TEST_UNICODE_LITERALS

import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord

from frb.dm import dm_ism_healpix_map
from frb import mw

def test_dm_ism_from_healpix():

    # Float
    b, l = 5., 50.

    dm_map = dm_ism_healpix_map.get_dm_map()
    dm_ism = dm_ism_healpix_map.dm_ism_from_healpix_map(l, b, dm_map)

    assert np.isclose(dm_ism, 321.538682)

    # NE2001
    #coord = SkyCoord(l=l * u.deg, b=b * u.deg, frame='galactic')
    #dm_ism_ne2001 = mw.ismDM(coord)


    # Array
    b = np.array([5., 10., 15.])
    l = np.array([50., 60., 70.])

    dm_ism = dm_ism_healpix_map.dm_ism_from_healpix_map(l, b, dm_map)
    
    assert np.isclose(dm_ism[0], 321.538682)