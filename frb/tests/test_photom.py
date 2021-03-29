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
from frb.surveys.catalog_utils import convert_mags_to_flux

def test_dust_correct():

    correct = photom.extinction_correction('GMOS_S_r', 0.138)
    assert np.isclose(correct, 1.3818590723917497)

def test_flux_conversion():

    # Create a dummy table that should get converted in a known way
    tab = Table()
    tab['DES_r'] = [20.]
    tab['DES_r_err'] = 0.5
    tab['WISE_W1'] = 20.
    tab['WISE_W1_err'] = 0.5
    tab['VISTA_Y'] = 20.
    tab['VISTA_Y_err'] = 0.5

    fluxunits = 'mJy'

    fluxtab = convert_mags_to_flux(tab, fluxunits)

    # Check fluxes
    assert np.isclose(fluxtab['DES_r'], 0.036307805), "Check AB flux conversion."
    assert np.isclose(fluxtab['WISE_W1'], 0.0030954), "Check WISE flux conversion."
    assert np.isclose(fluxtab['VISTA_Y'], 0.0208732), "Check VISTA flux conversion."

    # Check errors
    assert np.isclose(fluxtab['DES_r_err'], 0.02123618797770558), "Check AB flux error."
    assert np.isclose(fluxtab['WISE_W1_err'], 0.0018104783879441312), "Check WISE flux error."
    assert np.isclose(fluxtab['VISTA_Y_err'], 0.012208592584879318), "Check VISTA flux error."