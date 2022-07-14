# Module to run tests on galaxy modules
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np
from pkg_resources import resource_filename

from astropy.table import Table
from astropy import units 
from astropy.io import fits 
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

from frb import frb
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
    tab['VISTA_Y'] = 20.
    tab['VISTA_Y_err'] = 0.5

    fluxunits = 'mJy'

    fluxtab = convert_mags_to_flux(tab, fluxunits)

    # Check fluxes
    assert np.isclose(fluxtab['DES_r'], 0.036307805), "Check AB flux conversion."
    # WISE conversion is now done at the survey level
    #assert np.isclose(fluxtab['WISE_W1'], 0.0030954), "Check WISE flux conversion."
    assert np.isclose(fluxtab['VISTA_Y'], 0.0208732), "Check VISTA flux conversion."

    # Check errors
    assert np.isclose(fluxtab['DES_r_err'], 0.02123618797770558), "Check AB flux error."
    #assert np.isclose(fluxtab['WISE_W1_err'], 0.0018104783879441312), "Check WISE flux error."
    assert np.isclose(fluxtab['VISTA_Y_err'], 0.012208592584879318), "Check VISTA flux error."


def test_fractional_flux():
    isize = 5
    # FRB and HG
    frbname = 'FRB20180924B'
    frbdat = frb.FRB.by_name(frbname)
    # frbcoord = frbdat.coord
    hg = frbdat.grab_host()
    # Read cutout
    cutout_file = os.path.join(resource_filename('frb','tests'), 'files',
                               'FRB180924_cutout.fits')
    hdul = fits.open(cutout_file)

    hgcoord = hg.coord
    size = units.Quantity((isize, isize), units.arcsec)
    cutout = Cutout2D(hdul[0].data, hgcoord, size, wcs=WCS(hdul[0].header))

    # Run
    med_ff, sig_ff, f_weight = photom.fractional_flux(cutout, frbdat, hg)    

    assert np.isclose(sig_ff, 0.2906803236219953)


