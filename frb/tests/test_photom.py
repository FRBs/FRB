# Module to run tests on galaxy modules
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import importlib_resources

from astropy.table import Table
from astropy import units 
from astropy.io import fits 
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

from frb import frb
from frb.galaxies import photom
from frb.surveys.catalog_utils import convert_mags_to_flux, fill_masked


def test_dust_correct():

    correct = photom.extinction_correction('GMOS_S_r', 0.138)
    assert np.isclose(correct, 1.3936491887210187)

def test_flux_conversion():

    # Create a dummy table that should get converted in a known way
    tab = Table()
    tab['DES_r'] = [20.]
    tab['DES_r_err'] = 0.5

    fluxunits = 'mJy'

    fluxtab = convert_mags_to_flux(tab, fluxunits, exact_mag_err=False)

    # Check fluxes
    assert np.isclose(fluxtab['DES_r'], 0.036307805), "Check AB flux conversion."

    # Check errors
    assert np.isclose(fluxtab['DES_r_err'], 0.016720362110466937), "Check AB flux error."


def test_flux_conversion_flags():
    """Bad magnitudes and errors, and upper limits, are flagged consistently"""
    mags = [20., 24., 22., 36., -5., -99., -999.]
    errs = [0.5, 999., -99., 0.3, 0.1, -99., -999.]
    tab = Table()
    tab['DES_r'] = mags
    tab['DES_r_err'] = errs

    fluxtab = convert_mags_to_flux(tab, 'mJy')

    # A good measurement
    assert np.isclose(fluxtab['DES_r'][0], 0.036307805)
    # Upper limits (error = 999 in the host tables, -99 in the survey catalogs)
    #  keep their flux but not their error
    for ii in [1, 2]:
        assert fluxtab['DES_r'][ii] > 0
        assert fluxtab['DES_r_err'][ii] == -99.
    # Bad magnitudes are flagged, and so are their errors even if they look fine
    for ii in [3, 4, 5, 6]:
        assert fluxtab['DES_r'][ii] == -99.
        assert fluxtab['DES_r_err'][ii] == -99.


def test_flux_conversion_deep():
    """Faint HST/JWST photometry (mag > 30) is valid; placeholders are not"""
    tab = Table()
    tab['DES_r'] = [31., 33.5, 98., 99.99]
    tab['DES_r_err'] = [0.3, 0.2, 99., 99.99]

    fluxtab = convert_mags_to_flux(tab, 'mJy')

    # Deep photometry gets a flux and an error
    assert np.allclose(fluxtab['DES_r'][:2], 3630.7805e3 * 10**(-np.array([31., 33.5])/2.5))
    assert np.all(fluxtab['DES_r'][:2] > 0) and np.all(fluxtab['DES_r_err'][:2] > 0)
    # Placeholders are still flagged
    assert np.all(fluxtab['DES_r'][2:] == -99.) and np.all(fluxtab['DES_r_err'][2:] == -99.)


def test_flux_conversion_masked():
    """Masked entries (e.g. from merged survey catalogs) stay masked"""
    tab = Table(masked=True)
    tab['DES_r'] = np.ma.MaskedArray([20., 22.], mask=[False, True])
    tab['DES_r_err'] = np.ma.MaskedArray([0.5, 0.2], mask=[False, True])

    fluxtab = convert_mags_to_flux(tab, 'mJy')

    assert np.isclose(fluxtab['DES_r'][0], 0.036307805)
    assert fluxtab['DES_r'].mask[1] and fluxtab['DES_r_err'].mask[1]


def test_dust_correct_skips_missing():
    """-99. (surveys) and -999. (host tables) mean no measurement: no extinction correction"""
    tab = Table()
    tab['Name'] = ['FRB']
    tab['GMOS_S_r'] = [-99.]
    tab['GMOS_S_g'] = [-999.]
    tab['GMOS_S_i'] = [20.]
    code = photom.correct_photom_table(tab, 0.138, 'FRB')
    assert code == 0
    assert tab['GMOS_S_r'][0] == -99.
    assert tab['GMOS_S_g'][0] == -999.
    assert tab['GMOS_S_i'][0] < 20.  # brightened by the correction


def test_merge_photom_tables_narrow_strings():
    """Filling the merged photometry must not truncate the -999 of narrow string columns"""
    old = Table({'Name': ['A'], 'ra': [10.], 'dec': [20.],
                 'DES_r': [20.], 'DES_r_ref': ['x']})
    new = Table({'Name': ['B'], 'ra': [50.], 'dec': [-20.], 'WISE_W1': [18.]})

    merged = photom.merge_photom_tables(new, old)

    assert not merged.has_masked_values
    row_b = merged[merged['Name'] == 'B'][0]
    assert row_b['DES_r'] == -999.
    # DES_r_ref is a 1-character column holding a missing reference
    assert float(row_b['DES_r_ref']) == -999.
    assert float(merged[merged['Name'] == 'A'][0]['WISE_W1']) == -999.


def test_macquart_style_string_fill():
    """All-string tables (as in macquart_plot) fill to parsable text, whatever the width"""
    from astropy.table import vstack
    tab1 = Table({'FRB': ['FRB20121102A'], 'z': ['0.1'], 'DM': ['5']})
    tab2 = Table({'FRB': ['FRB20180924B'], 'z': ['0.32']})
    stacked = fill_masked(vstack([tab1, tab2]), -99.)
    assert not stacked.has_masked_values
    assert float(stacked['DM'][1]) == -99.
    assert np.all(stacked['DM'].astype(float) == [5., -99.])


def test_fractional_flux():
    isize = 5
    # FRB and HG
    frbname = 'FRB20180924B'
    frbdat = frb.FRB.by_name(frbname)
    # frbcoord = frbdat.coord
    hg = frbdat.grab_host()
    # Read cutout
    cutout_file = importlib_resources.files('frb.tests.files')/'FRB180924_cutout.fits'
    hdul = fits.open(cutout_file)

    hgcoord = hg.coord
    size = units.Quantity((isize, isize), units.arcsec)
    cutout = Cutout2D(hdul[0].data, hgcoord, size, wcs=WCS(hdul[0].header))

    # Run
    med_ff, sig_ff, f_weight = photom.fractional_flux(cutout, frbdat, hg)    

    assert np.isclose(sig_ff, 0.2906803236219953)


