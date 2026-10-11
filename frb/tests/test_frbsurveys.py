# Module to run tests on surveys
#  Most of these are *not* done with Travis yet
# TEST_UNICODE_LITERALS

import numpy as np
import astropy
import pytest
import os, warnings, importlib

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.io.fits.hdu.image import PrimaryHDU

from frb.surveys import survey_utils
from frb.surveys.panstarrs import _ps1metadata
from frb.surveys import catalog_utils as cu

from PIL import Image

from numpy import setdiff1d

remote_data = pytest.mark.skipif(os.getenv('FRB_GDB') is None,
                                 reason='test requires dev suite')

nedlvs = pytest.mark.skipif('NEDLVS' not in os.environ, reason='Test reqires NEDLVS environment variable to be set.')


def _assert_masked_photometry(table):
    mag_cols, err_cols = cu._detect_mag_cols(table)

    for mag_col, err_col in zip(mag_cols, err_cols):
        masked = np.asarray(table[mag_col] == -99.0)
        if np.any(masked) and err_col in table.colnames:
            assert np.all(table[err_col][masked] == -99.0)

        good_photom = table[mag_col][~masked]
        assert np.all((good_photom > 0) & (good_photom < cu.MAX_VALID_MAG))

        if err_col in table.colnames:
            good_err = table[err_col][table[err_col] != -99.0]
            assert np.all((good_err > 0) & (good_err < 5))


def test_mask_bad_photometry_honors_mask():
    """Masked photometry is filled with the sentinel whatever data is under the mask"""
    for under in [np.nan, 0., 12., 99.99]:
        tab = Table(masked=True)
        tab['ra'] = [1., 2., 3.]
        tab['dec'] = [1., 2., 3.]
        tab['w4mag'] = np.ma.MaskedArray([12., under, 13.], mask=[False, True, False])
        tab['w4sigm'] = np.ma.MaskedArray([0.2, under, 9.0], mask=[False, True, False])
        pdict = {'WISE_W4': 'w4mag', 'WISE_W4_err': 'w4sigm', 'ra': 'ra', 'dec': 'dec'}
        out = cu.clean_cat(tab, pdict, mask_photometry=True)
        # Good measurement untouched
        assert out['WISE_W4'][0] == 12. and out['WISE_W4_err'][0] == 0.2
        # Masked entries are sentinels and not masked
        assert out['WISE_W4'][1] == -99. and out['WISE_W4_err'][1] == -99.
        assert not np.ma.is_masked(out['WISE_W4'][1])
        # Good magnitude with bad error: magnitude kept (upper limit), error flagged
        assert out['WISE_W4'][2] == 13. and out['WISE_W4_err'][2] == -99.


def _masked_survey_table():
    # An ID and a photometry column with masked entries, and metadata
    tab = Table(masked=True)
    tab['ra'] = [1., 2., 3.]
    tab['dec'] = [1., 2., 3.]
    tab['id'] = np.ma.MaskedArray([10, 11, 12], mask=[False, True, False])
    tab['mag'] = np.ma.MaskedArray([20., 21., 22.], mask=[False, True, False])
    tab['mag_err'] = np.ma.MaskedArray([0.1, 0.2, 0.3], mask=[False, True, False])
    tab.meta['survey'] = 'TEST'
    pdict = {'T_ID': 'id', 'T_g': 'mag', 'T_g_err': 'mag_err', 'ra': 'ra', 'dec': 'dec'}
    return tab, pdict


def test_clean_cat_fill_mask():
    """fill_mask and mask_photometry combine as documented in clean_cat"""
    # Rename only: masks are kept
    tab, pdict = _masked_survey_table()
    out = cu.clean_cat(tab, pdict)
    assert out['T_ID'].mask[1] and out['T_g'].mask[1]

    # fill_mask alone fills every masked entry (the main behavior)
    tab, pdict = _masked_survey_table()
    out = cu.clean_cat(tab, pdict, fill_mask=-99.)
    assert not out.has_masked_values
    assert out['T_ID'][1] == -99 and out['T_g'][1] == -99. and out['T_g_err'][1] == -99.
    assert out['T_g'][0] == 20.
    assert out.meta['survey'] == 'TEST'

    # mask_photometry alone: photometry gets -99, other masked entries stay masked
    tab, pdict = _masked_survey_table()
    out = cu.clean_cat(tab, pdict, mask_photometry=True)
    assert out['T_g'][1] == -99. and out['T_g_err'][1] == -99.
    assert not np.ma.is_masked(out['T_g'][1])
    assert out['T_ID'].mask[1]

    # Both, with a non-default value: used for the photometry and everything else
    tab, pdict = _masked_survey_table()
    out = cu.clean_cat(tab, pdict, fill_mask=-999., mask_photometry=True)
    assert not out.has_masked_values
    assert out['T_g'][1] == -999. and out['T_g_err'][1] == -999. and out['T_ID'][1] == -999

    # An unmasked table is untouched
    tab = Table({'ra': [1., 2.], 'dec': [1., 2.], 'mag': [20., 21.]})
    out = cu.clean_cat(tab, {'T_g': 'mag'}, fill_mask=-99.)
    assert np.all(out['T_g'] == [20., 21.])


def _masked_col_table(values, mask, **kwargs):
    tab = Table(masked=True)
    tab['c'] = np.ma.MaskedArray(values, mask=mask)
    return tab


def test_fill_masked_strings():
    """String columns are widened, not truncated, when filled with a number"""
    # Widths narrower than, equal to and wider than the text '-99.0'
    for values, expected in [(['t', 'f'], '-99.0'), (['ab', 'cd'], '-99.0'),
                             (['STAR', 'QSO'], '-99.0'), (['GALAXY', 'x'], '-99.0'),
                             (['DESI J123456.7+012345', 'x'], '-99.0')]:
        tab = _masked_col_table(values, [False, True])
        out = cu.fill_masked(tab, -99.)
        assert out['c'][0] == values[0]
        assert out['c'][1] == expected
        assert not out.has_masked_values
        assert tab.has_masked_values  # Input is untouched
    # The sentinel survives a round trip through float() (what chk_fill relies on)
    out = cu.fill_masked(_masked_col_table(['abc', 'def'], [False, True]), -999.)
    assert float(out['c'][1]) == -999.
    # Bytes and object columns
    out = cu.fill_masked(_masked_col_table(np.array([b'STAR', b'QSO']), [False, True]), -99.)
    assert out['c'].dtype.kind == 'S'  # Still bytes (astropy shows the items as str)
    assert out['c'][1] == '-99.0' and out['c'][0] == 'STAR'
    out = cu.fill_masked(_masked_col_table(np.array(['a', 'bcd'], dtype=object), [False, True]), -99.)
    assert out['c'][1] == '-99.0' and isinstance(out['c'][1], str)
    # Override the text
    out = cu.fill_masked(_masked_col_table(['ab', 'cd'], [False, True]), -99., str_fill='')
    assert out['c'][1] == ''


def test_fill_masked_other_dtypes():
    # Numeric columns: filled, with units and metadata kept
    tab = Table(masked=True)
    tab['f'] = np.ma.MaskedArray([1., 2.], mask=[False, True])
    tab['i'] = np.ma.MaskedArray([1, 2], mask=[False, True])
    tab['f'].unit = units.mag
    tab.meta['survey'] = 'TEST'
    out = cu.fill_masked(tab, -99.)
    assert out['f'][1] == -99. and out['i'][1] == -99
    assert out['f'].unit == units.mag and out.meta['survey'] == 'TEST'
    assert not out.has_masked_values
    # Columns that cannot hold the sentinel stay masked, without errors or fake values
    tab = Table(masked=True)
    tab['b'] = np.ma.MaskedArray([True, False], mask=[False, True])
    tab['u'] = np.ma.MaskedArray(np.array([1, 2], dtype=np.uint8), mask=[False, True])
    tab['f'] = np.ma.MaskedArray([1., 2.], mask=[False, True])
    out = cu.fill_masked(tab, -99.)
    assert out['b'].mask[1] and out['u'].mask[1]
    assert out['f'][1] == -99. and not np.ma.is_masked(out['f'][1])
    # An unmasked table is returned as is
    tab = Table({'a': [1., 2.], 's': ['x', 'y']})
    out = cu.fill_masked(tab, -99.)
    assert np.all(out['a'] == [1., 2.]) and list(out['s']) == ['x', 'y']


def test_clean_cat_fill_mask_strings():
    tab = Table(masked=True)
    tab['ra'] = [1., 2.]
    tab['dec'] = [1., 2.]
    tab['type'] = np.ma.MaskedArray(['PSF', 'DEV'], mask=[False, True])
    out = cu.clean_cat(tab, {'T_type': 'type', 'ra': 'ra', 'dec': 'dec'}, fill_mask=-99.)
    assert out['T_type'][0] == 'PSF' and out['T_type'][1] == '-99.0'


def test_xmatch_and_merge_fill_value_strings():
    # A string column that is missing for one of the sources (and is narrow)
    des = Table({'ra': [10.0000, 10.0100], 'dec': [20.0000, 20.0100],
                 'DES_r': [21., 22.], 'DES_type': ['PSF', 'DEV']})
    wise = Table({'ra': [10.00001, 10.0200], 'dec': [20.00001, 20.0200],
                  'WISE_W1': [18., 19.]})
    merged = cu.xmatch_and_merge_cats(des, wise, tol=1*units.arcsec, fill_value=-999.)
    assert not merged.has_masked_values
    assert sorted(merged['DES_type']) == ['-999.0', 'DEV', 'PSF']
    assert np.sum(merged['DES_r'] == -999.) == 1


def test_ensure_empty_schema_dtypes():
    columns = ['ra', 'dec', 'ID', 'type', 'mag']
    out = cu.ensure_empty_schema(Table(), columns, dtypes={'ID': int, 'type': str})
    assert len(out) == 0 and out.colnames == columns
    assert out['ra'].dtype.kind == 'f' and out['mag'].dtype.kind == 'f'
    assert out['ID'].dtype.kind == 'i' and out['type'].dtype.kind == 'U'
    # Without dtypes everything is a float, as before
    out = cu.ensure_empty_schema(Table(), columns)
    assert all(out[col].dtype.kind == 'f' for col in columns)
    # Existing columns of the wrong kind are converted, others are kept, extras are dropped
    tab = Table({'ra': np.array([], dtype=float), 'ID': np.array([], dtype=float),
                 'type': np.array([], dtype=object), 'extra': np.array([], dtype=float)})
    tab['ra'].unit = units.deg
    out = cu.ensure_empty_schema(tab, columns, dtypes={'ID': 'uint64', 'type': str})
    assert out['ID'].dtype == np.dtype('uint64') and out['type'].dtype.kind == 'U'
    assert out['ra'].unit == units.deg and 'extra' not in out.colnames
    # A catalog with rows is left alone
    tab = Table({'ra': [1.], 'ID': [1.5]})
    out = cu.ensure_empty_schema(tab, ['ra', 'ID'], dtypes={'ID': int})
    assert out['ID'].dtype.kind == 'f' and len(out) == 1


@pytest.mark.parametrize('survey,modname', [
    ('DECaL', 'decals'), ('DES', 'des'), ('DELVE', 'delve'), ('NSC', 'nsc'),
    ('VISTA', 'vista'), ('SDSS', 'sdss'), ('GALEX', 'galex'),
    ('Pan-STARRS', 'panstarrs'), ('Euclid', 'euclid'), ('DESI', 'desi')])
def test_survey_schema_dtypes(survey, modname):
    """The declared dtypes of every survey refer to columns of its schema and are honored"""
    module = importlib.import_module(f'frb.surveys.{modname}')
    schema = getattr(module, 'photom', {}).get(survey) or module.spectrom[survey]
    dtypes = module.schema_dtypes[survey]
    assert len(dtypes) > 0 and set(dtypes) <= set(schema)
    empty = cu.ensure_empty_schema(Table(), list(schema), dtypes=dtypes)
    for col in schema:
        kind = np.dtype(dtypes[col]).kind if col in dtypes else 'f'
        assert empty[col].dtype.kind == kind, col


def test_mask_bad_photometry_deep():
    """Deep (HST/JWST) magnitudes survive the masking, placeholders do not"""
    tab = Table()
    tab['ra'] = [1., 2., 3., 4.]
    tab['dec'] = [1., 2., 3., 4.]
    tab['mag'] = [31., 34.9, 35.1, 99.99]
    tab['mag_err'] = [0.3, 0.2, 0.2, 99.99]
    pdict = {'F_x': 'mag', 'F_x_err': 'mag_err', 'ra': 'ra', 'dec': 'dec'}
    out = cu.clean_cat(tab, pdict, mask_photometry=True)
    assert np.allclose(out['F_x'], [31., 34.9, -99., -99.])
    assert np.allclose(out['F_x_err'], [0.3, 0.2, -99., -99.])


def _merge_tables():
    # One source in both tables, one in each table only
    des = Table({'ra': [10.0000, 10.0100], 'dec': [20.0000, 20.0100],
                 'DES_r': [21., 22.], 'DES_r_err': [0.1, 0.2]})
    wise = Table({'ra': [10.00001, 10.0200], 'dec': [20.00001, 20.0200],
                  'WISE_W1': [18., 19.], 'WISE_W1_err': [0.05, 0.06]})
    return des, wise


def test_xmatch_and_merge_masks_by_default():
    des, wise = _merge_tables()
    merged = cu.xmatch_and_merge_cats(des, wise, tol=1*units.arcsec)
    assert len(merged) == 3
    # Sources missing from one of the tables are masked, not given fake values
    assert np.sum(merged['DES_r'].mask) == 1
    assert np.sum(merged['WISE_W1'].mask) == 1
    # and they stay masked through the flux conversion
    flux = cu.convert_mags_to_flux(merged)
    assert np.sum(flux['DES_r'].mask) == 1 and np.sum(flux['WISE_W1_err'].mask) == 1


def test_xmatch_and_merge_fill_value():
    des, wise = _merge_tables()
    merged = cu.xmatch_and_merge_cats(des, wise, tol=1*units.arcsec, fill_value=-999.)
    assert not np.ma.is_masked(merged['DES_r'])
    assert np.sum(merged['DES_r'] == -999.) == 1
    assert np.sum(merged['WISE_W1_err'] == -999.) == 1
    # -999. is "no measurement": both flux and error are flagged
    flux = cu.convert_mags_to_flux(merged)
    missing = merged['DES_r'] == -999.
    assert np.all(flux['DES_r'][missing] == -99.) and np.all(flux['DES_r_err'][missing] == -99.)


def _assert_empty_catalog(survey_name, catalog_kwargs=None):
    catalog_kwargs = {} if catalog_kwargs is None else catalog_kwargs

    candidates = [
        (SkyCoord(l=0., b=0., unit='deg', frame='galactic').transform_to('icrs'), 1 * units.arcsec),
        (SkyCoord(0., 90., unit='deg', frame='icrs'), 1 * units.arcsec),
        (SkyCoord(0., -90., unit='deg', frame='icrs'), 1 * units.arcsec),
        (SkyCoord(l=0., b=0., unit='deg', frame='galactic').transform_to('icrs'), 0.5 * units.arcsec),
    ]

    rng = np.random.default_rng(seed=0)
    for _ in range(10):
        coord = SkyCoord(rng.uniform(0., 360.), rng.uniform(-90., 90.), unit='deg', frame='icrs')
        candidates.append((coord, 1 * units.arcsec))

    for coord, radius in candidates:
        try:
            survey = survey_utils.load_survey_by_name(survey_name, coord, radius)
            empty_tbl = survey.get_catalog(**catalog_kwargs)
        except Exception:
            continue

        if isinstance(empty_tbl, Table) and len(empty_tbl) == 0:
            assert isinstance(empty_tbl, Table)
            assert len(empty_tbl) == 0
            assert len(empty_tbl.colnames) > 2
            # Non-float columns (IDs, strings) have their proper dtype even if empty
            module = importlib.import_module(type(survey).__module__)
            for col, dtype in getattr(module, 'schema_dtypes', {}).get(survey.survey, {}).items():
                assert empty_tbl[col].dtype.kind == np.dtype(dtype).kind, col
            return empty_tbl

    pytest.fail(f'Could not find an empty catalog for {survey_name}')

@remote_data
def test_sdss():
    #coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    coord = SkyCoord(0, 0, unit=units.deg)
    search_r = 1 * units.arcmin
    img_r = 30 * units.arcsec
    # Instantiate
    sdss_srvy = survey_utils.load_survey_by_name('SDSS', coord, search_r)
    sdss_tbl = sdss_srvy.get_catalog()
    # Reasonable behavior
    assert isinstance(sdss_tbl, Table)
    assert len(sdss_tbl) == 73

    # Test masking
    _assert_masked_photometry(sdss_tbl)

    # Test empty table handling
    _assert_empty_catalog('SDSS')

    imghdu = sdss_srvy.get_image(imsize=img_r, band='r')
    assert isinstance(imghdu, PrimaryHDU)
    assert imghdu.data is not None
    assert imghdu.data.ndim == 2
    assert imghdu.data.shape == (76, 76)


@remote_data
def test_wise():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 1 * units.arcmin
    img_r = 10 * units.arcsec

    wise_srvy = survey_utils.load_survey_by_name('WISE', coord, search_r)
    wise_tbl = wise_srvy.get_catalog()
    #
    assert isinstance(wise_tbl, Table)
    assert len(wise_tbl) == 15
    _assert_masked_photometry(wise_tbl)

    _assert_empty_catalog('WISE')


    # Test canonical FITS path
    imghdu = wise_srvy.get_image(imsize=img_r, band="W1")
    assert isinstance(imghdu, PrimaryHDU)
    assert imghdu.data.shape == (5,5)

    # Deprecated compatibility alias
    with pytest.warns(DeprecationWarning):
        alias_hdu = wise_srvy.get_cutout(imsize=img_r, band="W1")
    assert isinstance(alias_hdu, PrimaryHDU)

# THIS TEST IS NOW BROKEN
'''
@remote_data
def test_psrcat():
    # Catalog
    coord = SkyCoord('J000604.8+183459', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    psrcat_srvy = survey_utils.load_survey_by_name('PSRCAT', coord, search_r)
    pulsars = psrcat_srvy.get_catalog()
    #
    assert isinstance(pulsars, Table)
    assert len(pulsars) == 1
'''


@remote_data
def test_des():
    # Catalog
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    search_r = 30 * units.arcsec
    img_r = 10 * units.arcsec

    des_srvy = survey_utils.load_survey_by_name('DES', coord, search_r)
    des_tbl = des_srvy.get_catalog(print_query=True)
    assert isinstance(des_tbl, Table)
    assert len(des_tbl) == 26
    _assert_masked_photometry(des_tbl)

    _assert_empty_catalog('DES')

    
    # Canonical FITS path
    imghdu = des_srvy.get_image(imsize=img_r, band="g")
    assert isinstance(imghdu, PrimaryHDU)
    assert imghdu.data.shape == (39,39)

    # Deprecated compatibility alias
    with pytest.warns(DeprecationWarning):
        data, hdr = des_srvy.get_cutout(imsize=img_r, band="g")
    assert data.shape == (39,39)

@remote_data
def test_desi():
    # Catalog
    coord = SkyCoord(0, 0, unit="deg")
    search_r = 0.3 * units.arcmin # Can't go below this with Noirlab for some reason. No error.
                                  #Just a constant number of entries returned. 0 arcmin does give 0 entries though.
    desi_srvy = survey_utils.load_survey_by_name('DESI', coord, search_r)
    desi_tbl = desi_srvy.get_catalog(print_query=True, exclude_stars=True, zcat_primary_only=True)
    assert isinstance(desi_tbl, Table)
    assert len(desi_tbl) == 3230

    _assert_empty_catalog('DESI', catalog_kwargs={'exclude_stars': True, 'zcat_primary_only': True})


def test_euclid():
    from astropy.io import fits
    coord = SkyCoord("17h51m07.4s +65d31m50.8s", frame='icrs')
    search_r = 60 * units.arcsec

    euclid_srvy = survey_utils.load_survey_by_name('Euclid', coord, search_r)
    euclid_tbl = euclid_srvy.get_catalog(check_spectra=True, timeout=30)

    assert isinstance(euclid_tbl, Table)
    assert len(euclid_tbl) == 155
    _assert_masked_photometry(euclid_tbl)
    assert euclid_tbl.meta['survey'] == 'Euclid'
    assert 'ra' in euclid_tbl.colnames
    assert 'dec' in euclid_tbl.colnames
    assert 'Euclid_has_spectrum' in euclid_tbl.colnames

    _assert_empty_catalog('Euclid', catalog_kwargs={'check_spectra': True, 'timeout': 30})


    # Canonical FITS path
    image, image_hdr = euclid_srvy.get_image(imsize=2*units.arcmin, timeout=30)
    assert isinstance(image, np.ndarray)
    assert isinstance(image_hdr, fits.Header)
    assert image.shape == (1207, 1207)

    # Deprecated compatibility alias
    with pytest.warns(DeprecationWarning):
        cutout, cutout_hdr = euclid_srvy.get_cutout(imsize=2*units.arcmin, timeout=30)
    assert isinstance(cutout, np.ndarray)
    assert isinstance(cutout_hdr, fits.Header)
    assert cutout.shape == (1207, 1207)


@remote_data
def test_nsc():
    # Catalog
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    search_r = 60 * units.arcsec
    img_r = 10 * units.arcsec

    nsc_srvy = survey_utils.load_survey_by_name('NSC', coord, search_r)
    nsc_tbl = nsc_srvy.get_catalog(print_query=True)
    #
    assert isinstance(nsc_tbl, Table)
    assert len(nsc_tbl) == 43
    _assert_masked_photometry(nsc_tbl)

    _assert_empty_catalog('NSC')


    # Canonical FITS path
    imghdu = nsc_srvy.get_image(imsize=img_r, band="g")
    assert isinstance(imghdu, PrimaryHDU)
    assert imghdu.data.shape == (38,38)

    # Deprecated compatibility alias
    with pytest.warns(DeprecationWarning):
        data, hdr = nsc_srvy.get_cutout(imsize=img_r, band="g")
    assert data.shape == (38,38)

@remote_data
def test_hsc():
    # Catalog
    coord = SkyCoord(0,0, unit="deg")
    search_r = 30 * units.arcsec

    hsc_srvy = survey_utils.load_survey_by_name('HSC', coord, search_r)
    hsc_tbl = hsc_srvy.get_catalog(print_query=True)
    #
    assert isinstance(hsc_tbl, Table)
    assert len(hsc_tbl) == 64

    _assert_masked_photometry(hsc_tbl)

    _assert_empty_catalog('HSC')


@remote_data
def test_delve():
    # Catalog
    coord = SkyCoord("J102922+012133", unit=(units.hourangle, units.deg))
    search_r = 120 * units.arcsec

    delve_srvy = survey_utils.load_survey_by_name('DELVE', coord, search_r)
    delve_tbl = delve_srvy.get_catalog(print_query=True)
    #
    assert isinstance(delve_tbl, Table)
    assert len(delve_tbl) == 178
    _assert_masked_photometry(delve_tbl)

    _assert_empty_catalog('DELVE')


    # No image service available for DELVE

@remote_data
def test_vista():
    # Catalog
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    #coord = SkyCoord('J210000-400000', unit=(units.hourangle, units.deg))
    search_r = 120 * units.arcsec
    img_r = 120 * units.arcsec

    vista_srvy = survey_utils.load_survey_by_name('VISTA', coord, search_r)
    vista_tbl = vista_srvy.get_catalog(print_query=True)
    #
    assert isinstance(vista_tbl, Table)
    assert len(vista_tbl) == 152
    _assert_masked_photometry(vista_tbl)

    _assert_empty_catalog('VISTA')

    # VSA image retrieval is service-dependent; if unavailable the method should fail gracefully.
    imghdu = vista_srvy.get_image(imsize=img_r, band='J', timeout=120)
    assert isinstance(imghdu, PrimaryHDU)
    assert imghdu.header.get('ESO INS FILT1 NAME') == 'J'
    assert imghdu.data.ndim == 2
    assert imghdu.data.shape == (354, 354)


@remote_data
def test_decals():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 60 * units.arcsec
    img_r = 10 * units.arcsec

    decal_srvy = survey_utils.load_survey_by_name('DECaL', coord, search_r)
    decal_tbl = decal_srvy.get_catalog(print_query=True)
    #
    assert isinstance(decal_tbl, Table)
    assert len(decal_tbl) == 69
    _assert_masked_photometry(decal_tbl)

    _assert_empty_catalog('DECaL')


    # Canonical FITS path
    imghdu = decal_srvy.get_image(imsize=img_r, band="g")
    assert isinstance(imghdu, PrimaryHDU)
    assert imghdu.data.shape == (39, 39)

    # Deprecated compatibility alias
    with pytest.warns(DeprecationWarning):
        data, hdr = decal_srvy.get_cutout(imsize=img_r, band="g")
    assert data.shape == (39, 39)



@remote_data
def test_first():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec
    #
    first_srvy = survey_utils.load_survey_by_name('FIRST', coord, search_r)
    first_tbl = first_srvy.get_catalog()
    #
    assert isinstance(first_tbl, Table)
    assert len(first_tbl) == 1

    _assert_empty_catalog('FIRST')

    # Imaging from SkyView
    imghdu = first_srvy.get_image(imsize=10*units.arcsec)
    assert isinstance(imghdu, PrimaryHDU)
    assert imghdu.data is not None
    assert imghdu.data.ndim == 2
    assert imghdu.data.shape == (6, 6)


@remote_data
def test_panstarrs():
    #Test get_catalog
    coord = SkyCoord(0., 0.,unit="deg")
    search_r = 120*units.arcsec
    ps_survey = survey_utils.load_survey_by_name('Pan-STARRS',coord,search_r)
    ps_table = ps_survey.get_catalog(photoz=True)

    assert isinstance(ps_table, Table)
    assert len(ps_table) == 161
    _assert_masked_photometry(ps_table)

    assert 'z_phot' in ps_table.colnames
    assert 'z_photErr' in ps_table.colnames

    _assert_empty_catalog('Pan-STARRS', catalog_kwargs={'photoz': True})

    #Test get_cutout
    # Default imsize for both methods below: 30 arcsec
    cutout, = ps_survey.get_cutout()
    assert isinstance(cutout,Image.Image)
    assert cutout.size == (120,120)

    #Test get_image
    imghdu = ps_survey.get_image()
    assert isinstance(imghdu,PrimaryHDU)
    assert imghdu.data.shape == (120,120)

    # Deprecated Pan-STARRS arg alias
    with pytest.warns(DeprecationWarning):
        imghdu_depr = ps_survey.get_image(filt='i')
    assert isinstance(imghdu_depr, PrimaryHDU)

    # Test getting metadata repeatedly to check caching
    for index in range(10):
        metadata = _ps1metadata()
        assert isinstance(metadata,Table)
        assert len(metadata) > 0
        assert np.all(np.isin(metadata.colnames, ['name', 'datatype', 'description']))

@nedlvs
def test_nedlvs():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcmin

    # Test get_catalog
    nedlvs_srvy = survey_utils.load_survey_by_name('NEDLVS', coord, search_r)
    nedlvs_tbl = nedlvs_srvy.get_catalog()
    assert isinstance(nedlvs_tbl, Table)
    # Remote NEDLVS content can grow over time; require at least the historical matches.
    assert len(nedlvs_tbl) == 3
    _assert_empty_catalog('NEDLVS')

@remote_data
def test_tully():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 90 * units.deg

    # Test get_catalog
    tully_srvy = survey_utils.load_survey_by_name('TullyGroupCat', coord, search_r)
    tully_tbl = tully_srvy.get_catalog(transverse_distance_cut=5*units.Mpc)
    assert isinstance(tully_tbl, Table)
    assert len(tully_tbl) == 6

@remote_data
def test_galex():
    coord = SkyCoord('J142532.38+120121.17', unit=(units.hourangle, units.deg))
    search_r = 240 * units.arcsec
    img_r = 60 * units.arcsec

    # Test get_catalog
    galex_srvy = survey_utils.load_survey_by_name('GALEX', coord, search_r)
    galex_tbl = galex_srvy.get_catalog()
    assert isinstance(galex_tbl, Table)
    assert len(galex_tbl) == 194
    _assert_masked_photometry(galex_tbl)

    _assert_empty_catalog('GALEX')

    imghdu = galex_srvy.get_image(imsize=img_r, band='NUV')
    assert isinstance(imghdu, PrimaryHDU)
    assert imghdu.data is not None
    assert imghdu.data.ndim == 2
    assert imghdu.data.shape == (40, 40)


@remote_data
def test_2mass():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 240 * units.arcsec
    img_r = 60 * units.arcsec

    # Test get_catalog
    mass_srvy = survey_utils.load_survey_by_name('2MASS', coord, search_r)
    mass_tbl = mass_srvy.get_catalog()
    assert isinstance(mass_tbl, Table)
    assert len(mass_tbl) == 41
    _assert_masked_photometry(mass_tbl)

    _assert_empty_catalog('2MASS')

    imghdu = mass_srvy.get_image(imsize=img_r, band='J')
    assert isinstance(imghdu, PrimaryHDU)
    assert imghdu.data is not None
    assert imghdu.data.ndim == 2
    assert imghdu.data.shape == (60, 60)


@remote_data
def test_in_which_survey():
    """
    To test if `survey_utils.in_which_survey` works.
    """
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    
    with warnings.catch_warnings(record=True) as allwarns:
        inside = survey_utils.in_which_survey(coord, optical_only=False)
    expected_dict = {'Pan-STARRS': True,
                     'WISE': True,
                     'SDSS': True,
                     'DES': False,
                     'DESI': False,
                     'DELVE': True,
                     'DECaL': True,
                     'Euclid': False,
                     'VISTA': False,
                     'NSC': True,
                     'HSC': False,
                     'NEDLVS': True,
                     '2MASS': True,
                     'GALEX': True,
                     'NVSS': False,
                     'FIRST': False,
                     'WENSS': False}

    for key in inside.keys():
        assert expected_dict[key] == inside[key], "{} did not match expectations.".format(key)
    
    # Test if warnings were produced the correct number of times.
    # Only for stable versions. For some reason, the 4.3dev version
    # returns empty table for the Heasarc surveys but 4.2 returns 1 or 2 objects.
    # Strange.
    if 'dev' not in astropy.__version__:
        warncount = 0
        for w in allwarns:
            if "Check location manually" in w.message.args[0]:
                warncount += 1
        assert warncount == 2

@remote_data
def test_search_all():
    """
    Test if survey_utils.search_all_surveys() works
    """
    # Small radius as it might fail when
    # merging 0 length catalogs
    radius = 5*units.arcsec
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    combined_cat = survey_utils.search_all_surveys(coord, radius=radius)
    assert len(combined_cat)==2


    # Nothing from NEDLVS and so not in the combined catalog
    colnames = ['ra', 'dec', 'separation',
                '2MASS_ID', '2MASS_h', '2MASS_h_err', '2MASS_j', '2MASS_j_err', '2MASS_k', '2MASS_k_err',
                'DECaL_ID', 'DECaL_brick', 'DECaL_g', 'DECaL_g_err', 'DECaL_r', 'DECaL_r_err', 'DECaL_type', 'DECaL_z', 'DECaL_z_err',
                'DELVE_ID', 'DELVE_g', 'DELVE_g_err', 'DELVE_i', 'DELVE_i_err', 'DELVE_r', 'DELVE_r_err', 'DELVE_z', 'DELVE_z_err',
                'DESI_ID', 'DESI_name', 'DESI_specsubtype', 'DESI_spectype', 'DESI_survey', 'DESI_z', 'DESI_z_err', 'DESI_z_warn', 'DESI_zcat_nspec', 'DESI_zcat_primary',
                'NSC_ID', 'NSC_VR', 'NSC_VR_err', 'NSC_Y', 'NSC_Y_err', 'NSC_g', 'NSC_g_err', 'NSC_i', 'NSC_i_err', 'NSC_r', 'NSC_r_err', 'NSC_u', 'NSC_u_err', 'NSC_z', 'NSC_z_err',
                'Pan-STARRS_ID', 'Pan-STARRS_g', 'Pan-STARRS_g_err', 'Pan-STARRS_i', 'Pan-STARRS_i_err', 'Pan-STARRS_r', 'Pan-STARRS_r_err', 'Pan-STARRS_y', 'Pan-STARRS_y_err', 'Pan-STARRS_z', 'Pan-STARRS_z_err',
                'SDSS_ID', 'SDSS_field', 'SDSS_g', 'SDSS_g_err', 'SDSS_i', 'SDSS_i_err', 'SDSS_r', 'SDSS_r_err', 'SDSS_u', 'SDSS_u_err', 'SDSS_z', 'SDSS_z_err',
                'WISE_W1', 'WISE_W1_err', 'WISE_W2', 'WISE_W2_err', 'WISE_W3', 'WISE_W3_err', 'WISE_W4', 'WISE_W4_err',
                'camcol', 'class', 'class_star', 'class_star_g', 'class_star_i', 'class_star_r', 'class_star_z', 'ebv', 'extinction_g', 'extinction_i', 'extinction_r', 'extinction_u', 'extinction_z', 'gPSFmag', 'gPSFmagErr', 'iPSFmag', 'iPSFmagErr', 'objInfoFlag', 'photo_z', 'photo_zerr', 'qualityFlag', 'rKronRad', 'rPSFmag', 'rPSFmagErr', 'rerun', 'run', 'source_id', 'survey', 'tmass_key', 'type', 'yPSFmag', 'yPSFmagErr', 'zPSFmag', 'zPSFmagErr', 'z_phot', 'z_photErr', 'z_phot_l68', 'z_phot_l95', 'z_phot_median', 'z_phot_u68', 'z_phot_u95', 'z_spec', 'z_spec_DECaL']
    assert len(setdiff1d(combined_cat.colnames, colnames))==0
    assert combined_cat['Pan-STARRS_ID'].mask[1] # This is a DECaLS source without a Pan-STARRS match.