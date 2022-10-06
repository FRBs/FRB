# Module to run tests on surveys
#  Most of these are *not* done with Travis yet
# TEST_UNICODE_LITERALS

import astropy
import pytest
import os, warnings

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.io.fits.hdu.image import PrimaryHDU

from frb.surveys import survey_utils
from PIL import Image

remote_data = pytest.mark.skipif(os.getenv('FRB_GDB') is None,
                                 reason='test requires dev suite')

@remote_data
def test_sdss():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec
    #
    sdss_srvy = survey_utils.load_survey_by_name('SDSS', coord, search_r)
    sdss_tbl = sdss_srvy.get_catalog()
    #
    assert isinstance(sdss_tbl, Table)
    assert len(sdss_tbl) == 2

@remote_data
def test_wise():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    wise_srvy = survey_utils.load_survey_by_name('WISE', coord, search_r)
    wise_tbl = wise_srvy.get_catalog()
    #
    assert isinstance(wise_tbl, Table)
    assert len(wise_tbl) == 1

    #Test get_image
    imghdu = wise_srvy.get_cutout(imsize=search_r, filter="W1")
    assert isinstance(imghdu,PrimaryHDU)
    assert imghdu.data.shape == (5,5)

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


@remote_data
def test_des():
    # Catalog
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    des_srvy = survey_utils.load_survey_by_name('DES', coord, search_r)
    des_tbl = des_srvy.get_catalog(print_query=True)
    #
    assert isinstance(des_tbl, Table)
    assert len(des_tbl) == 2
@remote_data
def test_nsc():
    # Catalog
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    nsc_srvy = survey_utils.load_survey_by_name('NSC', coord, search_r)
    nsc_tbl = nsc_srvy.get_catalog(print_query=True)
    #
    assert isinstance(nsc_tbl, Table)
    assert len(nsc_tbl) == 1

@remote_data
def test_vista():
    # Catalog
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    vista_srvy = survey_utils.load_survey_by_name('VISTA', coord, search_r)
    vista_tbl = vista_srvy.get_catalog(print_query=True)
    #
    assert isinstance(vista_tbl, Table)
    assert len(vista_tbl) == 1

@remote_data
def test_vista():
    # Catalog
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    vista_srvy = survey_utils.load_survey_by_name('DES', coord, search_r)
    vista_tbl = vista_srvy.get_catalog(print_query=True)
    #
    assert isinstance(vista_tbl, Table)
    assert len(vista_tbl) == 1



@remote_data
def test_decals():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    decal_srvy = survey_utils.load_survey_by_name('DECaL', coord, search_r)
    decal_tbl = decal_srvy.get_catalog(print_query=True)
    #
    assert isinstance(decal_tbl, Table)
    assert len(decal_tbl) == 3



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


@remote_data
def test_panstarrs():
    #Test get_catalog
    coord = SkyCoord(0, 0,unit="deg")
    search_r = 30*units.arcsec
    ps_survey = survey_utils.load_survey_by_name('Pan-STARRS',coord,search_r)
    ps_table = ps_survey.get_catalog()

    assert isinstance(ps_table, Table)
    assert len(ps_table) == 7

    #Test get_cutout
    cutout, = ps_survey.get_cutout()
    assert isinstance(cutout,Image.Image)
    assert cutout.size == (120,120)

    #Test get_image
    imghdu = ps_survey.get_image()
    assert isinstance(imghdu,PrimaryHDU)
    assert imghdu.data.shape == (120,120)

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
                    'DECaL': True,
                    'VISTA': False,
                    'NSC': True,
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
    colnames = ['Pan-STARRS_ID', 'ra', 'dec', 'objInfoFlag', 'qualityFlag', 'rKronRad',
                'gPSFmag','rPSFmag','iPSFmag','zPSFmag','yPSFmag','gPSFmagErr','rPSFmagErr',
                'iPSFmagErr','zPSFmagErr','yPSFmagErr','Pan-STARRS_g','Pan-STARRS_r',
                'Pan-STARRS_i','Pan-STARRS_z','Pan-STARRS_y','Pan-STARRS_g_err','Pan-STARRS_r_err',
                'Pan-STARRS_i_err','Pan-STARRS_z_err','Pan-STARRS_y_err','separation_1',
                'source_id','tmass_key','WISE_W1','WISE_W1_err','WISE_W2','WISE_W2_err',
                'WISE_W3','WISE_W3_err','WISE_W4','WISE_W4_err','SDSS_ID','run',
                'rerun','camcol','SDSS_field','type','SDSS_u','SDSS_g','SDSS_r','SDSS_i',
                'SDSS_z','SDSS_u_err','SDSS_g_err','SDSS_r_err','SDSS_i_err','SDSS_z_err',
                'extinction_u','extinction_g','extinction_r','extinction_i','extinction_z',
                'photo_z','photo_zerr','z_spec','separation_2','DECaL_ID','brick_primary',
                'DECaL_brick','gaia_pointsource','DECaL_g','DECaL_r','DECaL_z','DECaL_g_err',
                'DECaL_r_err','DECaL_z_err','NSC_ID','class_star','NSC_u','NSC_u_err',
                'NSC_g','NSC_g_err','NSC_r','NSC_r_err','NSC_i','NSC_i_err','NSC_z',
                'NSC_z_err','NSC_Y','NSC_Y_err','NSC_VR','NSC_VR_err']
    assert combined_cat.colnames==colnames
    assert combined_cat['Pan-STARRS_ID'][1] == -999.