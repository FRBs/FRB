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
    assert len(des_tbl) == 1


@remote_data
def test_decals():
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    decal_srvy = survey_utils.load_survey_by_name('DECaL', coord, search_r)
    decal_tbl = decal_srvy.get_catalog(print_query=True)
    #
    assert isinstance(decal_tbl, Table)
    assert len(decal_tbl) == 2



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
    assert len(ps_table) == 25

    #Test get_cutout
    cutout, = ps_survey.get_cutout()
    assert isinstance(cutout,Image.Image)
    assert cutout.size == (120,120)

    #Test get_image
    imghdu = ps_survey.get_image()
    assert isinstance(imghdu,PrimaryHDU)
    assert imghdu.data.shape == (120,120)

def test_in_which_survey():
    """
    To test if `survey_utils.in_which_survey` works.
    """
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    
    with warnings.catch_warnings(record=True) as allwarns:
        inside = survey_utils.in_which_survey(coord)

    expected_dict = {'SDSS': True,
                     'DES': False,
                     'NVSS': False,
                     'FIRST': False,
                     'WENSS': False,
                     'DECaL': True,
                     'WISE': True,
                     'Pan-STARRS': True}

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