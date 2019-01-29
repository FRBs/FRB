# Module to run tests on surveys
#  Most of these are *not* done with Travis yet
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units

from frb.surveys import survey_utils

def test_psrcat():
    try:
        import pulsars
    except ImportError:
        assert True
        return
    # Catalog
    coord = SkyCoord('J000604.8+183459', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    psrcat_srvy = survey_utils.load_survey_by_name('PSRCAT', coord, search_r)
    pulsars = psrcat_srvy.get_catalog()
    #
    assert isinstance(pulsars, Table)
    assert len(pulsars) == 1


def test_des():
    try:
        from dl import queryClient as qc, authClient as ac, helpers
    except ImportError:
        assert True
        return
    # Catalog
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    des_srvy = survey_utils.load_survey_by_name('DES', coord, search_r)
    des_tbl = des_srvy.get_catalog(print_query=True)
    #
    assert isinstance(des_tbl, Table)
    assert len(des_tbl) == 1


def test_decals():
    try:
        from dl import queryClient as qc, authClient as ac, helpers
    except ImportError:
        assert True
        return
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec

    decal_srvy = survey_utils.load_survey_by_name('DECaL', coord, search_r)
    decal_tbl = decal_srvy.get_catalog(print_query=True)
    #
    assert isinstance(decal_tbl, Table)
    assert len(decal_tbl) == 2


def test_first():
    try:
        from dl import queryClient as qc, authClient as ac, helpers
    except ImportError:
        assert True
        return
    coord = SkyCoord('J081240.68+320809', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec
    #
    first_srvy = survey_utils.load_survey_by_name('FIRST', coord, search_r)
    first_tbl = first_srvy.get_catalog()
    #
    assert isinstance(first_tbl, Table)
    assert len(first_tbl) == 1

