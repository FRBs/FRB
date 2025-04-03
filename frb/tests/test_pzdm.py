# Module to run tests on surveys
#  Most of these are *not* done with Travis yet
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

#import pytest
import os
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table

from linetools.spectra import xspectrum1d

from frb.galaxies import frbgalaxy, defs, utils
from frb.frb import FRB

from frb.galaxies import mag_dm
from frb.dm import prob_dmz


def test_pzdm_telescopes():
    
    telescope_dict = prob_dmz.telescope_dict
    #telescope_dict = {
    #    'DSA': 'DSA_pzdm.npy',
    #    'Parkes': 'parkes_mb_class_I_and_II_pzdm.npy',
    #    'CRAFT': 'CRAFT_class_I_and_II_pzdm.npy',
    #    'CRAFT_ICS_1300': 'CRAFT_ICS_1300_pzdm.npy',
    #    'CRAFT_ICS_892': 'CRAFT_ICS_892_pzdm.npy',
    #    'CRAFT_ICS_1632': 'CRAFT_ICS_1632_pzdm.npy',
    #    'FAST': 'FAST_pzdm.npy',
    #}

    # Load the CHIME grid
    sdict = prob_dmz.grab_repo_grid(telescope_dict['CHIME'])#'CHIME_pzdm.npz')
    PDM_z = sdict['pzdm']
    z = sdict['z']
    DM = sdict['DM']

    # Test
    assert len(z) == 500
    assert len(DM) == 1400
    assert PDM_z.shape == (500, 1400)
    
    # check the perfect grid
    sdict = prob_dmz.grab_repo_grid(telescope_dict['perfect'])
    PDM_z = sdict['PDM_z']
    z = sdict['z']
    DM = sdict['DM']
    
    # Test      
    assert len(z) == 200    
    assert len(DM) == 1000
    assert PDM_z.shape == (1000, 200)

    # check the other grids
    for key in telescope_dict.keys():
        if key == 'perfect':
            continue
        grid = prob_dmz.grab_repo_grid(telescope_dict[key])
        PDM_z = grid['pzdm']
        assert PDM_z.shape == (500, 1400)

    
    # test full run
    from frb.scripts.pzdm_mag import parser, main
    args = parser(["J081240.7+320809", "500", "--dm_host", "60", "--dm_mwhalo", "40",  "--telescope", "CHIME"])
    z_min, z_max, z_50, z_mode, frac_Lstar_min, frac_Lstar_max = main(args)
    assert isinstance(z_min, float)
    assert isinstance(z_max, float)
    assert isinstance(z_50, float)
    assert isinstance(z_mode, float)

    assert np.isclose(frac_Lstar_min, 0.00804, atol=1e-4)
    assert np.isclose(frac_Lstar_max, 2.88776, atol=1e-4)
    assert z_min < z_max
    assert 0 <= z_50 <= 1
    assert 0 <= z_mode <= 1
