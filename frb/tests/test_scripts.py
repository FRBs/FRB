# Module to run tests on instantiating FRB objects
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest

from frb.scripts import frb_summary
from frb.scripts import limiting_mag
from frb.scripts import pz_dm
from frb.scripts import tns


def test_frb_summary():
    pargs = frb_summary.parser(['180924'])
    frb_summary.main(pargs)

def test_frb_mag_limit():
    # Requires a file on disk that is too slow to generate in CI
    pargs = limiting_mag.parser(['J151849.52+122235.8', '200.', '23.'])
    Lmin, Lmax = limiting_mag.main(pargs)

    assert np.isclose(Lmax, 0.018052542432481264)

def test_frb_pz_dm():
    # Requires a file on disk that is too slow to generate in CI
    pargs = pz_dm.parser(['J151849.52+122235.8', '200.'])
    zmin, zmax, z_50, z_mode = pz_dm.main(pargs)

    assert np.isclose(zmin, 0.04020100502512563)
    assert np.isclose(zmax, 0.16080402010050251)
    assert np.isclose(z_50, 0.10050251256281408)
    assert np.isclose(z_mode, 0.12060301507537688)

def test_tns():
    tns.parser.units = 'deg'
    tns.parser.radius = 0.5
    tns.parser.dec = 19.17
    tns.parser.ra = 4.49
    tns.parser.skip = None
    tns.parser.outfile = 'outfile.txt'
    
    df = tns.main(tns.parser)

    assert df['prefix'][0] == 'FRB'
