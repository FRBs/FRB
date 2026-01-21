# Module to run tests on instantiating FRB objects
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest

from frb.scripts import frb_summary
from frb.scripts import pzdm_mag

remote_data = pytest.mark.skipif(os.getenv('FRB_GDB') is None,
                                 reason='test requires dev suite')

def test_frb_summary():
    pargs = frb_summary.parser(['180924'])
    frb_summary.main(pargs)

def test_frb_pzdm_mag():
    # Requires a file on disk that is too slow to generate in CI
    pargs = pzdm_mag.parser(['J151849.52+122235.8', '200.','--mag_limit', '23.'])
    zmin, zmax, z_50, z_mode, Lmin, Lmax = pzdm_mag.main(pargs)

    assert np.isclose(zmin, 0.04020100502512563)
    assert np.isclose(zmax, 0.16080402010050251)
    assert np.isclose(z_50, 0.10050251256281408)
    assert np.isclose(z_mode, 0.12060301507537688)
    assert np.isclose(Lmax, 0.023809, atol=1e-4)
