# Module to run tests on FRB calculations using DLAs
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import pytest
import os

from astropy.units import Quantity
from astropy import units as u

from frb.dlas import approx_avgDM
from frb.dlas import monte_DM

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

def test_approx_avgDM():
    DM = approx_avgDM(1.)
    assert isinstance(DM, Quantity)
    assert DM.unit == (u.pc/u.cm**3)
    assert np.isclose(DM.value, 0.00651401)
    # Array
    DMs = approx_avgDM(np.array([0., 1., 2.]))
    assert len(DMs) == 3
    # Error
    with pytest.raises(IOError):
        approx_avgDM(10.)

def test_monte_DM():
    zeval = np.array([0.,1.,2.])
    DMs = monte_DM(np.array(zeval))
    assert DMs.shape[0] == 100
    assert DMs.shape[1] == zeval.size
