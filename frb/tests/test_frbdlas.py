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
from frb.dlas import monte_tau


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
    """ Monte-carlo of DM values
    """
    zeval = np.array([0.,1.,2.])
    DMs = monte_DM(np.array(zeval))
    assert DMs.shape[0] == 100
    assert DMs.shape[1] == zeval.size


def test_monte_tau():
    """ Monte-carlo of temporal broadening
    """
    zeval = np.array([0.5,1.,2.])
    taus = monte_tau(np.array(zeval))
    assert taus.shape[1] == len(zeval)
