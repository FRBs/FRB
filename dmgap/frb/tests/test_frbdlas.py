# Module to run tests on FRB calculations using DLAs
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import pytest
import os

from astropy.units import Quantity
from astropy import units as u

from frb.dlas import approx_avgdm
from frb.dlas import monte_dm
from frb.dlas import monte_tau


def test_approx_avgdm():
    dm = approx_avgdm(1.)
    assert isinstance(dm, Quantity)
    assert dm.unit == (u.pc/u.cm**3)
    assert np.isclose(dm.value, 0.00651401)
    # Array
    dms = approx_avgdm(np.array([0., 1., 2.]))
    assert len(dms) == 3
    # Error
    with pytest.raises(IOError):
        approx_avgdm(10.)


def test_monte_dm():
    """ Monte-carlo of dm values
    """
    zeval = np.array([0.,1.,2.])
    dms = monte_dm(np.array(zeval))
    assert dms.shape[0] == 100
    assert dms.shape[1] == zeval.size


def test_monte_tau():
    """ Monte-carlo of temporal broadening
    """
    zeval = np.array([0.,1.,2.])
    taus = monte_tau(np.array(zeval))
    assert taus.shape[1] == len(zeval)
