# Module to run tests on FRB halo calculations

# TEST_UNICODE_LITERALS

import numpy as np
import pytest
from numpy.random import rand

try:
    import hmf_emulator
except ImportError:
    flg_aemHMF = False
else:
    flg_aemHMF = True

from astropy import units as u
from astropy.coordinates import SkyCoord

from frb.halos import models as halos
from frb.halos import hmf

dummy_xyz = np.reshape(np.array([10., 10., 10.]), (3,1))

def test_modified_NFW():
    # Init
    mNFW = halos.ModifiedNFW()
    # xyz
    xyz = (1-2*rand(3, 100)) * 100
    # rho
    rho = mNFW.rho_b(xyz)
    assert rho.size == 100
    assert rho.unit == u.g/u.cm**3
    # nH
    nH = mNFW.nH(xyz)
    assert rho.size == 100
    xyz0 = [100., 0., 0.]
    nH0 = mNFW.nH(xyz0)
    assert np.isclose(nH0, 0.00018175, rtol=1e-3)
    # ne
    ne = mNFW.ne(xyz)
    assert np.all(ne > nH)

#def test_mNFW_asymptote():
    # Init
mNFW = halos.ModifiedNFW(alpha=2, y0=2, f_hot=1., norm_by_asymptote=True)