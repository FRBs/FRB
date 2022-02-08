# Module to run tests on FRB halo calculations
from __future__ import print_function, absolute_import, division, unicode_literals

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

def test_frac_in_halos():
    # Imported (unlikely)?
    if not flg_aemHMF:
        assert True
        return
    # Init (this is a test, i.e. it is not used)
    hmf.init_hmf()
    # Run an example
    zvals = np.array([0., 0.1, 0.5, 1.])
    ratios = hmf.frac_in_halos(zvals, 1e11, 1e15)
    # Test
    np.testing.assert_allclose(
        ratios, 
        np.array([0.4626 , 0.44315, 0.3669, 0.2818]),rtol=1e-4)

def test_halo_incidence():
    # Imported (unlikely)?
    if not flg_aemHMF:
        assert True
        return
    # Run
    Navg = hmf.halo_incidence(1e12, 1., Mhigh=3e12)
    assert np.isclose(Navg, 1.1079098177338418)
    # Cumulative now
    zeval, Ncumul = hmf.halo_incidence(1e13, 1., Mhigh=3e13, cumul=True)
    assert zeval.size == 20
    assert np.isclose(Ncumul[-1], 0.5578253455474869, rtol=1e-4)

def test_YF17():
    yf17 = halos.YF17()
    ne = yf17.ne((0.,0.,20.))
    assert np.isclose(ne, 0.000881,atol=1e-6)

def test_MB04():
    mb04 = halos.MB04()
    ne = mb04.ne((0.,0.,20.))
    # Test
    assert np.isclose(ne, 0.0003937, rtol=1e-3)

def test_MB15():
    mb15 = halos.MB15()
    ne = mb15.ne((0.,0.,20.))
    # Test
    assert np.isclose(ne, 0.00016150865297256291)


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

def test_milky_way():
    Galaxy = halos.MilkyWay()
    assert np.isclose(Galaxy.M_halo.to('M_sun').value, 1.51356125e+12)

def test_m31():
    M31 = halos.M31()
    assert np.isclose(M31.coord.distance.to('kpc').value, 752.)
    # DM through the halo
    coord = SkyCoord('J004244.3+413009', unit=(u.hourangle, u.deg))
    DM = M31.DM_from_Galactic(coord)
    assert DM.unit == u.pc/u.cm**3
    assert np.isclose(DM.value, 80.606, rtol=1e-3)

def test_satellites():
    smc = halos.SMC()
    lmc = halos.LMC()
    m33 = halos.M33()

def test_ICM():
    icm = halos.ICM(log_Mhalo=14.5)
    ne = icm.ne(dummy_xyz)
    #
    assert np.isclose(ne, 0.012977, rtol=1e-3)
