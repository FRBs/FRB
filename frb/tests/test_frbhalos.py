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

from astropy import units as un
from astropy.coordinates import SkyCoord
from astropy import cosmology as acosmo

from frb.halos import models as halos
from frb.halos import hmf

dummy_xyz = np.reshape(np.array([10., 10., 10.]), (3,1))

@pytest.mark.skipif(flg_aemHMF, reason="hmf emulator not installed")
def test_frac_in_halos():
    # Init (this is a test, i.e. it is not used)
    hmf.init_hmf()
    # Run an example
    zvals = np.array([0., 0.1, 0.5, 1.])
    ratios = hmf.frac_in_halos(zvals, 1e11, 1e15)
    # Test
    np.testing.assert_allclose(
        ratios, 
        np.array([0.463629, 0.444242, 0.368071, 0.282909]),rtol=1e-4)

@pytest.mark.skipif(flg_aemHMF, reason="hmf emulator not installed")
def test_halo_incidence():
    # Run
    Navg = hmf.halo_incidence(1e12, 1., Mhigh=3e12)
    assert np.isclose(Navg, 1.1120959929493153)
    # Cumulative now
    zeval, Ncumul = hmf.halo_incidence(1e13, 1., Mhigh=3e13, cumul=True)
    assert zeval.size == 20
    assert np.isclose(Ncumul[-1], 0.5612214870531754, rtol=1e-4)

def test_YF17():
    yf17 = halos.YF17()
    ne = yf17.ne((0.,0.,20.))
    assert np.isclose(ne, 0.000881,atol=1e-6)

# New tests

@pytest.mark.parametrize("cosmo_name", acosmo.available)
def test_base_halo(cosmo_name):
    cosmo = getattr(acosmo, cosmo_name)
    Mv = 1e13 * un.Msun
    hal = halos.Halo(
        Mv,  # Mass
        1.3, # z
        cosmo=cosmo
    )
    assert np.isclose(hal.M_b, Mv * cosmo.Ob0 / cosmo.Om0)
    assert hal.dist_comov == cosmo.comoving_distance(1.3)

def test_nfw_virial_conventions():
    Mvir = 1e12 * un.Msun
    hal0 = halos.NFWHalo(Mvir, 0, conc=3.5)
    assert hal0.conc == 3.5
    assert np.isclose(hal0.r_vir, ((Mvir / (4/3 * np.pi * hal0.rho_vir))**(1/3)).to('kpc'))
    M200, c200 = hal0.mc_convert(200)
    hal1 = halos.NFWHalo(M200, hal0.redshift, conc=c200, del_c=200)
    with pytest.warns(UserWarning, match="Overdensity factor is"):
        assert np.isclose(hal0.r200, hal1.r200)
    assert np.isclose(hal1.r200, hal1.r_vir)


def test_mnfw():
    Mvir = 1e12 * un.Msun
    z = 1.3
    hal = halos.NewModifiedNFW(Mvir, z, conc=7.3, y0=2, alpha=2)
    print(hal.y0, hal.alpha, hal.r_vir)
    assert hal.conc == 7.3
    assert hal.y0 == 2


def test_mnfw_deprecated_calls():
    # Check that deprecated function calls yield warning
    Mvir = 1e12 * un.Msun
    z = 1.3
    hal = halos.NewModifiedNFW(Mvir, z, conc=7.3, y0=2, alpha=2)
    xyz = np.zeros((3, 10))
    xyz[0, :] = np.arange(1, 11)
    xyz_q = xyz * un.kpc
    with pytest.warns(UserWarning, match="Passing xyz cartesian coordinates to ne is deprecated"):
        hal.ne(xyz)
    hal.zero_inner_ne = 5.
    with pytest.warns(UserWarning, match="Passing xyz cartesian coordinates to ne is deprecat"):
        with pytest.warns(UserWarning, match="Attribute zero_inner_ne is deprecated"):
            hal.ne(xyz_q)
    with pytest.warns(UserWarning, match="Attribute `c` has been renamed"):
        hal.c == 7.3
    with pytest.warns(UserWarning, match="Attribute `z` has been renamed"):
        hal.z == z

#####
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
    assert rho.unit == g/un.cm**3
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
    coord = SkyCoord('J004244.3+413009', unit=(hourangle, un.deg))
    DM = M31.DM_from_Galactic(coord)
    assert DM.unit == pc/un.cm**3
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
