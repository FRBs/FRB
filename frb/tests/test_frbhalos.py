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
from astropy.tests.helper import assert_quantity_allclose

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

def test_halo_logM_depr_handler():
    # Confirm that halo init correctly identifies when it's given a log_Mhalo and
    # issues the appropriate warning.
    with pytest.warns(UserWarning, match="Since log_Mhalo is defined"):
        hal0 = halos.Halo(log_Mhalo=12.2)
    with pytest.warns(UserWarning, match="Since log_Mhalo is defined"):
        with pytest.warns(UserWarning, match="Interpreting first argument"):
            hal1 = halos.Halo(12.2)
    with pytest.warns(UserWarning, match="Since log_Mhalo is defined"):
        hal2 = halos.Halo(log_Mhalo=12.2, z=3.0)
    assert hal2.redshift == 3.0
    assert hal2.m_vir == hal1.m_vir


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

def test_dm_method():
    # Error case in parent class, successful call in ModifiedNFW
    hal = halos.Halo(1e14 * un.Msun, 1.3)
    with pytest.raises(NotImplementedError, match="ne function is not defined"):
        hal.dm(0, 15)
    hal = halos.NewModifiedNFW(1e14 * un.Msun, 1.3)
    hal.dm(15)

def test_dm_interp():
    # Test interpolation function and caching
    hal = halos.NewModifiedNFW(1e14 * un.Msun, 1.3)
    xvals = np.linspace(-3, 3, 15)
    dm_vals_05 = hal.dm_interp(xvals, 0.5 * hal.r_vir, xmax=3)
    assert hal._dm_vs_x is not None     # Confirm the interpolator was cached
    assert hal._b == 0.5                # And the impact parameter
    old_int = id(hal._dm_vs_x)
    dm_vals_06 = hal.dm_interp(xvals * hal.r_vir, 0.6 * hal.r_vir)
    assert hal._b == 0.6
    assert id(hal._dm_vs_x) != old_int  # Confirm the interpolator is replaced

    # Check interpolated against quad_vec integrated values
    rv = hal.r_vir.to_value("kpc")
    dm_comp_05 = np.zeros(dm_vals_05.size) * dm_vals_05.unit
    for ii in range(len(dm_comp_05)):
        dm_comp_05[ii] = hal.dm(0.5 * hal.r_vir, smin=-3, smax=xvals[ii])

    assert_quantity_allclose(dm_vals_05, dm_comp_05, atol=1e-3 * un.pc/un.cm**3)

def test_mnfw_mb_of_r():
    # Function computing baryon mass vs radius
    hal = halos.NewModifiedNFW(1e14 * un.Msun, 0.4)
    assert np.isclose(hal.mass_r(1, length_unit='virial'), hal.M_b, rtol=1e-3)


def test_against_old_mnfw():
    # Compare ne, dm, and RM values in NewModifiedNFW to original ModifiedNFW
    # TODO -- also call deprecated functions in NewModifiedNFW
    z = 0
    hal1 = halos.NewModifiedNFW(1e14 * un.Msun, z=0, r_max=1)
    hal0 = halos.ModifiedNFW(np.log10(hal1.m_vir.value), alpha=2, y0=2, c=hal1.conc, f_hot=1.0, z=z)
    radii = np.linspace(0.1, 2.3, 5) * hal1.r_vir
    xyz = np.zeros((3, 5)) * un.kpc
    xyz[0,:] = radii
    ne0 = hal0.ne(xyz.to_value('kpc')) / un.cm**3
    ne1 = hal1.ne(radii)
    assert_quantity_allclose(ne0, ne1)

    impact = 0.5 * hal1.r_vir
    dm0 = hal0.Ne_Rperp(impact, rmax=1)
    dm1 = hal1.dm(impact)
    assert np.isclose(dm0, dm1, atol=0.01 * un.pc/un.cm**3)

    Bpar = 0.5 * un.microgauss
    rm0 = hal0.RM_Rperp(impact, Bpar, rmax=1)
    rm1 = hal1.rm(impact, Bpar)
    assert np.isclose(rm0, rm1, atol = 0.01 * un.rad/un.m**2)


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
    mNFW = halos.NewModifiedNFW()
    # xyz
    xyz = (1-2*rand(3, 100)) * 100
    # rho
    rho = mNFW.rho_b(xyz)
    assert rho.size == 100
    assert rho.unit == un.g/un.cm**3
    # nH
    nH = mNFW.nH(xyz)
    assert rho.size == 100
    xyz0 = [[100.], [0.], [0.]]
    nH0 = mNFW.nH(xyz0).to_value("1/cm3")
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
    coord = SkyCoord('J004244.3+413009', unit=(un.hourangle, un.deg))
    DM = M31.DM_from_Galactic(coord)
    assert DM.unit == un.pc/un.cm**3
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
