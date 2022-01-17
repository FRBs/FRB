# Module to run tests on FRB calculations using DLAs
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import pytest

from astropy import units as u
from astropy.cosmology import Planck15, FlatLambdaCDM

from frb.dm import igm

def test_rhoMstar():
    rho_Mstar_full = igm.avg_rhoMstar(1., remnants=True)
    # Test
    assert rho_Mstar_full.unit == u.Msun/u.Mpc**3
    assert np.isclose(rho_Mstar_full.value, 4.65882439e+08)

def test_rhoISM():
    rhoISM = igm.avg_rhoISM(0.)
    # Test
    assert rhoISM.unit == u.Msun/u.Mpc**3
    assert np.isclose(rhoISM.value, 2.19389268e+08)

def test_igmDM():
    DM = igm.average_DM(1.)
    # Value and unit
    assert DM.unit == u.pc/u.cm**3
    assert np.isclose(DM.value, 923.074, rtol=0.001)
    # Cumulative
    DM_cum, _ = igm.average_DM(1., cumul=True)
    assert DM == DM_cum[-1]
    # Cross through HeII reionization
    DM4 = igm.average_DM(4.)
    assert np.isclose(DM4.value, 3542.598, rtol=0.001)


def test_z_from_DM():
    # Note this removes 100 DM units of 'nuisance'
    z = igm.z_from_DM(1000.*u.pc/u.cm**3)
    # Test
    assert np.isclose(z, 0.97557, rtol=0.001)

def test_igmDM_varyH0():
    DM = igm.average_DM(1., cosmo=Planck15)

    # New cosmology
    cosmo2 = FlatLambdaCDM(H0=Planck15.H0.value*2, 
                          Ob0=Planck15.Ob0, 
                          Om0=Planck15.Om0, 
                          Neff=Planck15.Neff,
                          Tcmb0=2.725)
    DM2 = igm.average_DM(1., cosmo=cosmo2)

    # Test
    assert DM2 > DM*2   # Not sure why it is so much higher
