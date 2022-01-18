# Module to run tests on turbulent scattering
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import pytest

from astropy import units as u
from astropy import constants as const

from frb.turb_scattering import Turbulence

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

# Default values
def_l0 = 1*u.AU
def_L0 = 0.001 * u.pc
def_ne = 1e-2 / u.cm**3
def_DL = 1 * u.kpc
def_zL = 1.


def test_init_turb():
    #
    turb = Turbulence(def_ne, def_l0, def_L0, def_zL)
    print(turb)
    assert turb.regime == 0
    # Try another beta
    turb2 = Turbulence(def_ne, def_l0, def_L0, def_zL, beta=3.)
    assert np.isclose(turb2.beta, 3.)


def test_temp_smearing():
    inu = 1 * u.GHz
    lobs = (const.c/inu).to('cm')
    # Galaxy at z=1
    turb = Turbulence(def_ne, def_l0, def_L0, def_zL, DL=def_DL, lobs=lobs)
    # tau
    zsource = 2.
    tau = turb.temporal_smearing(lobs, zsource)
    # Test
    assert tau.unit == u.ms
    assert np.isclose(tau.to('ms').value, 0.48865, rtol=1e-3)


def test_ang_broadening():
    # Galaxy at z=1
    turb = Turbulence(def_ne, def_l0, def_L0, def_zL)
    inu = 1 * u.GHz
    lobs = (const.c/inu).to('cm')
    # SM, rdiff
    turb.set_SM_obj(def_DL)
    turb.set_rdiff(lobs)
    # tau
    zsource = 2.
    theta = turb.angular_broadening(lobs, zsource)
    assert theta.unit == u.arcsec
    assert np.isclose(theta.to('arcsec').value, 9.26810535e-06)


