# Module to run tests on instantiating FRB objects
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np

from astropy import units

from frb.frb import GenericFRB, FRB

def test_generic():
    FRB1 = GenericFRB(0.6 * units.Jy, 350 * units.MHz, 500 * units.pc / units.cm ** 3)
    assert np.isclose(FRB1.DM.value, 500.)

def test_named():
    frb121102 = FRB('FRB121102', 'J053158.7+330852.5',
                         558.1 * units.pc / units.cm ** 3, z_frb=0.19273)
    # Error ellipse
    frb121102.set_ee(0.1, 0.1, 0., 95.)
    assert isinstance(frb121102.eellipse, dict)

    # Pulse -- These are made up
    frb121102.set_pulse(1*units.GHz,
        time_res=0.054*units.ms,
        t0=0.66*units.ms,
        Wi=1.1*units.ms,
        tscatt=0.041*units.ms,
        tscatt_err= 0.002 * units.ms,
        scatt_index=-3.84,
        scatt_index_err=0.77)

    assert np.isclose(frb121102.pulse['freq'].value, 1.)

    # Test writing
    frb121102.write_to_json()
    assert np.isclose(frb121102.pulse['freq'].value, 1.)

    # Test load
    tst = FRB.from_json('FRB121102.json')

    assert np.isclose(tst.pulse['freq'].value, 1.)
