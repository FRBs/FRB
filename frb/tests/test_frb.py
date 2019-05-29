# Module to run tests on instantiating FRB objects
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os

from astropy.table import Table
from astropy import units

from frb.frb import GenericFRB, FRB

def test_generic():
    FRB1 = GenericFRB(0.6 * units.Jy, 350 * units.MHz, 500 * units.pc / units.cm ** 3)
    # Set width
    FRB1.set_width('Wi', 3e-5 * units.s)
    assert FRB1.Wi.unit == units.s

def test_named():
    frb121102 = FRB('FRB121102', 'J053158.7+330852.5',
                         558.1 * units.pc / units.cm ** 3, z_frb=0.19273)
    # Error ellipse
    frb121102.set_ee(0.1, 0.1, 0., 95.)
    assert isinstance(frb121102.eellipse, dict)
    # Test writing
    frb121102.write_to_json()

    # Test load
    tst = FRB.from_json('FRB121102.json')

    # By name
    tst = FRB.by_name('FRB121102')