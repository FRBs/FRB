# Module to run tests on em.py module
#  Most of these are *not* done with Travis yet
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units

from frb import em


def test_em():
    # EM from Halpha
    em_121102 = em.em_from_halpha(6.8e-16 * units.erg / units.cm ** 2 / units.s / units.arcsec ** 2, 0.1927)
    # Test
    assert em_121102.unit == units.pc / units.cm**6
    assert np.isclose(em_121102.value, 668.58867698)


def test_dm_from_em():
    em_121102 = em.em_from_halpha(6.8e-16 * units.erg / units.cm ** 2 / units.s / units.arcsec ** 2, 0.1927)
    DM_s = em.dm_from_em(em_121102, 1 * units.kpc)
    # Test
    assert DM_s.unit == units.pc / units.cm**2
    assert np.isclose(DM_s.value, 408.52143)
