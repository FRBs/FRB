# Module to run tests on surveys
#  Most of these are *not* done with Travis yet
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units

from frb.galaxies import frbgalaxy, defs

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_frbhost():
    repeater_coord = SkyCoord('05h31m58.698s +33d8m52.59s', frame='icrs')  # Use as host here
    # Instantiate
    host121102 = frbgalaxy.FRBHost(repeater_coord.ra.value, repeater_coord.dec.value, '121102')
    # Add a few nebular lines  (Tengulkar+17)
    neb_lines = {}
    neb_lines['Ha'] = 0.652e-16
    neb_lines['Ha_err'] = 0.009e-16
    neb_lines['Ha_Al'] = 0.622
    #
    neb_lines['Hb'] = 0.118e-16
    neb_lines['Hb_err'] = 0.011e-16
    neb_lines['Hb_Al'] = 0.941
    AV = 2.42
    # Deal with Galactic extinction
    for key in neb_lines.keys():
        if '_err' in key:
            continue
        if 'Al' in key:
            continue
        # Ingest
        host121102.neb_lines[key] = neb_lines[key] * 10 ** (neb_lines[key + '_Al'] * AV / 2.5)
        host121102.neb_lines[key + '_err'] = neb_lines[key + '_err'] * 10 ** (neb_lines[key + '_Al'] * AV / 2.5)
    # Vette
    for key in host121102.neb_lines.keys():
        if '_err' in key:
            continue
        assert key in defs.valid_neb_lines
    # Morphology
    host121102.morphology['reff_ang'] = 0.41
    host121102.morphology['reff_ang_err'] = 0.06
    #
    host121102.morphology['n'] = 2.2
    host121102.morphology['n_err'] = 1.5
    #
    host121102.morphology['b/a'] = 0.25
    host121102.morphology['b/a'] = 0.13
    # Vette
    for key in host121102.morphology.keys():
        if '_err' in key:
            continue
        assert key in defs.valid_morphology
    # Write
    outfile = data_path('test_frbhost.json')
    host121102.write_to_json(outfile=outfile)


def test_read_frbhost():
    # This test will fail if the previous failed
    host121102 = frbgalaxy.FRBHost.from_json(data_path('test_frbhost.json'))
    # Test
    assert host121102.frb == '121102'
    assert np.isclose(host121102.morphology['b/a'], 0.13)
    
def test_by_name():
    host121102 = frbgalaxy.FRBHost.by_name('121102')
    assert host121102.frb == '121102'
    assert np.isclose(host121102.morphology['b/a'], 0.13)
