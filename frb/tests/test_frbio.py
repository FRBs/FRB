# Module to run tests on FRB I/O
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os

from astropy.table import Table

from frb.io import load_dla_fits

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

def test_load_dlafits():
    dla_fits = load_dla_fits()
    assert isinstance(dla_fits, dict)
    for key in ['lz', 'fN', 'nenH']:
        assert key in dla_fits.keys()

