# Module to run tests on instantiating FRBObs
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os

from astropy.table import Table

from frb.frbcat import FRBCat

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

def test_init():
    frbobs = FRBCat()
    assert isinstance(frbobs.frbcat, Table)
    assert isinstance(frbobs.uniq_frb, Table)
    # Specify file
    frbobs2 = FRBCat(frbcat_file='frbcat_2017-03-16.csv')
    assert len(frbobs2.frbcat) == 53
    assert len(frbobs2.uniq_frb) == 18

