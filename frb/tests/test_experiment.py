# Module to run tests on instantiating FRBObs
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os

from astropy.table import Table

from frb.experiment import Experiment

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

def test_init():
    chime = Experiment('chime')
    for key in ['FOV', 'G', 'np', 'Trec', 'Channels']:
        assert key in chime.data.keys()

