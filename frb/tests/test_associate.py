import os

import pandas

from astropath import priors
from frb.associate import frbassociate
from frb.associate import frbs

import pytest

remote_data = pytest.mark.skipif(os.getenv('FRB_GDB') is None,
                                 reason='test requires FRB data')

@remote_data
def test_individual():
    # This needs to be hidden
    orig_priors = priors.load_std_priors()
    config = getattr(frbs, 'FRB20180924B')
    frbA = frbassociate.run_individual(config, orig_priors['adopted'])

    # Test
    assert isinstance(frbA.candidates, pandas.DataFrame)
    assert frbA.candidates.iloc[0].P_Ox > 0.98

