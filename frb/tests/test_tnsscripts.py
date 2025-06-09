# Module to run tests on instantiating FRB objects
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest

from frb.scripts import tns

remote_data = pytest.mark.skipif(os.getenv('TNS_API_KEY') is None,
                                 reason='test requires dev suite')

remote_data = pytest.mark.skipif(os.getenv('TNS_BOT_ID') is None,
                                 reason='test requires dev suite')

remote_data = pytest.mark.skipif(os.getenv('TNS_BOT_NAME') is None,
                                 reason='test requires dev suite')

@remote_data
def test_tns():
    tns.check_tns_api_keywords()

    radius = 0.001
    theta = 0.0
    dec = 54.3116510708
    ra = 210.910674637
    skip = None
    single_obj = True
    outfile = 'outfile.txt'

    match = tns.main(filename=None,
                     name=None, 
                     ra=ra, 
                     dec=dec,
                     radius=radius, 
                     a=radius,
                     b=radius,
                     theta=theta,
                     single_obj=single_obj)

    assert 129698 in match.keys()
    assert match[129698]['objname']=='2023ixf'
    assert match[129698]['prefix']=='SN'
    assert np.abs(match[129698]['radeg']-ra)<radius
    assert np.abs(match[129698]['decdeg']-dec)<radius
    assert match[129698]['redshift']<0.01

if __name__=="__main__":
    test_tns()
