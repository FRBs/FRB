# Module to run tests on instantiating FRB objects
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest

from frb.scripts import frb_summary


def test_frb_summary():
    pargs = frb_summary.parser(['180924'])
    frb_summary.main(pargs)
