# Module to run tests on builds

# TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest

from frb.scripts import build
from frb.builds import build_hosts
from frb.galaxies import frbgalaxy
from frb.frb import FRB

remote_data = pytest.mark.skipif(os.getenv('FRB_GDB') is None,
                                 reason='test requires dev suite')

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@remote_data
def test_host_build():
    outfile = data_path('FRB20180924_host.json')
    if os.path.isfile(outfile):
        os.remove(outfile)

    # Requires a file on disk that is too slow to generate in CI
    pargs = build.parser(['Hosts', '--frb', 'FRB20180924'])
    frbs = pargs.frb.split(',')
    frbs = [ifrb.strip() for ifrb in frbs]

    out_path = data_path('')
    build_hosts.main(frbs, options=pargs.options, 
                             hosts_file=pargs.data_file,
                             lit_refs=pargs.lit_refs,
                             override=pargs.override,
                             out_path=out_path) 

    # Check
    frb20180924 = FRB.by_name('FRB20180924')
    host = frbgalaxy.FRBHost.from_json(frb20180924, outfile)

    # Clean up
    os.remove(outfile)
