# Module to run tests on surveys
#  Most of these are *not* done with Travis yet
# TEST_UNICODE_LITERALS

import pytest
import os
import shutil
import numpy as np

from astropy.table import Table

from frb.galaxies.frbgalaxy import FRBHost
from frb.galaxies import eazy as frbeazy
from frb.frb import FRB

from distutils.spawn import find_executable

eazy_exec = pytest.mark.skipif(find_executable('eazy') is None,
                                 reason='test requires galfit')

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def host_obj():
    # VLT
    photom = Table()
    photom['Name'] = ['G_TEST']
    photom['ra'] = 123.422
    photom['dec'] = 23.222
    # These are observed
    photom['LRISb_V'] = 25.86
    photom['LRISb_V_err'] = 0.25
    photom['GMOS_S_r'] = 23.61
    photom['GMOS_S_r_err'] = 0.15
    photom['LRISr_I'] = 23.09
    photom['LRISr_I_err'] = 0.1
    photom['NOT_z'] = 23.35
    photom['NOT_z_err'] = 0.3
    photom['NIRI_J'] = 21.75 + 0.91
    photom['NIRI_J_err'] = 0.2

    #

    host190613A = FRBHost(photom['ra'], photom['dec'], FRB.by_name('FRB20121102'))
    host190613A.parse_photom(photom)
    host190613A.name = 'G_TEST'

    return host190613A

@eazy_exec
def test_eazy(host_obj):
    if os.path.isdir(data_path('eazy')):
        shutil.rmtree(data_path('eazy'))
    os.mkdir(data_path('eazy'))
    # Generate
    frbeazy.eazy_input_files(host_obj.photom, data_path('eazy/input'),
                             host_obj.name,
                             data_path('eazy/output'),
                             templates='br07_default',
                             prior_filter='GMOS_S_r')
    # Test
    assert os.path.isfile(data_path('eazy/input/G_TEST.cat'))
    assert os.path.isfile(data_path('eazy/input/zphot.param.G_TEST'))
    assert os.path.isfile(data_path('eazy/input/zphot.translate.G_TEST'))

    # Run
    frbeazy.run_eazy(data_path('eazy/input'),
                     host_obj.name,
                     os.path.join(data_path('eazy/output'), 'logfile'))
    assert os.path.isfile(data_path('eazy/output/photz.zout'))

    # Read
    zgrid, pzi, prior = frbeazy.getEazyPz(-1, MAIN_OUTPUT_FILE='photz',
                                          OUTPUT_DIRECTORY=data_path('eazy/output'),
                                          CACHE_FILE='Same', binaries=None, get_prior=True)
    
    zphot, sig_zphot = frbeazy.eazy_stats(zgrid, pzi)
    assert np.isclose(zphot, 0.5929259648750858, rtol=1e-4)

    # Remove
    shutil.rmtree(data_path('eazy'))
