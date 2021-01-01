# Module to test 
# the galfit wrapper

import pytest
import os
import shutil
import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from pkg_resources import resource_filename

from frb.galaxies import galfit as glf
remote_data = pytest.mark.skipif(os.getenv('FRB_GDB') is None,
                                        reason='test requires dev suite')

def test_platescale():
    cutout_file = resource_filename('frb','tests/files/cutout_DES_i.fits')
    _, hdr = fits.getdata(cutout_file, header=True)
    wcs = WCS(hdr)
    platescale = glf.get_platescale(wcs)
    assert np.isclose(platescale, 0.263)

@remote_data
def test_run():
    cutout_file = resource_filename('frb','tests/files/cutout_DES_i.fits')
    psf_file = resource_filename('frb', 'tests/files/avg_DES_psf_i.fits')
    badpix = resource_filename('frb', 'tests/files/badpix.fits')
    outdir = resource_filename('frb', 'tests/files/galfit_out')
    return_val = glf.run(cutout_file, psf_file, outdir=outdir, badpix=badpix, r_e=4, n=2, pa=0, finesample=4)
    assert return_val==0
    assert os.path.isdir(outdir)
    outfile = os.path.join(outdir, 'out.fits')
    assert os.path.isfile(outfile)
    hdulist = fits.open(outfile)
    assert len(hdulist)==5
    result_tab = Table(hdulist[4].data)
    assert len(result_tab)==1
    assert np.isclose(result_tab['reff_ang'][0], 0.57071,rtol=1e-2, atol=1e-4)
    shutil.rmtree(outdir)

