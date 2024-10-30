# Module to test 
# the galfit wrapper

import pytest
import os
import shutil
import numpy as np
import shutil

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

import importlib_resources

from frb.frb import FRB

from frb.galaxies import galfit as glf
galfit_exec = pytest.mark.skipif(shutil.which('galfit') is None,
                                        reason='test requires galfit')

def test_platescale():
    cutout_file = importlib_resources.files('frb.tests.files')/'cutout_DES_i.fits'
    _, hdr = fits.getdata(cutout_file, header=True)
    wcs = WCS(hdr)
    platescale = glf.get_platescale(wcs)
    assert np.isclose(platescale, 0.263)

@galfit_exec
def test_run():
    cutout_file = str(importlib_resources.files('frb.tests.files')/'cutout_DES_i.fits')
    psf_file = str(importlib_resources.files('frb.tests.files')/'avg_DES_psf_i.fits')
    badpix = str(importlib_resources.files('frb.tests.files')/'badpix.fits')
    outdir = str(importlib_resources.files('frb.tests.files')/'galfit_out')
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

def test_parse_galfit():
    frb = FRB.by_name("FRB20121102A")
    host = frb.grab_host()
    galfit_outfile = importlib_resources.files('frb.tests.files')/'HG121102_galfit.fits'
    # Test two components
    host.parse_galfit(galfit_outfile,twocomponent=True)
    assert type(host.morphology['PA'])==np.ndarray
    assert len(host.morphology['PA'])==2
    # Test a single component
    host.parse_galfit(galfit_outfile, twocomponent=False)
    assert type(host.morphology['PA'])==np.float64