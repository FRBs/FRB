# TESTS FOR FINDER CHARTS
import pytest

import numpy as np
import os

import matplotlib
import shutil

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.wcs import WCS
from astropy.io import fits

from frb.figures import finder

remote_data = pytest.mark.remote_data

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_basic():
    if not shutil.which('latex'):
        pass
        return
    # Requires a local latex installation which travis doesn't have..
    # Load up an image
    hdul = fits.open(data_path('DES_r.fits'))
    header = hdul[0].header
    image = hdul[0].data
    wcs = WCS(header)

    # Make it
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    fig, ax = finder.generate(image, wcs, 'FRB 180924', primary_coord=coord,
                              vmnx=(-10., 200.), outfile=data_path('tst.png'))
    # Test
    assert isinstance(fig, matplotlib.figure.Figure)

    '''
    # Log -- Gives a memory crash!
    fig, ax = finder.generate(image, wcs, 'FRB 180924', log_stretch=True, primary_coord=coord,
                              vmnx=(-10., 200.), outfile=data_path('tst.png'))
    '''
    # Cutout
    fig, ax = finder.generate(image, wcs, 'FRB 180924', cutout=(coord, 20*units.arcsec),
                              primary_coord=coord,
                              vmnx=(-10., 200.), outfile=data_path('tst.png'))

