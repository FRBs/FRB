""" Module for RM calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import os

from pkg_resources import resource_filename

from astropy import units

import healpy as hp

def galactic_rm(coord):
    """
    Provide the Oppermann et al. 2014 -- https://arxiv.org/abs/1404.3701
    estimate for Galactic Farady RM and its uncertainty towards the input coordinate

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for the RM esimation

    Returns:
        Quantity, Quantity: RM and RM_err with units of rad/m^2

    """
    galactic_rm_file = resource_filename('frb', 'data/RM/opp14_foreground.fits')

    # Load
    rm_sky = hp.read_map(galactic_rm_file, hdu=4)
    sig_sky = hp.read_map(galactic_rm_file, hdu=6)
    nside = 128

    # Find the pixel
    pix = hp.ang2pix(nside, coord.galactic.l.value, coord.galactic.b.value, lonlat=True)

    # Return
    return rm_sky[pix]*units.rad/units.m**2, sig_sky[pix]*units.rad/units.m**2

