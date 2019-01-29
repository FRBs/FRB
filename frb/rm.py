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
    galactic_rm_file = resource_filename('frb', 'data/RM/opp14_foreground.fits')

    # Load
    rm_sky = hp.read_map(galactic_rm_file, hdu=4)
    sig_sky = hp.read_map(galactic_rm_file, hdu=6)
    nside = 128

    #
    pix = hp.ang2pix(nside, coord.galactic.l.to('rad').value,
                     coord.galactic.b.to('rad').value)
    # Return
    return rm_sky[pix], sig_sky[pix]

