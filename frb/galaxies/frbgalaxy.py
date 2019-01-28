""" Module for galaxies related to FRBs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

import warnings

from pkg_resources import resource_filename

from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.special import hyp2f1
from scipy.interpolate import interp1d

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy import constants
from astropy.table import Table


class FRBGalaxy(object):
    """

    """
    def __init__(self, ra, dec, frb):

        # Init
        self.coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        self.frb = frb

        # Main attributes
        self.redshift = {}
        self.photom = {}
        self.morphology = {}
        self.lines = {}
        self.derived = {}
        self.kinematics = {}

    def parse_photom(self, phot_tbl, max_off=1*units.arcsec, overwrite=True):
        phot_coord = SkyCoord(ra=phot_tbl['ra'], dec=phot_tbl['dec'], unit='deg')
        sep = self.coord.separation(phot_coord)
        imin = np.argmin(sep)
        # Satisfy minimum offset?
        if sep[imin] > max_off:
            print("No photometric sources within {} of the galaxy".format(max_off))
            return
        # Fill


class FRBHost(FRBGalaxy):

    def __init__(self, ra, dec, frb, z_frb=None, **kwargs):
        # Instantiate
        super(FRBHost, self).__init__(ra, dec, frb, **kwargs)

        # Load up FRB info from name

        # Optional
        if z_frb is not None:
            self.redshift['z_FRB'] = z_frb
