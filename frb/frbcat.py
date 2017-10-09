""" Module for FRBObs Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import warnings

from pkg_resources import resource_filename

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord

from pkg_resources import resource_filename

try:
    basestring
except NameError:  # For Python 3
    basestring = str

class FRBCat(object):
    """ Class to load up and provide FRB Observations in a simple package
    Designed to ingest tables from FRBCAT

    Parameters
    ----------

    Attributes
    ----------
    """

    def __init__(self, frbcat_file=None, verbose=True, **kwargs):
        """
        """
        # Init
        self.verbose = verbose
        self.load_cat(frbcat_file=frbcat_file)
        # Name, Creation date
        if verbose:
            pass
            #print("Created on {:s}".format(spdbu.hdf_decode(self.qcat.cat_attr['CREATION_DATE'])))
        # Return
        return

    def load_cat(self, frbcat_file=None):
        """ Load the catalog

        Parameters
        ----------
        db_file : str

        Returns
        -------

        """
        import glob
        path = resource_filename('frb', 'data/FRBs')
        if frbcat_file is None:
            fils = glob.glob(path + '/frbcat_*')
            #fils = glob.glob(frbdm.__path__[0]+'/data/FRBs/frbcat_*')
            fils.sort()
            infil = fils[-1]  # Expecting these are ordered by date
            self.frbcat_file = infil
        else:
            self.frbcat_file = path+frbcat_file
        # Read
        if 'csv' in self.frbcat_file:
            self.frbcat = Table.read(self.frbcat_file, format='csv')
        else:
            raise IOError("Not prepared for this file format")
        if self.verbose:
            print("Using {:s} for the FRB catalog".format(self.frbcat_file))
        # Muck with RA/DEC
        coord_list = []
        cname = []
        for row in self.frbcat:
            cname.append('{:s} {:s}'.format(row['RAJ'], row['DECJ']))
        self.coords = SkyCoord(cname, unit=(u.hourangle, u.deg))
        self.frbcat['RA'] = self.coords.ra.value
        self.frbcat['DEC'] = self.coords.dec.value
        # Restrict to unique sources
        uni, uidx = np.unique(self.frbcat['Name'], return_index=True)
        self.unidx = uidx
        self.uniq_frb = self.frbcat[uidx]
        self.nfrb = len(self.uniq_frb)

    def __repr__(self):
        txt = '<{:s}:  FRB Catalog with {:d} sources\n'.format(self.__class__.__name__,
                                                                 len(self.uniq_frb))
        # Surveys
        #txt += '   Loaded groups are {} \n'.format(self.groups)
        txt += '>'
        return (txt)
