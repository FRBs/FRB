""" Module for FRBObs Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import warnings
import os

from pkg_resources import resource_filename

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord



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
        path = resource_filename('frb', os.path.join('data','FRBs'))
        if frbcat_file is None:
            fils = glob.glob(os.path.join(path,'frbcat_*'))
            #fils = glob.glob(frbdm.__path__[0]+'/data/FRBs/frbcat_*')
            fils.sort()
            infil = fils[-1]  # Expecting these are ordered by date
            self.frbcat_file = infil
        else:
            self.frbcat_file = os.path.join(path,frbcat_file)
        # Read
        if 'csv' in self.frbcat_file:
            self.frbcat = Table.read(self.frbcat_file, format='csv')#, delimiter='#')
        else:
            raise IOError("Not prepared for this file format")
        if self.verbose:
            print("Using {:s} for the FRB catalog".format(self.frbcat_file))
        # Muck with RA/DEC
        if 'RAJ' in self.frbcat.keys(): # ORIGINAL
            orig=True
            cname = []
            for row in self.frbcat:
                cname.append('{:s} {:s}'.format(row['RAJ'], row['DECJ']))
            self.coords = SkyCoord(cname, unit=(u.hourangle, u.deg))
        elif 'rop_gl' in self.frbcat.keys(): # 2018
            orig=False
            self.coords = SkyCoord(l=self.frbcat['rop_gl'], b=self.frbcat['rop_gb'], unit='deg', frame='galactic')
        # Set
        self.frbcat['RA'] = self.coords.icrs.ra.value
        self.frbcat['DEC'] = self.coords.icrs.dec.value

        # Restrict to unique sources
        if orig:
            uni, uidx = np.unique(self.frbcat['Name'], return_index=True)
        else:
            # Hack for CVS problem
            keys = list(self.frbcat.keys())
            if 'frb_name' in keys[0]:
                self.frbcat.rename_column(keys[0], 'frb_name')
            uni, uidx = np.unique(self.frbcat['frb_name'], return_index=True)
            # Fix UTC
            for row in self.frbcat:
                row['utc'] = row['utc'].replace('/','-')
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
