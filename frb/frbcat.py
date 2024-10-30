""" Module for FRBObs Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

import importlib_resources

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
        path = importlib_resources.files('frb.data.FRBs')
        if frbcat_file is None:
            fils = glob.glob(str(path/'frbcat_*'))
            fils.sort()
            infil = fils[-1]  # Expecting these are ordered by date
            self.frbcat_file = infil
        else:
            self.frbcat_file = str(path/frbcat_file)
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

        # DM
        self.frbcat['DM'] = 0.
        self.frbcat['DM_err'] = -999.
        for kk, dms in enumerate(self.frbcat['rmp_dm']):
            ps = dms.split('&plusmn')
            self.frbcat['DM'][kk] = float(ps[0])
            if 'plusmn' in dms:
                self.frbcat['DM_err'][kk] = float(ps[1])

        # RM
        self.frbcat['RM'] = np.nan
        self.frbcat['RM_err'] = np.nan
        for kk, rms in enumerate(self.frbcat['rmp_rm']):
            if rms == 'null':
                continue
            ps = rms.split('&plusmn')
            self.frbcat['RM'][kk] = float(ps[0])
            if 'plusmn' in rms:
                self.frbcat['RM_err'][kk] = float(ps[1])


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
