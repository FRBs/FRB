#!/bin/env python3
import os
import os.path
from . import surveycoord
from astropy.table import Table
from frb.defs import frb_cosmo
from astropy.coordinates import SkyCoord
from astropy import units as u

import numpy as np


class NEDLVS(surveycoord.SurveyCoord):
    """
    A class to handle local NED Local Volume Sample queries.
    This requires the LVS table to be downloaded
    from https://ned.ipac.caltech.edu/NED::LVS/
    and linked via the environment variable NEDLVS
    """


    def __init__(self, coord, radius=90.*u.deg, cosmo=None, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        assert 'NEDLVS' in os.environ, "NEDLVS environment variable not set. Please download the LVS table from https://ned.ipac.caltech.edu/NED::LVS/ and set the environment variable NEDLVS to the path of the downloaded file."        
        self.survey = 'NEDLVS'
        self.datapath = os.environ['NEDLVS']
        if cosmo is None:
            self.cosmo = frb_cosmo
        else:
            self.cosmo = cosmo

        # Read in the data and store in memory
        self.datatab = Table.read(self.datapath)
        self.datatab['coord'] = SkyCoord(self.datatab['ra'], self.datatab['dec'], unit="deg")
        self.datatab['ang_sep'] = self.coord.separation(self.datatab['coord']).to('arcmin')

        # Set redshift distances using the cosmology of choice
        redshift_dist_sources = self.datatab['DistMpc_method']=='Redshift'
        self.datatab['DistMpc'][redshift_dist_sources] = self.cosmo.luminosity_distance(self.datatab['z'][redshift_dist_sources]).to('Mpc').value
        self.datatab['phys_sep'] = self.datatab['DistMpc']*u.Mpc*np.sin(self.datatab['ang_sep'].to('rad').value)
    
    def get_column_names(self):
        return self.datatab.colnames

    def get_catalog(self, z_lim=np.inf,
                    impact_par_lim=np.inf*u.Mpc,
                    query_fields=None,
                    print_query=False):
        """
        Get the catalog of objects within the given limits of redshift, impact parameter, and angular separation.
        Args:
            z_lim (float): The maximum redshift of the objects to include in the catalog.
            impact_par_lim (Quantity): The maximum impact parameter of the objects to include in the catalog.
            ang_sep_lim (Quantity): The maximum angular separation of the objects to include in the catalog.
            query_fields (list): The fields to include in the catalog. If None, the default fields are used.
        Returns:
            A table of objects within the given limits.
        """
        if query_fields is None:
            query_fields = ['objname', 'ra', 'dec', 'ebv', 'z', 'z_unc', 'z_tech', 'DistMpc', 'DistMpc_unc', 'DistMpc_method', 'Mstar', 'Mstar_unc', 'ang_sep', 'phys_sep'] 
        else:
            assert np.isin(query_fields, self.datatab.colnames).all(), "One or more of the requested fields is not in the NEDLVS table. Check the column names with get_column_names()."
        # ...
        if print_query:
            print(f"Querying NEDLVS for objects within {self.radius} of {self.coord} with z < {z_lim} and impact parameter < {impact_par_lim}.")
            print(f"Query fields: {query_fields}")
        distance_cut = self.datatab['DistMpc']<self.cosmo.luminosity_distance(z_lim).to('Mpc').value #Only need foreground objects
        valid_distances = self.datatab['DistMpc']>0 # Exclude weird sources with negative distances
        phys_sep_cut = self.datatab['phys_sep']<impact_par_lim # Impact param within limit

         # Make sure the earth is not between the FRB and the galaxy
        assert self.radius<=90*u.deg, "The radius of the search cone is too large. Please set it to a value less than 90 degrees."

        ang_sep_cut = self.datatab['ang_sep']<self.radius
        is_nearby_fg = valid_distances&distance_cut & phys_sep_cut & ang_sep_cut
        
        close_by = self.datatab[is_nearby_fg][query_fields]

        self.catalog = close_by
        return self.catalog