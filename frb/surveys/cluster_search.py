"""
A module to query for galaxy groups/clusters around a given FRB.
Currently has only the Tully cluster catalog but can be possibly extended for
other sources.
"""

from . import surveycoord
from frb.defs import frb_cosmo
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.vizier import Vizier

import numpy as np

class VizierCatalogSearch(surveycoord.SurveyCoord):
    """
    A class to query sources within a Vizier catalog.
    """


    def __init__(self, coord, radius = 90*u.deg, catalog = None, cosmo=None, **kwargs):
        # Initialize a SurveyCoord object
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        self.survey = None # Name
        self.catalog = catalog # Name of the Vizier table to draw from.
        self.coord = coord # Location around which to perform the search
        self.radius = radius.to('deg').value # Radius of cone search
        if cosmo is None: # Use the same cosmology as elsewhere in this repository unless specified.
            self.cosmo = frb_cosmo
        else:
            self.cosmo = cosmo
    

    def get_catalog(self, query_fields=None, transverse_distance_cut = 5*u.Mpc, **kwargs):
        pass

        

class TullyGroupCat(surveycoord.SurveyCoord):
    """
    A class to query sources within the Tully 2015
    group/cluster catalog.
    """


    def __init__(self, coord, radius = 90*u.deg, cosmo=None, **kwargs):
        # Initialize a SurveyCoord object
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        self.survey = 'Tully_2015' # Name
        self.catalog = "J/AJ/149/171/table5" # Name of the Vizier table to draw from.
        self.coord = coord # Location around which to perform the search
        self.radius = radius.to('deg').value # Radius of cone search
        if cosmo is None: # Use the same cosmology as elsewhere in this repository unless specified.
            self.cosmo = frb_cosmo
        else:
            self.cosmo = cosmo
    

    def get_catalog(self, query_fields=None, transverse_distance_cut = np.inf*u.Mpc, richness_cut = 5):
        """
        Get the catalog of objects
        Args:
            z_lim (float): The maximum redshift of the objects to include in the catalog.
            transverse_distance_cut (Quantity): The maximum impact parameter of the objects to include in the catalog.
            richness_cut (int): The minimum number of members in any group/cluster returned.
            query_fields (list): The fields to include in the catalog. If None, all fields are used.
        Returns:
            A table of objects within the given limits.
        """
        if query_fields is None:
            query_fields = ['**'] # Get all.
        # Query Vizier
        v = Vizier(catalog = self.catalog, columns=query_fields, row_limit= -1) # No row limit
        result = v.query_region(self.coord, radius=self.radius*u.deg)[0] # Just get the first (and only table here)
        result.rename_columns(['_RA.icrs', '_DE.icrs'], ['ra', 'dec']) # Rename the columns to match the SurveyCoord class
        
        # Convert distances from h^-1 Mpc to Mpc based on the cosmology being used.
        result['Dist'] /=self.cosmo.h

        # Apply a transverse distance cut
        angular_dist = self.coord.separation(SkyCoord(result['ra'], result['dec'], unit='deg')).to('rad').value
        transverse_dist = result['Dist']*np.sin(angular_dist)
        result = result[transverse_dist<transverse_distance_cut]
        result = result[result['Nmb']>=richness_cut]
        self.catalog = result

        return self.catalog