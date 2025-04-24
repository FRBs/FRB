"""
A module to query for galaxy groups/clusters around a given FRB.
Currently has only the Tully cluster catalog but can be possibly extended for
other sources.
"""

from . import surveycoord
from frb.defs import frb_cosmo
from astropy.coordinates import SkyCoord
from astropy import units as u

try:
    from astroquery.vizier import Vizier
except ImportError:
    print("Warning: You need to install astroquery to use the cluster searches...")

import numpy as np

class VizierCatalogSearch(surveycoord.SurveyCoord):
    """
    A class to query sources within a Vizier catalog.
    """


    def __init__(self, coord, radius = 90*u.deg, survey=None, viziercatalog = None, cosmo=None, **kwargs):
        # Initialize a SurveyCoord object
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        self.survey = survey # Name
        self.viziercatalog = viziercatalog # Name of the Vizier table to draw from.
        self.coord = coord # Location around which to perform the search
        self.radius = radius.to('deg').value # Radius of cone search
        if cosmo is None: # Use the same cosmology as elsewhere in this repository unless specified.
            self.cosmo = frb_cosmo
        else:
            self.cosmo = cosmo
    
    def clean_catalog(self, catalog):
        """
        This will be survey specific.
        """
        pass

    def _transverse_distance_cut(self, catalog, transverse_distance_cut, distance_column='Dist'):
        # Apply a transverse distance cut
        angular_dist = self.coord.separation(SkyCoord(catalog['ra'], catalog['dec'], unit='deg')).to('rad').value
        transverse_dist = catalog[distance_column]*np.sin(angular_dist)
        catalog = catalog[transverse_dist<transverse_distance_cut]
        return catalog

    def _get_catalog(self, query_fields=None, **kwargs):
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
        v = Vizier(catalog = self.viziercatalog, columns=query_fields, row_limit= -1, **kwargs) # No row limit
        result = v.query_region(self.coord, radius=self.radius*u.deg)[0] # Just get the first (and only table here)
        return result

# Tully 2015
class TullyGroupCat(VizierCatalogSearch):
    """
    A class to query sources within the Tully 2015
    group/cluster catalog.
    """


    def __init__(self, coord, radius = 90*u.deg, cosmo=None, **kwargs):
        # Initialize a SurveyCoord object
        super(TullyGroupCat, self).__init__(self, coord, radius,
                                            survey="Tully+2015",
                                            viziercatalog="J/AJ/149/171/table5",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog):

        catalog.rename_columns(['_RA.icrs', '_DE.icrs', 'Nmb'], ['ra', 'dec', 'Ngal']) # Rename the columns to match the SurveyCoord class
        
        # Convert distances from h^-1 Mpc to Mpc based on the cosmology being used.
        catalog['Dist'] /=self.cosmo.h

        return catalog
    

    def get_catalog(self, query_fields=None,
                    transverse_distance_cut = np.inf*u.Mpc, richness_cut = 5):
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
        result = super(TullyGroupCat, self)._get_catalog(query_fields=query_fields)

        result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(TullyGroupCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        result = result[result['Ngal']>=richness_cut]
        self.catalog = result

        return self.catalog
    
# Wen+2024
class WenGroupCat(VizierCatalogSearch):
    """
    A class to query sources within the Wen+2024
    group/cluster catalog.
    """


    def __init__(self, coord, radius = 90*u.deg, cosmo=None, **kwargs):
        # Initialize a SurveyCoord object
        super(WenGroupCat, self).__init__(self, coord, radius,
                                            survey="Wen+2024",
                                            viziercatalog="J/ApJS/272/39/table2",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog):

        catalog.rename_columns(['RAJ2000', 'DEJ2000', 'zcl'], ['ra', 'dec', 'z']) # Rename the columns to match the SurveyCoord class
        
        # Add a distance estimate in Mpc using the given cosmology
        catalog['Dist'] = self.cosmo.lookback_distance(catalog['z']).to('Mpc').value

        return catalog
    
    def get_catalog(self, query_fields=None,
                    transverse_distance_cut = np.inf*u.Mpc, richness_cut = 5):
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
        result = super(WenGroupCat, self)._get_catalog(query_fields=query_fields)

        result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(TullyGroupCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        result = result[result['Ngal']>=richness_cut]
        self.catalog = result
        return self.catalog
    
# Bahk and Hwang 2024 (Updated Planck+2015)
class UPClusterSZCat(VizierCatalogSearch):
    """
    A class to query sources within the Bahk and Hwang 2024
    group/cluster catalog.
    """


    def __init__(self, coord, radius = 90*u.deg, cosmo=None, **kwargs):
        # Initialize a SurveyCoord object
        super(UPClusterSZCat, self).__init__(self, coord, radius,
                                            survey="UPClusterSZ",
                                            viziercatalog="J/ApJS/272/7/table2",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog):

        catalog.rename_columns(['RAJ2000', 'DEJ2000'], ['ra', 'dec']) # Rename the columns to match the SurveyCoord class
        
        # Add a distance estimate in Mpc using the given cosmology
        catalog['Dist'] = self.cosmo.lookback_distance(catalog['z']).to('Mpc').value

        return catalog

    def get_catalog(self, query_fields=None,
                    transverse_distance_cut = np.inf*u.Mpc):
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
        result = super(UPClusterSZCat, self)._get_catalog(query_fields=query_fields)

        result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(UPClusterSZCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        return self.catalog
    
# Xu+2022 (ROSAT X ray cluster)

class ROSATXClusterCat(VizierCatalogSearch):
    """
    A class to query sources within the Xu+2022
    group/cluster catalog.
    """


    def __init__(self, coord, radius = 90*u.deg, cosmo=None, **kwargs):
        # Initialize a SurveyCoord object
        super(ROSATXClusterCat, self).__init__(self, coord, radius,
                                            survey="ROSATXCluster",
                                            viziercatalog="J/A+A/658/A59/table3",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog):

        catalog.rename_columns(['RAJ2000', 'DEJ2000'], ['ra', 'dec']) # Rename the columns to match the SurveyCoord class
        
        # Add a distance estimate in Mpc using the given cosmology
        catalog['Dist'] = self.cosmo.lookback_distance(catalog['z']).to('Mpc').value

        return catalog
    
    def get_catalog(self, query_fields=None,
                    transverse_distance_cut = np.inf*u.Mpc):
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
        result = super(ROSATXClusterCat, self)._get_catalog(query_fields=query_fields)

        result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(ROSATXClusterCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        return self.catalog
    
# Tempel+2018

class TempelClusterCat(VizierCatalogSearch):
    """
    A class to query sources within the Tempel+2018
    group/cluster catalog.
    """


    def __init__(self, coord, radius = 90*u.deg, cosmo=None, **kwargs):
        # Initialize a SurveyCoord object
        super(TempelClusterCat, self).__init__(self, coord, radius,
                                            survey="TempelCluster",
                                            viziercatalog="J/A+A/618/A81/2mrs_gr",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog):

        catalog.rename_columns(['RAJ2000', 'DEJ2000', 'zcmb'], ['ra', 'dec', 'z']) # Rename the columns to match the SurveyCoord class
        
        # Add a distance estimate in Mpc using the given cosmology
        catalog['Dist'] = self.cosmo.lookback_distance(catalog['z']).to('Mpc').value

        return catalog
    
    def get_catalog(self, query_fields=None,
                    transverse_distance_cut = np.inf*u.Mpc):
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
        result = super(TempelClusterCat, self)._get_catalog(query_fields=query_fields)

        result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(TempelClusterCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        return self.catalog