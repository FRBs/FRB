"""
A module to query for galaxy groups/clusters around a given FRB.
Currently has only the Tully cluster catalog but can be possibly extended for
other sources.
"""

from . import surveycoord
from frb.defs import frb_cosmo
from astropy.coordinates import Angle, SkyCoord
from astropy.cosmology import Cosmology
from astropy import units as u
from astropy.table import Table

try:
    from astroquery.vizier import Vizier
except ImportError:
    print("Warning: You need to install astroquery to use the cluster searches...")

import numpy as np

class VizierCatalogSearch(surveycoord.SurveyCoord):
    """
    A class to query sources within a Vizier catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search. Note that it is stored in
            ``self.radius`` as a float in degrees.
        survey (str, optional): Name of the survey.
        viziercatalog (str, optional): Name of the Vizier table to draw from.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 90*u.deg,
                 survey: str | None = None, viziercatalog: str | None = None,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(VizierCatalogSearch, self).__init__(coord, radius, **kwargs)
        self.survey = survey # Name
        self.viziercatalog = viziercatalog # Name of the Vizier table to draw from.
        self.coord = coord # Location around which to perform the search
        self.radius = radius.to('deg').value # Radius of cone search
        if cosmo is None: # Use the same cosmology as elsewhere in this repository unless specified.
            self.cosmo = frb_cosmo
        else:
            self.cosmo = cosmo
    
    def clean_catalog(self, catalog: Table) -> Table:
        """
        This will be survey specific.

        Child classes rename the columns of the Vizier table to the
        standard names ('ra', 'dec', 'z', ...) and add a 'Dist' column.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.

        Returns:
            astropy.table.Table: The cleaned table.
        """
        pass

    def _transverse_distance_cut(self, catalog: Table,
                                 transverse_distance_cut: u.Quantity,
                                 distance_column: str = 'Dist') -> Table:
        """
        Apply a transverse distance cut.

        Args:
            catalog (astropy.table.Table): Table with 'ra' and 'dec' columns
                and a column of distances (in Mpc, if it has no units).
            transverse_distance_cut (astropy.units.Quantity): The maximum
                impact parameter (in length units) of the objects to keep.
            distance_column (str, optional): Name of the column with the
                radial distance.

        Returns:
            astropy.table.Table: The rows of ``catalog`` with transverse
            distances from ``self.coord`` below the cut.

        """
        angular_dist = self.coord.separation(SkyCoord(catalog['ra'], catalog['dec'], unit='deg')).to('rad').value
        radial_dist = catalog[distance_column]
        if getattr(radial_dist, 'unit', None) is None:
            radial_dist = radial_dist * u.Mpc
        else:
            radial_dist = radial_dist.to(u.Mpc)
        transverse_dist = radial_dist * np.sin(angular_dist)
        catalog = catalog[transverse_dist<transverse_distance_cut]
        return catalog

    def _get_catalog(self, query_fields: list[str] | None = None,
                     **kwargs) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            **kwargs: Passed to ``astroquery.vizier.Vizier``.

        Returns:
            astropy.table.Table: The (uncleaned) Vizier table of objects within the search radius.
            If no objects found, returns an empty table with fields ['ra','dec', and 'z'].
        """
        if query_fields is None:
            query_fields = ['**'] # Get all.

        # Query Vizier
        v = Vizier(catalog = self.viziercatalog, columns=query_fields, row_limit= -1, **kwargs) # No row limit
        result = v.query_region(self.coord, radius=self.radius*u.deg)
        if len(result) == 0:
            print("No objects found within the given radius.")
            return Table(names = ('ra','dec','z'))
        else:
            result = result[0]
        return result

# Tully 2015
class TullyGroupCat(VizierCatalogSearch):
    """
    A class to query sources within the Tully 2015
    group/cluster catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 90*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(TullyGroupCat, self).__init__(coord, radius,
                                            survey="Tully+2015",
                                            viziercatalog="J/AJ/149/171/table5",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog: Table) -> Table:
        """
        Rename the columns of the Vizier table to the standard names
        and add a distance estimate.

        The distances are converted from h^-1 Mpc to Mpc with the cosmology first.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.
                Modified in place.

        Returns:
            astropy.table.Table: The cleaned table, with 'ra', 'dec', 'z'
            (where applicable) and 'Dist' (angular diameter distance in Mpc,
            using ``self.cosmo``) columns.

        """
        catalog.rename_columns(['_RA.icrs', '_DE.icrs', 'Nmb'], ['ra', 'dec', 'Ngal']) # Rename the columns to match the SurveyCoord class
        
        # Convert distances from h^-1 Mpc to Mpc based on the cosmology being used.
        catalog['Dist'] /=self.cosmo.h
        rec_velocity = catalog['Dist'].value*self.cosmo.H0.value
        c_kms = 299792.458
        redshift = ((1+rec_velocity/c_kms)/(1-rec_velocity/c_kms))**0.5-1
        catalog['Dist'] = self.cosmo.angular_diameter_distance(redshift).value
        return catalog
    

    def get_catalog(self, query_fields: list[str] | None = None,
                    transverse_distance_cut: u.Quantity = np.inf*u.Mpc,
                    richness_cut: int = 5) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            transverse_distance_cut (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.
            richness_cut (int, optional): The minimum number of members in any group/cluster returned.

        Returns:
            astropy.table.Table: A table of objects within the given limits.
        """
        result = super(TullyGroupCat, self)._get_catalog(query_fields=query_fields)
        if len(result) > 0:
            result = self.clean_catalog(result)

            # Apply a transverse distance cut
            if transverse_distance_cut<np.inf*u.Mpc:
                result = super(TullyGroupCat, self)._transverse_distance_cut(result, transverse_distance_cut)
            result = result[result['Ngal']>=richness_cut]
        self.catalog = result
        
        # Normalize and validate
        self.validate_catalog()

        return self.catalog
    
# Wen+2024
class WenGroupCat(VizierCatalogSearch):
    """
    A class to query sources within the Wen+2024
    group/cluster catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 0.2*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(WenGroupCat, self).__init__(coord, radius,
                                            survey="Wen+2024",
                                            viziercatalog="J/ApJS/272/39/table2",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog: Table) -> Table:
        """
        Rename the columns of the Vizier table to the standard names
        and add a distance estimate.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.
                Modified in place.

        Returns:
            astropy.table.Table: The cleaned table, with 'ra', 'dec', 'z'
            (where applicable) and 'Dist' (angular diameter distance in Mpc,
            using ``self.cosmo``) columns.

        """
        try:
            catalog.rename_columns(['RAJ2000'],['ra'])
        except KeyError:
            assert 'ra' in catalog.keys()
        try:
            catalog.rename_columns(['DEJ2000'],['dec'])
        except KeyError:
            assert 'dec' in catalog.keys()
        try:
            catalog.rename_columns(['zCl'],['z'])
        except KeyError:
            assert 'z' in catalog.keys()

        
        # Add a distance estimate in Mpc using the given cosmology
        catalog['Dist'] = self.cosmo.angular_diameter_distance(catalog['z']).value

        return catalog
    
    def get_catalog(self, query_fields: list[str] | None = None,
                    transverse_distance_cut: u.Quantity = np.inf*u.Mpc,
                    richness_cut: int = 5) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            transverse_distance_cut (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.
            richness_cut (int, optional): The minimum number of members in any group/cluster returned.

        Returns:
            astropy.table.Table: A table of objects within the given limits.
        """
        result = super(WenGroupCat, self)._get_catalog(query_fields=query_fields)
        if len(result) > 0:
            result = self.clean_catalog(result)

            # Apply a transverse distance cut
            if transverse_distance_cut<np.inf*u.Mpc:
                result = super(TullyGroupCat, self)._transverse_distance_cut(result, transverse_distance_cut)
            result = result[result['Ngal']>=richness_cut]
        self.catalog = result
        
        # Normalize and validate
        self.validate_catalog()
        
        return self.catalog
    
# Bahk and Hwang 2024 (Updated Planck+2015)
class UPClusterSZCat(VizierCatalogSearch):
    """
    A class to query sources within the Bahk and Hwang 2024
    group/cluster catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 90*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(UPClusterSZCat, self).__init__(coord, radius,
                                            survey="UPClusterSZ",
                                            viziercatalog="J/ApJS/272/7/table2",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog: Table) -> Table:
        """
        Rename the columns of the Vizier table to the standard names
        and add a distance estimate.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.
                Modified in place.

        Returns:
            astropy.table.Table: The cleaned table, with 'ra', 'dec', 'z'
            (where applicable) and 'Dist' (angular diameter distance in Mpc,
            using ``self.cosmo``) columns.

        """
        if len(catalog) > 0:
            try:
                catalog.rename_columns(['RAJ2000', 'DEJ2000'], ['ra', 'dec']) # Rename the columns to match the SurveyCoord class
            except KeyError:
                print(catalog.keys())

                assert 'ra' in catalog.keys() and 'dec' in catalog.keys()
            
            # Add a distance estimate in Mpc using the given cosmology
            catalog['Dist'] = self.cosmo.angular_diameter_distance(catalog['z']).value

        return catalog

    def get_catalog(self, query_fields: list[str] | None = None,
                    transverse_distance_cut: u.Quantity = np.inf*u.Mpc) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            transverse_distance_cut (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.

        Returns:
            astropy.table.Table: A table of objects within the given limits.
        """
        result = super(UPClusterSZCat, self)._get_catalog(query_fields=query_fields)
        if len(result) > 0:
            result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(UPClusterSZCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        
        # Normalize and validate
        self.validate_catalog()
        
        return self.catalog
    
# Xu+2022 (ROSAT X ray cluster)

class ROSATXClusterCat(VizierCatalogSearch):
    """
    A class to query sources within the Xu+2022
    group/cluster catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 90*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(ROSATXClusterCat, self).__init__(coord, radius,
                                            survey="ROSATXCluster",
                                            viziercatalog="J/A+A/658/A59/table3",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog: Table) -> Table:
        """
        Rename the columns of the Vizier table to the standard names
        and add a distance estimate.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.
                Modified in place.

        Returns:
            astropy.table.Table: The cleaned table, with 'ra', 'dec', 'z'
            (where applicable) and 'Dist' (angular diameter distance in Mpc,
            using ``self.cosmo``) columns.

        """
        if len(catalog) > 0:
            catalog.rename_columns(['RAJ2000', 'DEJ2000'], ['ra', 'dec']) # Rename the columns to match the SurveyCoord class
            
            # Add a distance estimate in Mpc using the given cosmology
            catalog['Dist'] = self.cosmo.angular_diameter_distance(catalog['z']).value

        return catalog
    
    def get_catalog(self, query_fields: list[str] | None = None,
                    transverse_distance_cut: u.Quantity = np.inf*u.Mpc) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            transverse_distance_cut (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.

        Returns:
            astropy.table.Table: A table of objects within the given limits.
        """
        result = super(ROSATXClusterCat, self)._get_catalog(query_fields=query_fields)
        if len(result) > 0:
            result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(ROSATXClusterCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        
        # Normalize and validate
        self.validate_catalog()
        
        return self.catalog
    
# Tempel+2018

class TempelClusterCat(VizierCatalogSearch):
    """
    A class to query sources within the Tempel+2018
    group/cluster catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 90*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(TempelClusterCat, self).__init__(coord, radius,
                                            survey="TempelCluster",
                                            viziercatalog="J/A+A/618/A81/2mrs_gr",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog: Table) -> Table:
        """
        Rename the columns of the Vizier table to the standard names
        and add a distance estimate.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.
                Modified in place.

        Returns:
            astropy.table.Table: The cleaned table, with 'ra', 'dec', 'z'
            (where applicable) and 'Dist' (angular diameter distance in Mpc,
            using ``self.cosmo``) columns.

        """
        if len(catalog) > 0:
            catalog.rename_columns(['RAJ2000', 'DEJ2000', 'zcmb'], ['ra', 'dec', 'z']) # Rename the columns to match the SurveyCoord class
            
            # Add a distance estimate in Mpc using the given cosmology
            catalog['Dist'] = self.cosmo.angular_diameter_distance(catalog['z']).to('Mpc').value

        return catalog
    
    def get_catalog(self, query_fields: list[str] | None = None,
                    transverse_distance_cut: u.Quantity = np.inf*u.Mpc) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            transverse_distance_cut (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.

        Returns:
            astropy.table.Table: A table of objects within the given limits.
        """
        result = super(TempelClusterCat, self)._get_catalog(query_fields=query_fields)
        if len(result) > 0:
            result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(TempelClusterCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        
        # Normalize and validate
        self.validate_catalog()
        
        return self.catalog


# Klein+ 2023 RASS-MCMF (Rosat All Sky Survey Multi Component Matched Filter)

class RASSClusterCat(VizierCatalogSearch):
    """
    A class to query sources within the Klein+ 2023
    group/cluster catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 90*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(RASSClusterCat, self).__init__(coord, radius,
                                            survey="RASSCluster",
                                            viziercatalog="J/MNRAS/526/3757/catalog",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog: Table) -> Table:
        """
        Rename the columns of the Vizier table to the standard names
        and add a distance estimate.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.
                Modified in place.

        Returns:
            astropy.table.Table: The cleaned table, with 'ra', 'dec', 'z'
            (where applicable) and 'Dist' (angular diameter distance in Mpc,
            using ``self.cosmo``) columns.

        """
        if len(catalog) > 0:
            catalog.rename_columns(['RAJ2000', 'DEJ2000','zsp1'], ['ra', 'dec','z']) # Rename the columns to match the SurveyCoord class
            
            # Add a distance estimate in Mpc using the given cosmology
            catalog['Dist'] = self.cosmo.angular_diameter_distance(catalog['z']).value

        return catalog
    
    def get_catalog(self, query_fields: list[str] | None = None,
                    transverse_distance_cut: u.Quantity = np.inf*u.Mpc) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            transverse_distance_cut (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.

        Returns:
            astropy.table.Table: A table of objects within the given limits.
        """
        result = super(RASSClusterCat, self)._get_catalog(query_fields=query_fields)
        if len(result) > 0:
            result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(RASSClusterCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        
        # Normalize and validate
        self.validate_catalog()
        
        return self.catalog

    
# Rykoff+ 2014 Clusters identified in SDSS8 using red mapper algorithm

class RedMapperClusterCat(VizierCatalogSearch):
    """
    A class to query sources within the Rykoff+ 2014
    group/cluster catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 90*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(RedMapperClusterCat, self).__init__(coord, radius,
                                            survey="RedMapperCluster",
                                            viziercatalog="J/ApJ/785/104/table1",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog: Table) -> Table:
        """
        Rename the columns of the Vizier table to the standard names
        and add a distance estimate.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.
                Modified in place.

        Returns:
            astropy.table.Table: The cleaned table, with 'ra', 'dec', 'z'
            (where applicable) and 'Dist' (angular diameter distance in Mpc,
            using ``self.cosmo``) columns.

        """
        if len(catalog) > 0:
            catalog.rename_columns(['RAJ2000', 'DEJ2000','zspec'], ['ra', 'dec','z']) # Rename the columns to match the SurveyCoord class
            
            # Add a distance estimate in Mpc using the given cosmology
            redshift = catalog['z']
            redshift[redshift<0] = catalog['zlambda'][redshift<0]
            catalog['Dist'] = self.cosmo.angular_diameter_distance(redshift).value

        return catalog
    
    def get_catalog(self, query_fields: list[str] | None = None,
                    transverse_distance_cut: u.Quantity = np.inf*u.Mpc) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            transverse_distance_cut (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.

        Returns:
            astropy.table.Table: A table of objects within the given limits.
        """
        result = super(RedMapperClusterCat, self)._get_catalog(query_fields=query_fields)
        if len(result) > 0:
            result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(RedMapperClusterCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        
        # Normalize and validate
        self.validate_catalog()
        
        return self.catalog

    
# Klein+ 2024 Clusters identified in ACT data release 5

class ACTDR5ClusterCat(VizierCatalogSearch):
    """
    A class to query sources within the Klein+ 2024
    group/cluster catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 90*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(ACTDR5ClusterCat, self).__init__(coord, radius,
                                            survey="ACTDR5Cluster",
                                            viziercatalog="J/A+A/690/A322/catalog",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog: Table) -> Table:
        """
        Rename the columns of the Vizier table to the standard names
        and add a distance estimate.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.
                Modified in place.

        Returns:
            astropy.table.Table: The cleaned table, with 'ra', 'dec', 'z'
            (where applicable) and 'Dist' (angular diameter distance in Mpc,
            using ``self.cosmo``) columns.

        """
        if len(catalog) > 0:
            catalog.rename_columns(['RAJ2000', 'DEJ2000', 'z1C'], ['ra', 'dec','z']) # Rename the columns to match the SurveyCoord class
            
            # Add a distance estimate in Mpc using the given cosmology
            redshift = catalog['z']
            catalog['Dist'] = self.cosmo.angular_diameter_distance(redshift).value

        return catalog
    
    def get_catalog(self, query_fields: list[str] | None = None,
                    transverse_distance_cut: u.Quantity = np.inf*u.Mpc) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            transverse_distance_cut (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.

        Returns:
            astropy.table.Table: A table of objects within the given limits.
        """
        result = super(ACTDR5ClusterCat, self)._get_catalog(query_fields=query_fields)
        if len(result) > 0:
            result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(ACTDR5ClusterCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        
        # Normalize and validate
        self.validate_catalog()
        
        return self.catalog

    
# Kluge+ 2024 Clusters identified in eRosita x-ray all sky survey

class ERASSClusterCat(VizierCatalogSearch):
    """
    A class to query sources within the Kluge+ 2024
    group/cluster catalog.

    Args:
        coord (astropy.coordinates.SkyCoord): Location around which to
            perform the search.
        radius (astropy.coordinates.Angle or astropy.units.Quantity, optional):
            Radius of the cone search.
        cosmo (astropy.cosmology.Cosmology, optional): Cosmology used for
            distances. Defaults to ``frb.defs.frb_cosmo``.
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """


    def __init__(self, coord: SkyCoord, radius: Angle | u.Quantity = 90*u.deg,
                 cosmo: Cosmology | None = None, **kwargs):
        # Initialize a SurveyCoord object
        super(ERASSClusterCat, self).__init__(coord, radius,
                                            survey="ERASSCluster",
                                            viziercatalog="J/A+A/688/A210/tablee1",
                                            cosmo=cosmo,  **kwargs)
        
    def clean_catalog(self, catalog: Table) -> Table:
        """
        Rename the columns of the Vizier table to the standard names
        and add a distance estimate.

        Args:
            catalog (astropy.table.Table): Table returned by Vizier.
                Modified in place.

        Returns:
            astropy.table.Table: The cleaned table, with 'ra', 'dec', 'z'
            (where applicable) and 'Dist' (angular diameter distance in Mpc,
            using ``self.cosmo``) columns.

        """
        if len(catalog) > 0:
            catalog.rename_columns(['RAJ2000', 'DEJ2000', 'Bestz'], ['ra', 'dec','z']) # Rename the columns to match the SurveyCoord class
            
            # Add a distance estimate in Mpc using the given cosmology
            redshift = catalog['z']
            catalog['Dist'] = self.cosmo.angular_diameter_distance(redshift).value

        return catalog
    
    def get_catalog(self, query_fields: list[str] | None = None,
                    transverse_distance_cut: u.Quantity = np.inf*u.Mpc) -> Table:
        """
        Get the catalog of objects

        Args:
            query_fields (list of str, optional): The fields to include in the catalog. If None, all fields are used.
            transverse_distance_cut (astropy.units.Quantity, optional): The maximum impact parameter of the objects to include in the catalog.

        Returns:
            astropy.table.Table: A table of objects within the given limits.
        """
        result = super(ERASSClusterCat, self)._get_catalog(query_fields=query_fields)
        if len(result) > 0:
            result = self.clean_catalog(result)

        # Apply a transverse distance cut
        if transverse_distance_cut<np.inf*u.Mpc:
            result = super(ERASSClusterCat, self)._transverse_distance_cut(result, transverse_distance_cut)
        self.catalog = result
        
        # Normalize and validate
        self.validate_catalog()
        
        return self.catalog
