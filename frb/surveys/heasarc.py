""" Surveys to be accessed through the HEASARC interface (via astroquery"""

from astropy.table import Table
from astropy import units

try:
    from astroquery.heasarc import Heasarc
except ImportError:
    print("Warning:  You need astroquery installed to use the surveys from HEASARC and SkyView")

from frb.surveys import surveycoord
from frb.surveys import catalog_utils
from frb.surveys.skyview import SkyView_Survey


class HEASARC_Survey(surveycoord.SurveyCoord):
    """
        Class to handle queries on the HEASARC survey.
        Uses `astroquery` for searching the Heasarc SQL database.

    
    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate
        mission (str): Mission served by HEASAR for the data searches
    
    """
    def __init__(self, coord, radius, mission, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        #
        self.survey = None
        self.mission = mission
        # Instantiate astroquery object
        self.heasarc = Heasarc()

    def get_catalog(self):
        """
        Grab a catalog of sources around the input coordinate to the search radius

        
        Returns:
            astropy.table.Table:  Catalog of sources returned
        """
        try:
            catalog = self.heasarc.query_region(self.coord,
                                                catalog=self.mission,
                                                radius=self.radius)
        except (ValueError, TypeError):  # No table found
            self.catalog = catalog_utils.ensure_empty_schema(Table(), ['ra', 'dec'])
        else:
            # Clean
            if len(catalog)!=0:
                if "RA" in catalog.colnames:
                    catalog.rename_column("RA", "ra")
                if "DEC" in catalog.colnames:
                    catalog.rename_column("DEC", "dec")
                for key in ['ra', 'dec']:
                    catalog[key].unit = units.deg
                # Sort
                self.catalog = catalog_utils.sort_by_separation(catalog,
                                                                self.coord,
                                                                radec=('ra', 'dec'))
            else:
                self.catalog = catalog
        if len(self.catalog) == 0:
            self.catalog = catalog_utils.ensure_empty_schema(self.catalog, ['ra', 'dec'])
        # Add meta, etc.
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey
        # Validate
        self.validate_catalog()
        # Return
        return self.catalog


class NVSS_Survey(HEASARC_Survey, SkyView_Survey):
    """ Uses SkyView an HEASARC to get both images and catalogs for the VLA NVSS survey at 1.4 GHz.
    """

    def __init__(self, coord, radius, **kwargs):
        HEASARC_Survey.__init__(self, coord, radius, 'nvss', **kwargs)
        SkyView_Survey.__init__(self, coord, radius, 'nvss', **kwargs)
        self.survey = 'NVSS'


class FIRST_Survey(HEASARC_Survey, SkyView_Survey):
    """ Uses SkyView an HEASARC to get both images and catalogs for the VLA FIRST survey at 1.4 GHz.
    """
    def __init__(self, coord, radius, **kwargs):
        HEASARC_Survey.__init__(self, coord, radius, 'first', **kwargs)
        SkyView_Survey.__init__(self, coord, radius, 'first', **kwargs)
        self.survey = 'FIRST'


class WENSS_Survey(HEASARC_Survey, SkyView_Survey):
    """ Uses SkyView an HEASARC to get both images and catalogs for the WSRT northern sky survey at 325 MHz.
    """
    def __init__(self, coord, radius, **kwargs):
        HEASARC_Survey.__init__(self, coord, radius, 'wenss', **kwargs)
        SkyView_Survey.__init__(self, coord, radius, 'wenss', **kwargs)
        self.survey = 'WENSS'
