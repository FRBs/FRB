"""
Slurp data from the Vizier archive. This is
a super-class for the PAN-STARRS and WISE
survey classes.
"""

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table

import warnings

try:
    from astroquery.vizier import Vizier
except ImportError:
    warnings.warn("Warning: You need to install astroquery to use the survey tools...")

from frb.surveys import surveycoord,catalog_utils,images

class Vizier_Catalog(surveycoord.SurveyCoord):
    """
    A class to access all the catalogs hosted on the
    Vizier database. Inherits from SurveyCoord. This
    is a super class not meant for use by itself and
    instead meant to instantiate specific children
    classes like PAN-STARRS_Survey
    """
    def __init__(self,coord,radius,**kwargs):
        surveycoord.SurveyCoord.__init__(self,coord,radius,**kwargs)

        #Define photometric band names
        self.Survey = ''
    
    def get_catalog(self,catalog=None,query_fields=['*'],timeout=120,print_query=False):
        """
        Query a catalog in the VizieR database for
        photometry.

        Args:
            catalog: str
                Catalog name to search.
            query_fields: list, optional
                A list of query fields to
                get. Gets all fields by default.
            timeout: float, optional
                Query timeout in sec. Exits with an error
                if the query time exceeds this value.
                Default value of 120 s.
            print_query: bool, optional
                If true, prints the SQL query used
                on screen.
        
        Returns:
            catalog: astropy.table.Table
                Contains all query results
        """
        
        