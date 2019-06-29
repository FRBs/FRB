"""
Slurp data from the Vizier archive. Mainly
for Pan-STARRS
"""

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table
from ..galaxies.defs import PanSTARRS_bands

import warnings

try:
    from astroquery.vizier import Vizier
except ImportError:
    warnings.warn("Warning: You need to install astroquery to use the survey tools...")

from frb.surveys import surveycoord,catalog_utils,images

#Just to check if the catalog argument is valid
#More can be added if necessary.
_ALLOWED_catalogs = ["Pan-STARRS","GALEX","2MASX"]
_DEFAULT_catalog = "Pan-STARRS"
_DEFAULT_bands = PanSTARRS_bands

# Define the data model for DES data
def _gen_photom(catalog,bands):
    photom = {}
    photom['PA'] = {}
    if catalog == "Pan-STARRS":
        for band in bands:
            photom[catalog][catalog+'_{:s}'.format(band)] = '{:s}mag'.format(band.lower())
            photom[catalog][catalog+'_{:s}'.format(band)] = 'e_{:s}mag'.format(band.lower())
            photom[catalog][catalog+'_ID'] = 'objID'
    else:
        for band in bands:
            photom[catalog][catalog+'_{:s}'.format(band)] = '{:s}'.format(band.lower())
            photom[catalog][catalog+'_{:s}'.format(band)] = 'e_{:s}'.format(band.lower())
            photom[catalog][catalog+'_ID'] = 'objid'
    photom[catalog]['ra'] = 'ra'
    photom[catalog]['dec'] = 'dec'
    photom[catalog][catalog+'_field'] = 'field'
    return photom

class Vizier_Database(surveycoord.SurveyCoord):
    """
    A class to access all the catalogs hosted on the
    Vizier database. Inherits from SurveyCoord. This
    is a super class not meant for use by itself and
    instead meant to instantiate specific children
    classes like PAN-STARRS_Survey
    """
    def __init__(self,coord,radius,catalog=_DEFAULT_catalog,**kwargs):
        surveycoord.SurveyCoord.__init__(self,coord,radius,**kwargs)

        self.Survey = catalog
    
    def get_catalog(self,catalog=_DEFAULT_catalog,query_fields=None,timeout=120,print_query=False):
        """
        Query a catalog in the VizieR database for
        photometry.

        Args:
            catalog: str, optional
                Catalog name to search. PanSTARRS by default
            query_fields: list, optional
                A list of query fields to
                get.
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
        assert catalog in _ALLOWED_catalogs, "Searches from only these catalogs supported: {:s}".format(str(_ALLOWED_catalogs))

        if catalog is _DEFAULT_catalog:
            query_fields = ['objID','RAJ2000','DEJ2000','f_objID','Qual']
            query_fields +=['{:s}mag'.format(band) for band in _DEFAULT_bands]
            query_fields +=['e_{:s}mag'.format(band) for band in _DEFAULT_bands]
        vclient = Vizier(columns=query_fields)
        tablelist = vclient.query_region(self.coord,radius=self.radius,
                                            timeout=timeout,catalog=catalog)
        if len(tablelist)==0:
            self.catalog = Table()
            self.catalog.meta['radius'] = self.radius
            self.catalog.meta['survey'] = self.survey
            # Validate
            self.validate_catalog()
            return self.catalog.copy()
        photom_catalog = tablelist[0]
        pdict = _gen_photom(catalog,bands)
        photom_catalog = catalog_utils.clean_cat(photom_catalog,)
        self.catalog = catalog_utils.sort_by_separation(photom_catalog, self.coord,
                                                        radec=('ra','dec'), add_sep=True)
        # Meta
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        #Validate
        self.validate_catalog()

        #Return
        return self.catalog.copy()
        

        
        