"""
Slurp data from GALEX catalog using the MAST API.

"""

import numpy as np

from astropy import units as u,utils as astroutils
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table, join
from ..galaxies.defs import GALEX_bands
from astroquery.mast import Catalogs

from frb.surveys import surveycoord,catalog_utils,images

import os

# Define the data model for GALEX data
photom = {}
photom['GALEX'] = {}
for band in GALEX_bands:
    photom["GALEX"]["GALEX"+'_{:s}'.format(band)] = '{:s}_mag'.format(band.lower())
    photom["GALEX"]["GALEX"+'_{:s}_err'.format(band)] = '{:s}_magerr'.format(band.lower())
    photom["GALEX"]["GALEX_ID"] = 'objID'
photom["GALEX"]['ra'] = 'ra'
photom["GALEX"]['dec'] = 'dec'

# Define the default set of query fields
# See: http://www.galex.caltech.edu/researcher/files/mcat_columns_long.txt
# for additional Fields
_DEFAULT_query_fields = ['distance_arcmin','objID','survey','ra','dec','e_bv']
_DEFAULT_query_fields +=['{:s}_mag'.format(band) for band in GALEX_bands]
_DEFAULT_query_fields +=['{:s}_magerr'.format(band) for band in GALEX_bands]

class GALEX_Survey(surveycoord.SurveyCoord):
    """
    A class to access all the catalogs hosted on the
    MAST database. Inherits from SurveyCoord. This
    is a super class not meant for use by itself and
    instead meant to instantiate specific children
    classes like GALEX_Survey
    """
    def __init__(self,coord,radius,**kwargs):
        surveycoord.SurveyCoord.__init__(self,coord,radius,**kwargs)

        self.Survey = "GALEX"
    
    def get_catalog(self,query_fields=None):
        """
        Query a catalog in the MAST GALEX database for
        photometry.

        Args:
            query_fields: list, optional
                A list of query fields to
                get in addition to the
                default fields.
        
        Returns:
            catalog: astropy.table.Table
                Contains all query results
        """
        #assert self.radius <= 0.5*u.deg, "Cone serches have a maximum radius"
        # url = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/{:s}/{:s}".format(release,table)
        if query_fields is None:
            query_fields = _DEFAULT_query_fields
        else:
            query_fields = _DEFAULT_query_fields+query_fields
        
        data = {}
        data['ra'] = self.coord.ra.value
        data['dec'] = self.coord.dec.value
        data['radius'] = self.radius.to(u.deg).value
        data['columns'] = query_fields
        data['format'] = 'csv'

        ret = Catalogs.query_region(self.coord, radius=self.radius, 
                                       catalog="GALEX")
        print(ret)

        # ret = Table(ret[np.argmin(ret['distance_arcmin'])])
        # .to_pandas().sort_values(by='distance_arcmin')[query_fields]

        pdict = photom['GALEX'].copy()
        
        photom_catalog = catalog_utils.clean_cat(ret,pdict)

        # Remove duplicate entries.
        photom_catalog = catalog_utils.remove_duplicates(photom_catalog, "GALEX_ID")

        self.catalog = catalog_utils.sort_by_separation(photom_catalog, self.coord,
                                                        radec=('ra','dec'), add_sep=True)
        # Meta
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        #Validate
        self.validate_catalog()

        #Return
        return self.catalog.copy()