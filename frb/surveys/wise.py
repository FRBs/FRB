"""WISE Survey"""

from IPython import embed

import numpy as np
from astropy import units, io, utils

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils

# Dependencies
try:
    from pyvo.dal import sia
except ImportError:
    print("Warning:  You need to install pyvo to retrieve DES images")
    _svc = None
else:
    _DEF_ACCESS_URL = "https://datalab.noao.edu/sia/des_dr1"
    _svc = sia.SIAService(_DEF_ACCESS_URL)

# Define the data model for DES data
photom = {}

# DES-WISE
photom['WISE'] = {}
DES_WISE_bands = ['W1', 'W2', 'W3', 'W4']
for band in DES_WISE_bands:
    photom['WISE']['{:s}'.format(band)] = '{:s}mpro'.format(band.lower())
    photom['WISE']['{:s}_err'.format(band)] = '{:s}sigmpro'.format(band.lower())
photom['WISE']['ra'] = 'ra'
photom['WISE']['dec'] = 'dec'


class WISE_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the WISE survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate

    """

    def __init__(self, coord, radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'WISE'
        self.bands = ['W1', 'W2', 'W3', 'W4']
        self.svc = sia.SIAService("https://datalab.noao.edu/sia/allwise")
        self.qc_profile = "default"
        self.database = "allwise.source"

    def get_catalog(self, query=None, query_fields=None, print_query=False):
        """
        Grab a catalog of sources around the input coordinate to the search radius

        Args:
            query: Not used
            query_fields (list, optional): Over-ride list of items to query
            print_query (bool): Print the SQL query generated

        Returns:
            astropy.table.Table:  Catalog of sources returned.  Includes WISE
            photometry for matched sources.
        """
        # Main DES query
        main_cat = super(WISE_Survey, self).get_catalog(query_fields=query_fields, print_query=print_query)
        if len(main_cat) == 0:
            return main_cat
        main_cat = catalog_utils.clean_cat(main_cat, photom['WISE'], fill_mask=-999.)

        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog

    def _gen_cat_query(self,query_fields=None):
        """
        Generate SQL query for catalog search

        self.query is modified in place

        Args:
            query_fields (list):  Override the default list for the SQL query

        """
        if query_fields is None:
            object_id_fields = ['source_id','ra','dec','tmass_key']
            mag_fields = ['w1mpro', 'w2mpro', 'w3mpro', 'w4mpro',
                          'w1sigmpro', 'w2sigmpro', 'w3sigmpro', 'w4sigmpro', 'ph_qual',
                          'moon_lev']
            query_fields = object_id_fields+mag_fields

        self.query = dlsurvey._default_query_str(query_fields, self.database, self.coord, self.radius)

    def _select_best_img(self,imgTable,verbose,timeout=120):
        """
        Select the best band for a cutout

        Args:
            imgTable: Table of images
            verbose (bool):  Print status
            timeout (int or float):  How long to wait before timing out, in seconds

        Returns:
            HDU: header data unit for the downloaded image

        """
        row = imgTable[np.argmax(imgTable['exptime'].data.data.astype('float'))] # pick image with longest exposure time
        url = row['access_url'].decode()
        if verbose:
            print ('downloading deepest stacked image...')

        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat

