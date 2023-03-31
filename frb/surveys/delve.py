"""DELVE survey"""

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
    _DEF_ACCESS_URL = "https://datalab.noao.edu/sia/delve_dr2"
    _svc = sia.SIAService(_DEF_ACCESS_URL)

# Define the data model for DELVE data
# See https://datalab.noirlab.edu/query.php?name=delve_dr2.objects for
# table schema
photom = {}
photom['DELVE'] = {}
photom['DELVE']['DELVE_ID'] = 'quick_object_id'
photom['DELVE']['ra'] = 'ra'
photom['DELVE']['dec'] = 'dec'
photom['DELVE']['ebv'] = 'ebv' # Schegel, Finkbeiner, Davis (1998)
DELVE_bands = ['g', 'r', 'i', 'z']
for band in DELVE_bands:
    photom['DELVE'][f'DELVE_{band}'] = f'mag_auto_{band}' #mag
    photom['DELVE'][f'DELVE_{band}_err'] = f'magerr_auto_{band}' #magerr
    photom['DELVE'][f'class_star_{band}'] = f'extended_class_{band}' #morphology class



class DELVE_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DELVE survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate

    """

    def __init__(self, coord, radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DELVE'
        self.bands = DELVE_bands
        self.svc = sia.SIAService("https://datalab.noao.edu/sia/delve_dr2")
        self.qc_profile = "default"
        self.database = "delve_dr2.objects"

    def _parse_cat_band(self,band):
        """
        Internal method to generate the bands for grabbing
        a cutout image

        For DES, nothing much is necessary.

        Args:
            band (str): Band desired

        Returns:
            list, list, str:  Table columns, Column values, band string for cutout

        """
        table_cols = ['proctype','prodtype']
        col_vals = ['Stack','image']

        return table_cols, col_vals, band

    def get_catalog(self, query=None, query_fields=None, print_query=False,**kwargs):
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
        main_cat = super(DELVE_Survey, self).get_catalog(query_fields=query_fields, print_query=print_query,**kwargs)
        if len(main_cat) == 0:
            main_cat = catalog_utils.clean_cat(main_cat,photom['DELVE'])
            return main_cat
        main_cat = catalog_utils.clean_cat(main_cat, photom['DELVE'])
        #import pdb; pdb.set_trace()
        for col in main_cat.colnames:
            if main_cat[col].dtype==float:
                mask = np.isnan(main_cat[col])+(main_cat[col]==99.99)
                main_cat[col] = np.where(~mask, main_cat[col], -999.0)
        
        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog

    def _gen_cat_query(self,query_fields=None, qtype='main'):
        """
        Generate SQL Query for catalog search

        self.query is modified in place

        Args:
            query_fields (list):  Override the default list for the SQL query

        """
        if query_fields is None:
            query_fields = []
            # Main query
            if qtype == 'main':
                for key,value in photom['DELVE'].items():
                    query_fields += [value]
                database = self.database
            else:
                raise IOError("Bad qtype")

        self.query = dlsurvey._default_query_str(query_fields, database,self.coord,self.radius)
        # Return
        return self.query

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

