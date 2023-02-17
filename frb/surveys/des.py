"""DES Survey"""

import numpy as np
from astropy import units, io, utils

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils
from frb.surveys import defs

from IPython import embed

# Dependencies
try:
    from pyvo.dal import sia
except ImportError:
    print("Warning:  You need to install pyvo to retrieve DES images")
    _svc = None
else:
    _svc = sia.SIAService(defs.NOIR_DEF_ACCESS_URL+'des_dr2')

# Define the data model for DES data
photom = {}
photom['DES'] = {}
DES_bands = ['g', 'r', 'i', 'z', 'Y']
for band in DES_bands:
    photom['DES']['DES_{:s}'.format(band)] = 'mag_auto_{:s}'.format(band.lower())
    photom['DES']['DES_{:s}_err'.format(band)] = 'magerr_auto_{:s}'.format(band.lower())
photom['DES']['DES_ID'] = 'coadd_object_id'
photom['DES']['ra'] = 'ra'
photom['DES']['dec'] = 'dec'
photom['DES']['DES_tile'] = 'tilename'
photom['DES']['star_flag_r'] = "class_star_r"
photom['DES']['star_flag_err'] = "spreaderr_model_r"

class DES_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DECaL survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate

    """

    def __init__(self, coord, radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DES'
        self.bands = ['g', 'r', 'i', 'z', 'y']
        self.svc = _svc
        self.qc_profile = "default"
        self.database = "des_dr2.main"

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

    def get_catalog(self, query=None, query_fields=None, 
                    print_query=False, **kwargs):
        """
        Grab a catalog of sources around the input coordinate to the search radius

        Args:
            query: Not used
            query_fields (list, optional): Over-ride list of items to query
            print_query (bool): Print the SQL query generated
            grab_wise (bool): Attempt to grab WISE data too.
                This is not recommended.

        Returns:
            astropy.table.Table:  Catalog of sources returned.  Includes WISE
            photometry for matched sources.
        """
        # Main DES query
        main_cat = super(DES_Survey, self).get_catalog(query_fields=query_fields, print_query=print_query,**kwargs)
        if len(main_cat) == 0:
            main_cat = catalog_utils.clean_cat(main_cat,photom['DES'])
            return main_cat
        main_cat = catalog_utils.clean_cat(main_cat, photom['DES'])
        ## Finish
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
                for key,value in photom['DES'].items():
                    query_fields += [value]
                database = self.database
            else:
                raise IOError("Bad qtype")
        else:
            if qtype == 'main':
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
        url = row['access_url']
        if verbose:
            print ('downloading deepest stacked image...')

        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat

