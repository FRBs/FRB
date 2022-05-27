"""VISTA catalog"""

import numpy as np
from astropy import units, io, utils

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils
from frb.galaxies.defs import VISTA_bands

# Dependencies
try:
    from pyvo.dal import sia
except ImportError:
    print("Warning:  You need to install pyvo to retrieve VISTA images")
    _svc = None
else:
    _DEF_ACCESS_URL = "https://datalab.noao.edu/sia/vhs_dr5"
    _svc = sia.SIAService(_DEF_ACCESS_URL)

# Define the data model for DES data
photom = {}
photom['VISTA'] = {}
photom['VISTA']['VISTA_ID'] = 'sourceid'
photom['VISTA']['ra'] = 'ra2000'
photom['VISTA']['dec'] = 'dec2000'
for band in VISTA_bands:
    photom['VISTA']['VISTA_{:s}'.format(band)] = '{:s}petromag'.format(band.lower())
    photom['VISTA']['VISTA_{:s}_err'.format(band)] = '{:s}petromagerr'.format(band.lower())



class VISTA_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DECaL survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate

    """

    def __init__(self, coord, radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'VISTA'
        self.bands = VISTA_bands
        self.svc = _svc
        self.qc_profile = "default"
        self.database = "vhs_dr5.vhs_cat_v3"

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
                for key,value in photom['VISTA'].items():
                    query_fields += [value]
                database = self.database
            else:
                raise IOError("Bad qtype")

        self.query = dlsurvey._default_query_str(query_fields, database,self.coord,self.radius)

        # Because they HAD to include the epoch in the colname.
        self.query = self.query.replace('ra,dec,','ra2000,dec2000,')
        # Return
        return self.query

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
        if query==None:
            self.query = self._gen_cat_query(query_fields=query_fields)
        else:
            self.query = query
        main_cat = super(VISTA_Survey, self).get_catalog(query=self.query, print_query=print_query,
                                                         photomdict=photom['VISTA'],**kwargs)
        if len(main_cat) == 0:
            main_cat = catalog_utils.clean_cat(main_cat,photom['VISTA'])
            return main_cat
        main_cat = catalog_utils.clean_cat(main_cat, photom['VISTA'])
        for col in main_cat.colnames:
            if main_cat[col].dtype==float:
                mask = np.isnan(main_cat[col])+(main_cat[col]==99.99)
                main_cat[col] = np.where(~mask, main_cat[col], -999.0)

        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog


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

