"""DES Survey"""

import pdb

import numpy as np
from astropy.table import Table
from astropy import units, io, utils
from pyvo.dal import sia

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils

photom = {}
#DES
photom['DES'] = {}
DES_bands = ['g', 'r', 'i', 'z', 'Y']
for band in DES_bands:
    photom['DES']['DES_{:s}'.format(band)] = 'mag_auto_{:s}'.format(band.lower())
    photom['DES']['DES_{:s}_err'.format(band)] = 'magerr_auto_{:s}'.format(band.lower())
photom['DES']['DES_ID'] = 'coadd_object_id'
photom['DES']['ra'] = 'ra'
photom['DES']['dec'] = 'dec'
photom['DES']['DES_tile'] = 'tilename'

# Dependencies
try:
    from dl import queryClient as qc, authClient as ac, helpers
except ImportError:
    print("Warning:  You need to install datalab to query DES data")
else:
    token = ac.login('anonymous')
    try:
        from pyvo.dal import sia
    except ImportError:
        print("Warning:  You need to install pyvo to retrieve DES images")
    else:
        _DEF_ACCESS_URL = "https://datalab.noao.edu/sia/des_dr1"
        _svc = sia.SIAService(_DEF_ACCESS_URL)

# DES-WISE
photom['DES-WISE'] = {}
DES_WISE_bands = ['W1', 'W2', 'W3', 'W4']
for band in DES_WISE_bands:
    photom['DES-WISE']['{:s}'.format(band)] = '{:s}mpro'.format(band.lower())
    photom['DES-WISE']['{:s}_err'.format(band)] = '{:s}sigmpro'.format(band.lower())
photom['DES-WISE']['DES_ID'] = 'coadd_object_id'
photom['DES-WISE']['DES_ra'] = 'des_ra'
photom['DES-WISE']['DES_dec'] = 'des_dec'
photom['DES-WISE']['WISE_ra'] = 'ra'
photom['DES-WISE']['WISE_dec'] = 'dec'


class DES_Survey(dlsurvey.DL_Survey):

    def __init__(self,coord,radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self,coord,radius, **kwargs)
        self.survey = 'DES'
        self.bands = ['g', 'r', 'i', 'z', 'Y']
        self.svc = sia.SIAService("https://datalab.noao.edu/sia/des_dr1")
        self.qc_profile = "default"
        self.database = "des_dr1.main"

    def _parse_cat_band(self,band):
        """
        For DES, nothing much is necessary.
        """
        table_cols = ['proctype','prodtype']
        col_vals = ['Stack','image']

        return table_cols, col_vals, band

    def get_catalog(self, query=None, query_fields=None, print_query=False):
        # Main query
        main_cat = super().get_catalog(query_fields=query_fields, print_query=print_query)
        if len(main_cat) == 0:
            return main_cat
        main_cat = catalog_utils.clean_cat(main_cat, photom['DES'])


        # WISE
        wise_query = self._gen_cat_query(qtype='wise')
        wise_cat = super().get_catalog(query=wise_query, print_query=print_query)
        wise_cat = catalog_utils.clean_cat(wise_cat, photom['DES-WISE'], fill_mask=-999.)

        # Match em up
        if len(wise_cat) > 0:
            idx = catalog_utils.match_ids(wise_cat['DES_ID'], main_cat['DES_ID'])
            # Fill me
            for band in DES_WISE_bands:
                main_cat['{:s}'.format(band)] = -999.
                main_cat['{:s}'.format(band)][idx] = wise_cat['{:s}'.format(band)]
                main_cat['{:s}_err'.format(band)] = -999.
                main_cat['{:s}_err'.format(band)][idx] = wise_cat['{:s}_err'.format(band)]

        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog

    def _gen_cat_query(self,query_fields=None, qtype='main'):
        """
        Generate SQL Query for catalog search
        """
        if query_fields is None:
            query_fields = []
            # Main query
            if qtype == 'main':
                for key,value in photom['DES'].items():
                    query_fields += [value]
                database = self.database
            elif qtype == 'wise':
                for key,value in photom['DES-WISE'].items():
                    query_fields += [value]
                database = "des_dr1.des_allwise"
            else:
                raise IOError("Bad qtype")

        self.query = dlsurvey._default_query_str(query_fields, database,self.coord,self.radius)
        # Return
        return self.query

    def _select_best_img(self,imgTable,verbose,timeout=120):
        row = imgTable[np.argmax(imgTable['exptime'].data.data.astype('float'))] # pick image with longest exposure time
        url = row['access_url'].decode()
        if verbose:
            print ('downloading deepest stacked image...')

        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat

