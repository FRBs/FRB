"""DES Survey"""

import numpy as np
from astropy import units, io, utils
from frb.surveys import dlsurvey
from pyvo.dal import sia

class DES_Survey(dlsurvey.DL_Survey):
    
    def __init__(self,coord,radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self,coord,radius, **kwargs)
        self.survey = 'DES'
        self.bands = ['g', 'r', 'i', 'z', 'Y']
        self.svc = sia.SIAService("https://datalab.noao.edu/sia/des_dr1")
        self.qc_profile = "default"

    def _parse_cat_band(self,band):
        """
        For DES, nothing much is necessary.
        """
        table_cols = ['proctype','prodtype']
        col_vals = ['Stack','image']

        return table_cols, col_vals, band

    def _gen_cat_query(self,query_fields=None):
        """
        Generate SQL Query for catalog search
        """
        if query_fields is None:
            object_id_fields = ['coadd_object_id', 'ra','dec','tilename']
            mag_fields = ['mag_auto_{:s}'.format(band) for band in self.bands]
            mag_err_fields = ['magerr_auto_{:s}'.format(band) for band in self.bands]
            query_fields = object_id_fields+mag_fields+mag_err_fields
        
        database = "des_dr1.main"
        self.query = dlsurvey._default_query_str(query_fields,database,self.coord,self.radius)

    def _select_best_img(self,imgTable,verbose,timeout=120):
        row = imgTable[np.argmax(imgTable['exptime'].data.data.astype('float'))] # pick image with longest exposure time
        url = row['access_url'].decode()
        if verbose:
            print ('downloading deepest stacked image...')
        
        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat

