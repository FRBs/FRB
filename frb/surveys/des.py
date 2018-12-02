"""DES Survey"""

import numpy as np
from astropy import units, io, utils
from frb.surveys.dlsurvey import DL_Survey
from pyvo.dal import sia

class DES_Survey(DL_Survey):
    
    def __init__(self,coord,radius, **kwargs):
        DL_Survey.__init__(self,coord,radius, **kwargs)
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
        self.query = "SELECT "
        if query_fields is None:
            object_id_fields = ['coadd_object_id', 'ra','dec','tilename']
            mag_fields = ['mag_auto_{:s}'.format(band) for band in self.bands]
            mag_err_fields = ['magerr_auto_{:s}'.format(band) for band in self.bands]
            query_fields = object_id_fields+mag_fields+mag_err_fields
        
        for field in query_fields:
            self.query += "{:s},".format(field)
        
        self.query = self.query[:-1] #remove the last comma
        self.query += "\nFROM des_dr1.main"
        self.query += "\nWHERE q3c_radial_query(ra,dec,{:f},{:f},{:f})".format(
            self.coord.ra.value,self.coord.dec.value,self.radius.to(units.deg).value)
    
    def _select_best_img(self,imgTable,verbose,timeout=120):
        row = imgTable[np.argmax(imgTable['exptime'].data.data.astype('float'))] # pick image with longest exposure time
        url = row['access_url'].decode()
        if verbose:
            print ('downloading deepest stacked image...')
        
        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat

