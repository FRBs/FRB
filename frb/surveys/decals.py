"""DECaLS"""

import numpy as np
from astropy import units, io, utils
from frb.surveys.dlsurvey import DL_Survey
from pyvo.dal import sia

class DECaL_Survey(DL_Survey):

    def __init__(self,coord,radius, **kwargs):
        DL_Survey.__init__(self,coord,radius, **kwargs)
        self.survey = 'DECaL'
        self.bands = ['g', 'r', 'z']
        self.svc = sia.SIAService("https://datalab.noao.edu/sia/ls_dr7")
        self.qc_profile = "default"
    
    def _parse_cat_band(self,band):
        if band is  'g':
            bandstr = "g DECam SDSS c0001 4720.0 1520.0"
        elif band is 'r':
            bandstr = "r DECam SDSS c0002 6415.0 1480.0	"
        elif band is 'z':
            bandstr = "z DECam SDSS c0004 9260.0 1520.0	"
        table_cols = ['prodtype']
        col_vals = ['image']
        return table_cols, col_vals, bandstr
    
    def _gen_cat_query(self,query_fields=None):
        """
        Generate SQL query for catalog search
        """
        self.query = "SELECT "
        if query_fields is None:
            object_id_fields = ['decals_id','brick_primary','brickid','ra','dec']
            mag_fields = ['mag_g','mag_r','mag_z','mag_w1','mag_w2','mag_w3','mag_w4']
            snr_fields = ['snr_g','snr_r','snr_z','snr_w1','snr_w2','snr_w3','snr_w4']
            query_fields = object_id_fields+mag_fields+snr_fields
        
        for field in query_fields:
            self.query += "{:s},".format(field)
        
        self.query = self.query[:-1] #remove the last comma
        self.query += "\nFROM ls_dr7.tractor"
        self.query += "\nWHERE q3c_radial_query(ra,dec,{:f},{:f},{:f})".format(
            self.coord.ra.value,self.coord.dec.value,self.radius.to(units.deg).value)

    def _select_best_img(self,imgTable,verbose):
        """
        Just one image here so no problem
        """
        row = imgTable[0]
        if verbose:
            print ('downloading deepest stacked image...')
        
        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat