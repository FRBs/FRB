"""
DataLab survey class. Gets data from any survey
available through the NOAO datalab-client.
"""
import pdb

import numpy as np
import warnings
from astropy.table import Table
from astropy import units, io, utils
try:
    from dl import queryClient as qc, authClient as ac, helpers
except ImportError:
    print("Warning:  You need to install data")

from frb.surveys import surveycoord

class DL_Survey(surveycoord.SurveyCoord):
    def __init__(self, coord, radius, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        
        #Define photmetric band names.
        self.token = ac.login('anonymous')
        self.bands = None
        #Instantiate sia service 
        self.svc = None
        #Generate query
        self.qc_profile = None
    
    def _parse_cat_band(self,band):
        pass
    
    def _gen_cat_query(self,query_fields=None):
        pass
    
    def _select_best_img(self,imgTable,verbose,timeout=120):
        pass

    
    def get_catalog(self,query_fields=None):
        """
        Get catalog sources around the given coordinates
        within a radius.
        Returns
        -------
        cat: astropy Table
            Table of objects obtained from the 
            SQL query.
        """
        qc.set_profile(self.qc_profile)

        self._gen_cat_query(query_fields)
        result = qc.query(self.token,sql=self.query)
        cat = helpers.convert(result)
        # TODO:: Suppress the print output from convert
        # TODO:: Dig into why the heck it doesn't want to natively
        #        output to a table when it was clearly intended to with 'outfmt=table'
        # Finish
        self.catalog = Table.from_pandas(cat)
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey
        # Validate
        self.validate_catalog()
        # Return
        return self.catalog.copy()
    
    def get_image(self,imsize,band,timeout=120,verbose=False):
        """
        Get images from the catalog if available
        for a given fov and band.
        """
        ra = self.coord.ra.value
        dec = self.coord.dec.value
        fov = imsize.to(units.deg).value
        
        if band.lower() not in self.bands:
            raise TypeError("Allowed filters (case-insensitive) for {:s} photometric bands are {}".format(self.survey,self.bands))

        table_cols,col_vals,bandstr = self._parse_cat_band(band)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            imgTable = self.svc.search((ra,dec), (fov/np.cos(dec*np.pi/180), fov), verbosity=2).to_table()
        if verbose:
            print("The full image list contains", len(imgTable), "entries")
        
        #Select band
        selection = imgTable['obs_bandpass'].astype(str)==bandstr

        #Select images in that band
        for column, value in zip(table_cols,col_vals):
            selection = selection & ((imgTable[column].astype(str)==value))
        imgTable = imgTable[selection]

        if(len(imgTable)>0):
            imagedat = self._select_best_img(imgTable,verbose,timeout)
            image = imagedat[0]
        else:
            print('No image available')
            image = None
        return image
    
    def get_cutout(self, imsize,band=None):
        """
        Get cutout 
        """
        if "r" in self.bands:
            band = "r"
        elif band is None:
            band = self.bands[-1]
            raise UserWarning("Retrieving cutout in {:s} band".format(band))
        
        self.cutout = self.get_image(imsize, band).data
        return self.cutout







        

        


        








    