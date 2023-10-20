"""DECaLS"""

import numpy as np
from astropy import units, io, utils
from astropy.table import Table

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils
from frb.surveys import defs

import pandas as pd

# Dependencies
try:
    from pyvo.dal import sia
except ImportError:
    print("Warning:  You need to install pyvo to retrieve DECaL images")
    _svc = None
else:
    _svc = sia.SIAService(defs.NOIR_DEF_ACCESS_URL+'ls_dr8')

# Define the Photometric data model for DECaL
photom = {}
photom['DECaL'] = {}
DECaL_bands = ['g', 'r', 'z']
for band in DECaL_bands:
    if "W" not in band:
        bandstr = 'DECaL_'+band
    else:
        bandstr = 'WISE_'+band
    photom['DECaL'][bandstr] = 'mag_{:s}'.format(band.lower())
    photom['DECaL'][bandstr+"_err"] = 'snr_{:s}'.format(band.lower())
photom['DECaL']['DECaL_ID'] = 'ls_id'
photom['DECaL']['ra'] = 'ra'
photom['DECaL']['dec'] = 'dec'
photom['DECaL']['DECaL_brick'] = 'brickid'
photom['DECaL']['gaia_pointsource'] = 'gaia_pointsource'

class DECaL_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DECaL survey
    
    Child of DL_Survey which uses datalab to access NOAO
    
    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate
        
    """

    def __init__(self, coord, radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DECaL'
        self.bands = ['g', 'r', 'z']
        self.svc = _svc # sia.SIAService("https://datalab.noao.edu/sia/ls_dr7")
        self.qc_profile = "default"
        self.database = "ls_dr8.tractor"
        self.default_query_fields = list(photom['DECaL'].values())

    def get_catalog(self, query=None, query_fields=None, print_query=False,exclude_gaia=False,**kwargs):
        """
        Grab a catalog of sources around the input coordinate to the search radius
        
        Args:
            query: SQL query
            query_fields (list, optional): Over-ride list of items to query
            exclude_gaia (bool,optional): If the field 'gaia_pointsource' is present and is 1,
                                         remove those objects from the output catalog.
            print_query (bool): Print the SQL query generated 

        Returns:
            astropy.table.Table:  Catalog of sources returned

        """
        # Query
        main_cat = super(DECaL_Survey, self).get_catalog(query=query,
                                                         query_fields=query_fields,
                                                         print_query=print_query,**kwargs)
        main_cat = Table(main_cat,masked=True)
        #
        for col in main_cat.colnames:
            main_cat[col].mask = pd.isnull(main_cat[col])
            # main_cat[col].mask = np.isnan(main_cat[col])
        #Convert SNR to mag error values.
        snr_cols = [colname for colname in main_cat.colnames if "snr" in colname]
        for col in snr_cols:
            main_cat[col].mask = main_cat[col]<0
            main_cat[col] = 2.5*np.log10(1+1/main_cat[col])
        
        main_cat = main_cat.filled(-99.0)
        #Remove gaia objects if necessary
        if exclude_gaia:
            self.catalog = main_cat[main_cat['gaia_pointsource']==0]
        else:
            self.catalog = main_cat
        # Clean
        main_cat = catalog_utils.clean_cat(main_cat, photom['DECaL'])
        self.validate_catalog()
        # Return
        return self.catalog

    def _parse_cat_band(self, band):
        """
        Internal method to generate the bands for grabbing
        a cutout image
        
        Args:
            band (str): Band desired 

        Returns:
            list, list, str:  Table columns, Column values, band string for cutout

        """
        if band == 'g':
            bandstr = "g DECam SDSS c0001 4720.0 1520.0"
        elif band == 'r':
            bandstr = "r DECam SDSS c0002 6415.0 1480.0"
        elif band == 'z':
            bandstr = "z DECam SDSS c0004 9260.0 1520.0"
        table_cols = ['prodtype']
        col_vals = ['image']
        return table_cols, col_vals, bandstr
