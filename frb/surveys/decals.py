"""DECaLS"""

import numpy as np
from astropy import units, io, utils
from astropy.table import Table

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils
from frb.surveys import defs

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
photom['DECaL']['DECaL_type'] = 'type' # Replaces `gaia_pointsource` from DR8.

class DECaL_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DECaL survey
    
    Child of DL_Survey which uses datalab to access NOAO
    
    Args:
        coord (SkyCoord): Coordinate for surveying around
        radius (Angle): Search radius around the coordinate
        
    """

    def __init__(self, coord, radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DECaL'
        self.bands = ['g', 'r', 'z']
        self.svc = _svc # sia.SIAService("https://datalab.noao.edu/sia/ls_dr8")
        self.qc_profile = "default"
        self.database = "ls_dr10.tractor"
        self.default_query_fields = list(photom['DECaL'].values())

    def get_catalog(self, query=None, query_fields=None, print_query=False,exclude_stars=False,**kwargs):
        """
        Grab a catalog of sources around the input coordinate to the search radius
        
        Args:
            query: SQL query
            query_fields (list, optional): Over-ride list of items to query
            exclude_gaia (bool,optional): If the field 'type' is present and is 'PSF',
                                         remove those objects from the output catalog.
            print_query (bool): Print the SQL query generated 

        Returns:
            astropy.table.Table:  Catalog of sources returned
                Can be empty

        """
        # Query
        if query is None:
            query = super(DECaL_Survey, self)._gen_cat_query(query_fields=query_fields, qtype='main')
            # include photo_z
            query = query.replace("SELECT", "SELECT z_phot_median, z_spec, survey, z_phot_l68, z_phot_u68, z_phot_l95, z_phot_u95,")
            query = query.replace("ls_id", "t.ls_id")
            query = query.replace("brickid", "t.brickid")
            query = query.replace(f"FROM {self.database}\n", f"FROM {self.database} as t LEFT JOIN {self.database.split('.')[0]}.photo_z AS p ON t.ls_id=p.ls_id\n")
        self.query = query
        main_cat = super(DECaL_Survey, self).get_catalog(query=self.query,
                                                        print_query=print_query,**kwargs)
        main_cat = Table(main_cat,masked=True)
        if len(main_cat)==0:
            return main_cat 
        #
        for col in main_cat.colnames:
            # Skip strings
            if main_cat[col].dtype not in [float, int]:
                continue
            else:
                try:
                    main_cat[col].mask = np.isnan(main_cat[col])
                except:
                    import pdb; pdb.set_trace()
        
        #Convert SNR to mag error values.
        snr_cols = [colname for colname in main_cat.colnames if "snr" in colname]
        for col in snr_cols:
            main_cat[col].mask = main_cat[col]<0
            main_cat[col] = 2.5*np.log10(1+1/main_cat[col])
        
        main_cat = main_cat.filled(-99.0)
        #Remove gaia objects if necessary
        if exclude_stars and 'type' in main_cat.colnames:
            self.catalog = main_cat[main_cat['DECaL_type']=='PSF']
        elif exclude_stars and 'type' not in main_cat.colnames:
            print("Warning: 'type' not found in catalog, cannot exclude stars.")
            self.catalog = main_cat
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