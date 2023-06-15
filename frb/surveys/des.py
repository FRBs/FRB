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
    _svc = sia.SIAService(defs.NOIR_DEF_ACCESS_URL+'des_dr1')

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
photom['DES']['class_star_r'] = "class_star_r"
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
        self.default_query_fields = list(photom['DES'].values())

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
        main_cat = super(DES_Survey, self).get_catalog(query=query,
                                                       query_fields=query_fields,
                                                       print_query=print_query,**kwargs)
        if len(main_cat) == 0:
            main_cat = catalog_utils.clean_cat(main_cat,photom['DES'])
            return main_cat
        main_cat = catalog_utils.clean_cat(main_cat, photom['DES'])
        ## Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog

