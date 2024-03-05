"""DELVE survey"""

import numpy as np
from astropy import units, io, utils

from frb.surveys import dlsurvey, defs
from frb.surveys import catalog_utils

# Dependencies
try:
    from pyvo.dal import sia
except ImportError:
    print("Warning:  You need to install pyvo to retrieve DELVE images")
    _svc = None
else:
    _svc = sia.SIAService(defs.NOIR_DEF_ACCESS_URL+'delve_dr1')

# Define the data model for DELVE data
# See https://datalab.noirlab.edu/query.php?name=delve_dr2.objects for
# table schema
photom = {}
photom['DELVE'] = {}
photom['DELVE']['DELVE_ID'] = 'quick_object_id'
photom['DELVE']['ra'] = 'ra'
photom['DELVE']['dec'] = 'dec'
photom['DELVE']['ebv'] = 'ebv' # Schegel, Finkbeiner, Davis (1998)
DELVE_bands = ['g', 'r', 'i', 'z']
for band in DELVE_bands:
    photom['DELVE'][f'DELVE_{band}'] = f'mag_auto_{band}' #mag
    photom['DELVE'][f'DELVE_{band}_err'] = f'magerr_auto_{band}' #magerr
    photom['DELVE'][f'class_star_{band}'] = f'class_star_{band}' #morphology class


class DELVE_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DELVE survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate

    """

    def __init__(self, coord, radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DELVE'
        self.bands = DELVE_bands
        self.svc = sia.SIAService("https://datalab.noao.edu/sia/delve_dr2")
        self.qc_profile = "default"
        self.database = "delve_dr2.objects"
        self.default_query_fields = list(photom['DELVE'].values())

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
        main_cat = super(DELVE_Survey, self).get_catalog(query=query,
                                                         query_fields=query_fields,
                                                         print_query=print_query,**kwargs)
        if len(main_cat) == 0:
            main_cat = catalog_utils.clean_cat(main_cat,photom['DELVE'])
            return main_cat
        main_cat = catalog_utils.clean_cat(main_cat, photom['DELVE'])
        #import pdb; pdb.set_trace()
        for col in main_cat.colnames:
            if main_cat[col].dtype in [int, float]:
                mask = np.isnan(main_cat[col])+(main_cat[col]==99.99)+(main_cat[col]==99)
                main_cat[col] = np.where(~mask, main_cat[col], -999.0)
        
        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog

