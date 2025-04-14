"""DESI"""

import numpy as np
from astropy import units, io, utils
from astropy.table import Table

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils
from frb.surveys import defs

# Define the spectrometric data model for DESI
spectrom = {}
spectrom['DESI'] = {}
spectrom['DESI']['DESI_ID'] = 'targetid'
spectrom['DESI']['ra'] = 'mean_fiber_ra'
spectrom['DESI']['dec'] = 'mean_fiber_dec'
spectrom['DESI']['DESI_name'] = 'desiname' # This is Just JXXXXXX.XX+YYYYYY.YY at J2000 epoch.
spectrom['DESI']['DESI_spectype'] = 'spectype' # STAR, GALAXY or QSO
spectrom['DESI']['DESI_specsubtype'] = 'subtype' # Futher classification of spectype.
spectrom['DESI']['DESI_survey'] = 'survey' # BGS, LRG or ELG.
spectrom['DESI']['DESI_z'] = 'z' # redshift
spectrom['DESI']['DESI_z_err'] = 'zerr' # redshift error
spectrom['DESI']['DESI_z_warn'] = 'zwarn' # redrock warning flags? Need to look up what they signify.
spectrom['DESI']['DESI_zcat_primary'] = 'zcat_primary' # In case there are multiple entries with this object, use this bool to choose the "preferred" z.
spectrom['DESI']['DESI_zcat_nspec'] = 'zcat_nspec' # Number of coadded spectra.



class DESI_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DESI survey
    
    Child of DL_Survey which uses datalab to access NOAO
    
    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate
        
    """

    def __init__(self, coord, radius, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DESI'
        self.qc_profile = "default"
        self.database = "desi_dr1.zpix"
        self.default_query_fields = list(spectrom['DESI'].values())

    def get_catalog(self, query=None, query_fields=None, print_query=False, 
                    exclude_stars=False, zcat_primary_only=True,**kwargs):
        """
        Grab a catalog of sources around the input coordinate to the search radius
        
        Args:
            query: SQL query
            query_fields (list, optional): Over-ride list of items to query
            exclude_stars (bool,optional): If the field 'spectype' is present and is 'STAR',
                                         remove those objects from the output catalog.
            print_query (bool): Print the SQL query generated
            zcat_primary_only (bool): If True, only return objects with zcat_primary=True

        Returns:
            astropy.table.Table:  Catalog of sources returned
                Can be empty

        """
        # Query
        if query is None:
            query = super(DESI_Survey, self)._gen_cat_query(query_fields=query_fields,qtype='main',
                                                            ra_col = spectrom['DESI']['ra'],
                                                            dec_col = spectrom['DESI']['dec'])
        self.query = query
        main_cat = super(DESI_Survey, self).get_catalog(query=self.query,
                                                        print_query=print_query, photomdict=spectrom['DESI'],**kwargs)
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
        
        main_cat = main_cat.filled(-99.0)
        #Remove gaia objects if necessary
        if exclude_stars and 'DESI_spectype' in main_cat.colnames:
            main_cat = main_cat[main_cat['DESI_spectype']!='STAR']
        elif exclude_stars and 'spectype' not in main_cat.colnames:
            print("Warning: 'DESI_spectype' not found in catalog, cannot exclude stars.")

        # Clean
        self.catalog = catalog_utils.clean_cat(main_cat, spectrom['DESI'], fill_mask=-99.0)

        # Only zcat_primary?
        if zcat_primary_only and 'DESI_zcat_primary' in self.catalog.colnames:
            self.catalog = self.catalog[self.catalog['DESI_zcat_primary'] == 't']
        elif zcat_primary_only and 'DESI_zcat_primary' not in self.catalog.colnames:
            print("Warning: 'DESI_zcat_primary' not stored in catalog, cannot filter by zcat_primary.")

        # Return
        return self.catalog