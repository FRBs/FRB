"""DESI"""

from typing import NoReturn

import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils

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

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['DESI'] = {'DESI_ID': int, 'DESI_name': str, 'DESI_spectype': str,
                         'DESI_specsubtype': str, 'DESI_survey': str,
                         'DESI_z_warn': int, 'DESI_zcat_primary': str,
                         'DESI_zcat_nspec': int}



class DESI_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DESI survey
    
    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.dlsurvey.DL_Survey`
            (e.g. ``verbose``)

    """

    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DESI'
        self.qc_profile = "default"
        self.database = "desi_dr1.zpix"
        self.default_query_fields = list(spectrom['DESI'].values())

    def get_catalog(self, query: str | None = None,
                    query_fields: list[str] | None = None,
                    print_query: bool = False, exclude_stars: bool = False,
                    zcat_primary_only: bool = True, **kwargs) -> Table:
        """
        Grab a catalog of sources around the input coordinate to the search radius

        Args:
            query (str, optional): SQL query. If None, it is generated
                from the default fields and ``query_fields``.
            query_fields (list of str, optional): Additional items to query
                on top of the default fields
            print_query (bool, optional): Print the SQL query generated
            exclude_stars (bool, optional): If the field 'spectype' is present and is 'STAR',
                remove those objects from the output catalog.
            zcat_primary_only (bool, optional): If True, only return objects with zcat_primary=True
            **kwargs: Passed to :meth:`frb.surveys.dlsurvey.DL_Survey.get_catalog`
                (e.g. ``timeout``)

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
            main_cat = catalog_utils.clean_cat(main_cat, spectrom['DESI'])
            main_cat = catalog_utils.ensure_empty_schema(main_cat, list(spectrom['DESI'].keys()),
                                                       dtypes=schema_dtypes['DESI'])
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
        
        main_cat = catalog_utils.fill_masked(main_cat, -99.0)
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

        # Normalize and validate
        self.validate_catalog()
        
        # Return
        return self.catalog

    def get_image(self, **kwargs) -> NoReturn:
        """
        Not implemented for DESI, which only provides spectroscopic data.

        Args:
            **kwargs: Ignored.

        Raises:
            NotImplementedError: Always.

        """
        raise NotImplementedError("Cutout retrieval not implemented for DESI. This class is meant to purely retreive spectroscopic data. For imaging, use the DeCAL_Survey class instead.")