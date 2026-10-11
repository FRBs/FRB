"""DELVE survey"""

from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table

from frb.surveys import dlsurvey, defs
from frb.surveys import catalog_utils

# Dependencies
try:
    from pyvo.dal import sia
except ImportError:
    print("Warning:  You need to install pyvo to retrieve DELVE images")
    _svc = None
else:
    _svc = sia.SIAService("https://datalab.noao.edu/sia/delve_dr2")

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

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['DELVE'] = {'DELVE_ID': int}


class DELVE_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DELVE survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.dlsurvey.DL_Survey`
            (e.g. ``verbose``)

    """

    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DELVE'
        self.bands = DELVE_bands
        self.svc = _svc
        self.qc_profile = "default"
        self.database = "delve_dr2.objects"
        self.default_query_fields = list(photom['DELVE'].values())

    def get_catalog(self, query: str | None = None,
                    query_fields: list[str] | None = None,
                    print_query: bool = False, **kwargs) -> Table:
        """
        Grab a catalog of sources around the input coordinate to the search radius

        Args:
            query (str, optional): SQL query. If None, it is generated
                from the default fields and ``query_fields``.
            query_fields (list of str, optional): Additional items to query
                on top of the default fields
            print_query (bool, optional): Print the SQL query generated
            **kwargs: Passed to :meth:`frb.surveys.dlsurvey.DL_Survey.get_catalog`
                (e.g. ``timeout``)

        Returns:
            astropy.table.Table:  Catalog of sources returned, with the
            columns renamed to the FRB survey names.  Can be empty.
        """
        # Main DES query
        main_cat = super(DELVE_Survey, self).get_catalog(query=query,
                                                         query_fields=query_fields,
                                                         print_query=print_query,**kwargs)
        main_cat = catalog_utils.clean_cat(main_cat, photom['DELVE'], mask_photometry=True)
        # Empty catalogs get the standard columns (no-op otherwise)
        main_cat = catalog_utils.ensure_empty_schema(main_cat, list(photom['DELVE'].keys()),
                                                     dtypes=schema_dtypes['DELVE'])
        
        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog

