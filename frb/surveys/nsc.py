"""NOIRLab source catalog"""

from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table

from frb.surveys import dlsurvey, defs
from frb.surveys import catalog_utils

# Dependencies
try:
    from pyvo.dal import sia
except ImportError:
    print("Warning:  You need to install pyvo to retrieve DES images")
    _svc = None
else:
    _svc = sia.SIAService(defs.NOIR_DEF_ACCESS_URL+'nsa')

# Define the data model for DES data
photom = {}
photom['NSC'] = {}
photom['NSC']['NSC_ID'] = 'id'
photom['NSC']['ra'] = 'ra'
photom['NSC']['dec'] = 'dec'
photom['NSC']['class_star'] = 'class_star'
NSC_bands = ['u','g', 'r', 'i', 'z', 'Y', 'VR']
for band in NSC_bands:
    photom['NSC']['NSC_{:s}'.format(band)] = '{:s}mag'.format(band.lower())
    photom['NSC']['NSC_{:s}_err'.format(band)] = '{:s}rms'.format(band.lower())

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['NSC'] = {'NSC_ID': str}

class NSC_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the NSC survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.dlsurvey.DL_Survey`
            (e.g. ``verbose``)

    """

    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'NSC'
        self.bands = NSC_bands
        self.svc = _svc
        self.qc_profile = "default"
        self.database = "nsc_dr2.object"
        self.default_query_fields = list(photom['NSC'].values())

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
        main_cat = super(NSC_Survey, self).get_catalog(query=query,
                                                       query_fields=query_fields,
                                                       print_query=print_query,**kwargs)
        if len(main_cat) == 0:
            main_cat = catalog_utils.clean_cat(main_cat, photom['NSC'], mask_photometry=True)
            main_cat = catalog_utils.ensure_empty_schema(main_cat, list(photom['NSC'].keys()),
                                                       dtypes=schema_dtypes['NSC'])
            return main_cat
        main_cat = catalog_utils.clean_cat(main_cat, photom['NSC'], mask_photometry=True)
        
        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog

