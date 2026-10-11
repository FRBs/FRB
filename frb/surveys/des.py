"""DES Survey"""

from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table

from frb.surveys import dlsurvey
from frb.surveys import catalog_utils
from frb.surveys import defs

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

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['DES'] = {'DES_ID': int, 'DES_tile': str}


class DES_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DES survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.dlsurvey.DL_Survey`
            (e.g. ``verbose``)

    """

    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DES'
        self.bands = ['g', 'r', 'i', 'z', 'y']
        self.svc = _svc
        self.qc_profile = "default"
        self.database = "des_dr2.main"
        self.default_query_fields = list(photom['DES'].values())

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
        main_cat = super(DES_Survey, self).get_catalog(query=query,
                                                       query_fields=query_fields,
                                                       print_query=print_query,**kwargs)
        if len(main_cat) == 0:
            main_cat = catalog_utils.clean_cat(main_cat, photom['DES'], mask_photometry=True)
            main_cat = catalog_utils.ensure_empty_schema(main_cat, list(photom['DES'].keys()),
                                                       dtypes=schema_dtypes['DES'])
            return main_cat
        main_cat = catalog_utils.clean_cat(main_cat, photom['DES'], mask_photometry=True)
        ## Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog

