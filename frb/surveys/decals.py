"""DECaLS"""

import warnings

import numpy as np
from astropy.coordinates import Angle, SkyCoord
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

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['DECaL'] = {'DECaL_ID': int, 'DECaL_brick': int, 'DECaL_type': str}

class DECaL_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the DECaL survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.dlsurvey.DL_Survey`
            (e.g. ``verbose``)

    """

    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'DECaL'
        self.bands = ['g', 'r', 'z']
        self.svc = _svc # sia.SIAService("https://datalab.noao.edu/sia/ls_dr8")
        self.qc_profile = "default"
        self.database = "ls_dr10.tractor"
        self.default_query_fields = list(photom['DECaL'].values())

    def get_catalog(self, query: str | None = None,
                    query_fields: list[str] | None = None,
                    print_query: bool = False, exclude_stars: bool = False,
                    **kwargs) -> Table:
        """
        Grab a catalog of sources around the input coordinate to the search radius

        Args:
            query (str, optional): SQL query. If None, it is generated
                from the default fields and ``query_fields``, and joined
                with the DECaLS photo-z table.
            query_fields (list of str, optional): Additional items to query
                on top of the default fields
            print_query (bool, optional): Print the SQL query generated
            exclude_stars (bool, optional): If the field 'type' is present and is 'PSF',
                remove those objects from the output catalog.
            **kwargs: Passed to :meth:`frb.surveys.dlsurvey.DL_Survey.get_catalog`
                (e.g. ``timeout``)

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
            main_cat = catalog_utils.clean_cat(main_cat, photom['DECaL'], mask_photometry=True)
            self.catalog = catalog_utils.ensure_empty_schema(main_cat, list(photom['DECaL'].keys()),
                                                             dtypes=schema_dtypes['DECaL'])
            self.validate_catalog()
            return self.catalog
        #
        for col in main_cat.colnames:
            # Skip strings
            if main_cat[col].dtype not in [float, int]:
                continue
            else:
                main_cat[col].mask = np.isnan(main_cat[col])
        
        #Convert SNR to mag error values.
        snr_cols = [colname for colname in main_cat.colnames if "snr" in colname]
        for col in snr_cols:
            main_cat[col].mask = main_cat[col]<0
            main_cat[col] = 2.5*np.log10(1+1/main_cat[col])
        # Clean
        main_cat = catalog_utils.clean_cat(main_cat, photom['DECaL'], mask_photometry=True)
        # Remove point sources if necessary
        if exclude_stars and 'DECaL_type' in main_cat.colnames:
            main_cat = main_cat[main_cat['DECaL_type'] != 'PSF']
        elif exclude_stars:
            warnings.warn("'DECaL_type' not found in catalog, cannot exclude stars.")
        self.catalog = main_cat
        self.validate_catalog()
        # Return
        return self.catalog

    def _parse_cat_band(self, band: str) -> tuple[list[str], list[str], str]:
        """
        Internal method to generate the bands for grabbing
        a cutout image

        Args:
            band (str): Band desired; 'g', 'r' or 'z'

        Returns:
            tuple: Three items:

                - list of str: Names of the image table columns to select on
                - list of str: Values the columns must have
                - str: Band string for cutout

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