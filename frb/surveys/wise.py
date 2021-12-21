"""WISE Survey"""

import numpy as np

from astropy import units, io, utils
from astropy.table import Table

from frb.surveys import surveycoord
from frb.surveys import catalog_utils
from frb.surveys import defs


from IPython import embed

try:
    from pyvo.dal import TAPService, sia
except ImportError:
    print("Warning: You need to install pyvo to use the survey tools...")
    _svc = None
else:
    _DEF_ACCESS_URL = "https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd"
    _svc = sia.SIAService(_DEF_ACCESS_URL)

# Define the data model for WISE data
photom = {}

# Define the default set of query fields
# See: https://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_1a.html
# for additional Fields
photom['WISE'] = {}
WISE_bands = ['W1', 'W2', 'W3', 'W4']
_DEFAULT_query_fields = ['source_id','ra','dec','tmass_key']
for band in WISE_bands:
    magcol = '{:s}mag'.format(band.lower()) # The "standard" aperture measurement
    errcol = '{:s}sigm'.format(band.lower())
    photom['WISE']['WISE_{:s}'.format(band)] = magcol
    photom['WISE']['WISE_{:s}_err'.format(band)] = errcol
    _DEFAULT_query_fields.append(magcol)
    _DEFAULT_query_fields.append(errcol)
photom['WISE']['ra'] = 'ra'
photom['WISE']['dec'] = 'dec'

class WISE_Survey(surveycoord.SurveyCoord):
    """
    Class to handle queries on the WISE survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate

    """

    def __init__(self, coord, radius, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        self.survey = 'WISE'
        self.bands = WISE_bands
        self.service = TAPService('https://irsa.ipac.caltech.edu/TAP')
        self.query = None
        self.database = "allwise_p3as_psd"

    def get_catalog(self, query=None, query_fields=_DEFAULT_query_fields, 
                    print_query=False, system='AB'):
        """
        Grab a catalog of sources around the input coordinate to the search radius

        Args:
            query: Not used
            query_fields (list, optional): Over-ride list of items to query
            print_query (bool): Print the SQL query generated
            system (str): Magnitude system ['AB', 'Vega']

        Returns:
            astropy.table.Table:  Catalog of sources returned.  Includes WISE
            photometry for matched sources.

            Magnitudes are in AB by default
        """
        # Main WISE query
        if query is None:
            self._gen_cat_query(query_fields)
        if print_query:
            print(self.query)
        main_cat = self.service.run_async(self.query).to_table()
        main_cat.meta['radius'] = self.radius
        main_cat.meta['survey'] = self.survey
        if len(main_cat) == 0:
            main_cat = catalog_utils.clean_cat(main_cat, photom['WISE'], fill_mask=-999.)
            return main_cat
        
        main_cat = catalog_utils.clean_cat(main_cat, photom['WISE'], fill_mask=-999.)

        # Convert to AB mag
        if system == 'AB':
            fnu0 = {'WISE_W1':309.54,
                    'WISE_W2':171.787,
                    'WISE_W3':31.674,
                    'WISE_W4':8.363}
            for item in ['1','2','3','4']:
                filt = 'WISE_W'+item
                main_cat[filt] -= 2.5*np.log10(fnu0[filt]/3631.)
        elif system == 'Vega':
            pass

        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog.copy()
    
    def get_cutout(self, imsize, filter, timeout=120):
        """
        Download an image from IRSA
        Args:
            imsize(Quantity): Size of the cutout in angular units.
            filter(str): One of "W1", "W2", "W3" or "W4"
            timeout(float): Number of seconds to wait to hear a response from
                the IRSA SIA server.
        Returns:
            imghdu(fits.HDU): Fits HDU with image
        """
        assert filter.upper() in self.bands, "Invalid filter name "+filter
        # First get a table with the image coadd_id
        meta_url = "https://irsa.ipac.caltech.edu/ibe/search/wise/allsky/4band_p3am_cdd?POS={:f},{:f}".format(self.coord.ra.value, self.coord.dec.value)
        img_metatab = Table.read(utils.data.download_file(meta_url, cache=True, show_progress=False, timeout=timeout), format="ascii")
        coadd_id = img_metatab['coadd_id'][img_metatab['band']==int(filter[1])][0]

        # Now generate the image url
        img_url = "https://irsa.ipac.caltech.edu/ibe/data/wise/allsky/4band_p3am_cdd/{:s}/{:s}/{:s}/{:s}-{:s}-int-3.fits".format(coadd_id[:2], coadd_id[:4], coadd_id, coadd_id, filter.lower())
        img_url +="?center={:f},{:f}&size={:f}pix".format(self.coord.ra.value, self.coord.dec.value, imsize.to("arcsec").value/2.75)
        self.cutout = io.fits.open(utils.data.download_file(img_url,cache=True,show_progress=False,timeout=timeout))[0]
        self.cutout_size = imsize
        return self.cutout.copy()
        
    def _gen_cat_query(self,query_fields=_DEFAULT_query_fields):
        """
        Generate ADQL query for catalog search

        self.query is modified in place

        Args:
            query_fields (list):  Override the default list for the SQL query

        """
        query_field_str = ""
        for field in query_fields:
            query_field_str += " {:s},".format(field)
        # Remove last comma
        query_field_str = query_field_str[:-1]
        self.query = """SELECT{:s}
        FROM {:s}
        WHERE CONTAINS(POINT('ICRS',ra, dec), CIRCLE('ICRS',{:f},{:f},{:f}))=1""".format(query_field_str,self.database,self.coord.ra.value,
                                self.coord.dec.value,self.radius.to(units.deg).value)
        return self.query