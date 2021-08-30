"""WISE Survey"""

from IPython import embed

import numpy as np
from astropy import units, io, utils

from frb.surveys import surveycoord
from frb.surveys import catalog_utils
from frb.surveys import images
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.table import Table

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
# See: hhttps://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_1a.html
# for additional Fields
photom['WISE'] = {}
WISE_bands = ['W1', 'W2', 'W3', 'W4']
_DEFAULT_query_fields = ['source_id','ra','dec','tmass_key']
for band in WISE_bands:
    magcol = '{:s}mag'.format(band.lower())
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

    def get_catalog(self, query=None, query_fields=_DEFAULT_query_fields, print_query=False):
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
        # Main WISE query
        if query is None:
            self._gen_cat_query(query_fields)
        if print_query:
            print(self.query)
        main_cat = self.service.run_async(self.query).to_table()
        main_cat.meta['radius'] = self.radius
        main_cat.meta['survey'] = self.survey
        if len(main_cat) == 0:
            return main_cat
        main_cat = catalog_utils.clean_cat(self.catalog, photom['WISE'], fill_mask=-999.)

        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog.copy()
    
    def get_image(self, imsize, filter):
        return
        

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

    def _select_best_img(self,imgTable,verbose,timeout=120):
        """
        Select the best band for a cutout

        Args:
            imgTable: Table of images
            verbose (bool):  Print status
            timeout (int or float):  How long to wait before timing out, in seconds

        Returns:
            HDU: header data unit for the downloaded image

        """
        row = imgTable[np.argmax(imgTable['exptime'].data.data.astype('float'))] # pick image with longest exposure time
        url = row['access_url'].decode()
        if verbose:
            print ('downloading deepest stacked image...')

        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat

def get_url(coord, imsize=30., scale=0.396127, grid=False, label=False, invert=False):
    """
    Generate the SDSS URL for an image retrieval

    Args:
        coord (astropy.coordiantes.SkyCoord): Center of image
        imsize: float, optional
          Image size (rectangular) in arcsec and without units
        scale (float, optional):
        grid (bool, optional):
        label (bool, optional):
        invert (bool, optional):

    Returns:
        str:  URL for the image

    """

    # Pixels
    npix = round(imsize/scale)
    xs = npix
    ys = npix

    # Generate the http call
    #name1='http://skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/'
    name1 = 'http://skyservice.pha.jhu.edu/DR12/ImgCutout/'
    name='getjpeg.aspx?ra='

    name+=str(coord.ra.value) 	#setting the ra (deg)
    name+='&dec='
    name+=str(coord.dec.value)	#setting the declination
    name+='&scale='
    name+=str(scale) #setting the scale
    name+='&width='
    name+=str(int(xs))	#setting the width
    name+='&height='
    name+=str(int(ys)) 	#setting the height

    #------ Options
    options = ''
    if grid is True:
        options+='G'
    if label is True:
        options+='L'
    if invert is True:
        options+='I'
    if len(options) > 0:
        name+='&opt='+options

    name+='&query='

    url = name1+name
    return url