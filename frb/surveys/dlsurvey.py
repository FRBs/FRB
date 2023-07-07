"""
DataLab survey class. Gets data from any survey
available through the NOAO datalab-client.
"""
import numpy as np
import warnings
from astropy import units, io, utils
import warnings

from frb.surveys import catalog_utils

# Dependencies
try:
    from dl import queryClient as qc, authClient as ac
    from dl.helpers.utils import convert
except:
    print("Warning:  datalab-client is not installed or will not properly connect")


try:
    from pyvo.dal import DALFormatError
except ImportError:
    print("Warning:  You need to install pyvo to retrieve images")
    _svc = None

from frb.surveys import surveycoord

class DL_Survey(surveycoord.SurveyCoord):
    """
    A survey class for all databases hosted
    by NOIR's DataLab. Inherits from SurveyCoord
    """
    def __init__(self, coord, radius, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        
        #Define photmetric band names.
        self.token = ac.login('anonymous')
        self.bands = []
        #Instantiate sia service 
        self.svc = None
        #Generate query
        self.query = None
        self.qc_profile = None
        self.default_query_fields = None
    
    def _parse_cat_band(self,band):
        """
        Internal method to generate the bands for grabbing
        a cutout image

        For most child classes, nothing much is necessary. Only
        gets modified in DECaLS. 

        Args:
            band (str): Band desired

        Returns:
            list, list, str:  Table columns, Column values, band string for cutout

        """
        table_cols = ['proctype','prodtype']
        col_vals = ['Stack','image']

        return table_cols, col_vals, band

    def _gen_cat_query(self,query_fields=None, qtype='main'):
        """
        Generate SQL Query for catalog search

        self.query is modified in place

        Args:
            query_fields (list):  Override the default list for the SQL query

        """
        if self.default_query_fields is None:
            raise IOError("DLSurvey child incorrectly instantiated.  Missing default_query_fields")
        if query_fields is None:
            # Main query
            if qtype == 'main':
                query_fields = self.default_query_fields
                database = self.database
            else:
                raise IOError("Bad qtype")
        else:
            if qtype == 'main':
                assert isinstance(query_fields, list), "query_fields must be a list"
                query_fields = np.union1d(self.default_query_fields, query_fields)
                database = self.database
            else:
                raise IOError("Bad qtype")

        self.query = _default_query_str(query_fields, database,self.coord,self.radius)
        # Return
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
        # Get one with maximum zero point.
        row = imgTable[np.argmax(imgTable['magzero'].data.data.astype('float'))] # pick image with longest exposure time
        url = row['access_url']
        if verbose:
            print ('downloading deepest stacked image...')

        imagedat = io.fits.open(utils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))
        return imagedat

    def get_catalog(self, query=None, query_fields=None, print_query=False,timeout=120, photomdict=None):
        """
        Get catalog sources around the given coordinates
        within self.radius.
        
        Args:
            query (str, optional): SQL query to generate the catalog
            query_fields (list, optional): Over-ride list of items to query
            print_query (bool): Print the SQL query generated 
        
        Returns:
            astropy.table.Table:  Catalog of sources obtained from the SQL query.
        """
        qc.set_profile(self.qc_profile)
        # Generate the query
        if query is None:
            self._gen_cat_query(query_fields)
            query = self.query
        if print_query:
            print(query)
        # Do it while silencing print statements
        result = qc.query(self.token, sql=query,timeout=timeout)
        self.catalog = convert(result,outfmt="table")

        if photomdict:
            self.catalog = catalog_utils.clean_cat(self.catalog, photomdict)
        
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey
        # Validate
        self.validate_catalog()
        # Return
        return self.catalog.copy()
    
    def get_image(self, imsize, band, timeout=120, verbose=False):
        """
        Get images from the catalog if available
            for a given fov and band.

        Args:
            imsize (Quantity): FOV for the desired image
            band (str): Band for the image (e.g. 'r')
            timeout (int, optional): Time to wait in seconds before timing out
            verbose (bool, optional):

        Returns:
            HDU: Image header data unit

        """
        if self.svc is None:
            raise RuntimeError("svc attribute cannot be None. Have you installed pyvo?")
        
        if band.lower() not in self.bands and band not in self.bands:
            raise TypeError("Allowed filters (case-insensitive) for {:s} photometric bands are {}".format(self.survey,self.bands))

        table_cols, col_vals, bandstr = self._parse_cat_band(band)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        try:
            imgTable = self.svc.search(self.coord, imsize, verbosity=2).to_table()
        except DALFormatError:
            warnings.warn_explicit(f"Image cannot be retrieved. Invalid base URL?: {self.svc._baseurl}.",
                               category=RuntimeWarning, filename="FRB/frb/surveys/dlsurvey.py", lineno=114)
            return None
        if verbose:
            print("The full image list contains", len(imgTable), "entries")
        
        #Select band
        selection = imgTable['obs_bandpass'].astype(str)==bandstr

        #from IPython import embed; embed(header='117')
        #Select images in that band
        for column, value in zip(table_cols,col_vals):
            selection = selection & ((imgTable[column].astype(str)==value))
        imgTable = imgTable[selection]

        if(len(imgTable)>0):
            imagedat = self._select_best_img(imgTable,verbose=True,timeout=timeout)
            img_hdu = imagedat[0]
        else:
            print('No image available')
            img_hdu = None
        return img_hdu
    
    def get_cutout(self, imsize, band=None):
        """
        Get cutout (and header)

        Args:
            imsize (Quantity): e.g 10*units.arcsec
            band (str): e.g. 'r'

        Returns:
            ndarray, Header: cutout image, cutout image header

        """
        self.cutout_size = imsize

        if band is None:
            if "r" in self.bands:
                band = "r"
            elif band is None:
                band = self.bands[-1]
                warnings.warn("Retrieving cutout in {:s} band".format(band))

        img_hdu = self.get_image(imsize, band)
        if img_hdu is not None:
            self.cutout = img_hdu.data
            self.cutout_hdr = img_hdu.header
        else:
            self.cutout = None
            self.cutout_hdr = None
        return self.cutout, self.cutout_hdr


def _default_query_str(query_fields, database, coord, radius):
    """
    Generates default query string for a catalog search.
    
    Args:
        query_fields (list of str): A list of query fields to
            retrieve from the database
        database (str): Name of the database
        coord (astropy.coordinates.SkyCoord): Central coordinate of the search
        radius (astropy.units.Quantity or Angle): Search radius
        
    Returns:
        str: A query to be fed to datalab's SQL client
        
    """
    query_field_str = ""
    for field in query_fields:
        query_field_str += " {:s},".format(field)
    # Remove last comma
    query_field_str = query_field_str[:-1]
    default_query = """SELECT{:s}
    FROM {:s}
    WHERE q3c_radial_query(ra,dec,{:f},{:f},{:f})
    """.format(query_field_str,database,coord.ra.value,
                            coord.dec.value,radius.to(units.deg).value)
    return default_query
