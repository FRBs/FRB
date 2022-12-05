"""
Slurp data from Pan-STARRS catalog using the MAST API.
A lot of this code has been directly taken from
http://ps1images.stsci.edu/ps1_dr2_api.html

"""

import numpy as np

from astropy import units as u,utils as astroutils
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table
from ..galaxies.defs import PanSTARRS_bands

from .images import grab_from_url

import warnings
import requests

from frb.surveys import surveycoord,catalog_utils,images

from IPython import embed

#TODO: It's potentially viable to use the same code for other
#catalogs in the VizieR database. Maybe a generalization wouldn't
#be too bad in the future.

# Define the data model for Pan-STARRS data
photom = {}
photom['Pan-STARRS'] = {}
for band in PanSTARRS_bands:
    # Pre 180301 paper
    #photom["Pan-STARRS"]["Pan-STARRS"+'_{:s}'.format(band)] = '{:s}PSFmag'.format(band.lower())
    #photom["Pan-STARRS"]["Pan-STARRS"+'_{:s}_err'.format(band)] = '{:s}PSFmagErr'.format(band.lower())
    photom["Pan-STARRS"]["Pan-STARRS"+'_{:s}'.format(band)] = '{:s}KronMag'.format(band.lower())
    photom["Pan-STARRS"]["Pan-STARRS"+'_{:s}_err'.format(band)] = '{:s}KronMagErr'.format(band.lower())
    photom["Pan-STARRS"]["Pan-STARRS_ID"] = 'objID'
photom["Pan-STARRS"]['ra'] = 'raStack'
photom["Pan-STARRS"]['dec'] = 'decStack'
photom["Pan-STARRS"]["Pan-STARRS_field"] = 'field'

# Define the default set of query fields
# See: https://outerspace.stsci.edu/display/PANSTARRS/PS1+StackObjectView+table+fields
# for additional Fields
_DEFAULT_query_fields = ['objID','raStack','decStack','objInfoFlag','qualityFlag', 
                         'rKronRad']#, 'rPSFMag', 'rKronMag']
_DEFAULT_query_fields +=['{:s}PSFmag'.format(band) for band in PanSTARRS_bands]
_DEFAULT_query_fields +=['{:s}PSFmagErr'.format(band) for band in PanSTARRS_bands]
_DEFAULT_query_fields +=['{:s}KronMag'.format(band) for band in PanSTARRS_bands]
_DEFAULT_query_fields +=['{:s}KronMagErr'.format(band) for band in PanSTARRS_bands]

class Pan_STARRS_Survey(surveycoord.SurveyCoord):
    """
    A class to access all the catalogs hosted on the
    MAST database. Inherits from SurveyCoord. This
    is a super class not meant for use by itself and
    instead meant to instantiate specific children
    classes like PAN-STARRS_Survey
    """
    def __init__(self,coord,radius,**kwargs):
        surveycoord.SurveyCoord.__init__(self,coord,radius,**kwargs)

        self.Survey = "Pan_STARRS"
    
    def get_catalog(self,query_fields=None,release="dr2",
                    table="stack",print_query=False,
                    use_psf=False):
        """
        Query a catalog in the MAST Pan-STARRS database for
        photometry.

        Args:
            query_fields: list, optional
                A list of query fields to
                get in addition to the
                default fields.
            release: str, optional
                "dr1" or "dr2" (default: "dr2").
                Data release version.
            table: str, optional
                "mean","stack" or "detection"
                (default: "stack"). The data table to
                search within.
            use_psf: bool, optional
                If True, use PSFmag instead of KronMag
        
        Returns:
            catalog: astropy.table.Table
                Contains all query results
        """
        #assert self.radius <= 0.5*u.deg, "Cone serches have a maximum radius"
        #Validate table and release input
        _check_legal(table,release)
        url = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/{:s}/{:s}.csv".format(release,table)
        if query_fields is None:
            query_fields = _DEFAULT_query_fields
        else:
            query_fields = _DEFAULT_query_fields+query_fields
        
        #Validate columns
        _check_columns(query_fields,table,release)
        data = {}
        data['ra'] = self.coord.ra.value
        data['dec'] = self.coord.dec.value
        data['radius'] = self.radius.to(u.deg).value
        data['columns'] = query_fields
        if print_query:
            print(url)
        ret = requests.get(url,params=data)
        ret.raise_for_status()
        if len(ret.text)==0:
            self.catalog = Table()
            self.catalog.meta['radius'] = self.radius
            self.catalog.meta['survey'] = self.survey
            # Validate
            self.validate_catalog()
            return self.catalog.copy()
        photom_catalog = Table.read(ret.text,format="ascii.csv")
        pdict = photom['Pan-STARRS'].copy()

        # Allow for PSF
        if use_psf:
            for band in PanSTARRS_bands:
                pdict["Pan-STARRS"+'_{:s}'.format(band)] = '{:s}PSFmag'.format(band.lower())
                pdict["Pan-STARRS"+'_{:s}_err'.format(band)] = '{:s}PSFmagErr'.format(band.lower())
        
        photom_catalog = catalog_utils.clean_cat(photom_catalog,pdict)

        #Remove bad positions because Pan-STARRS apparently decided
        #to flag some positions with large negative numbers. Why even keep
        #them?
        bad_ra = (photom_catalog['ra']<0)+(photom_catalog['ra']>360)
        bad_dec = (photom_catalog['dec']<-90)+(photom_catalog['dec']>90)
        bad_pos = bad_ra+bad_dec # bad_ra OR bad_dec
        photom_catalog = photom_catalog[~bad_pos]
        
        # Remove duplicate entries.
        photom_catalog = catalog_utils.remove_duplicates(photom_catalog, "Pan-STARRS_ID")

        
        self.catalog = catalog_utils.sort_by_separation(photom_catalog, self.coord,
                                                        radec=('ra','dec'), add_sep=True)
        # Meta
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        #Validate
        self.validate_catalog()

        #Return
        return self.catalog.copy()

    def get_cutout(self,imsize=30*u.arcsec,filt="irg",output_size=None):
        """
        Grab a color cutout (PNG) from Pan-STARRS

        Args:
            imsize (Quantity):  Angular size of image desired
            filt (str): A string with the three filters to be used
            output_size (int): Output image size in pixels. Defaults
                                to the original cutout size.
        Returns:
            PNG image, None (None for the header).
        """
        assert len(filt)==3, "Need three filters for a cutout."
        #Sort filters from red to blue
        filt = filt.lower() #Just in case the user is cheeky about the filter case.
        reffilt = "yzirg"
        idx = np.argsort([reffilt.find(f) for f in filt])
        newfilt = ""
        for i in idx:
            newfilt += filt[i]
        #Get image url
        url = _get_url(self.coord,imsize=imsize,filt=newfilt,output_size=output_size,color=True,imgformat='png')
        self.cutout = images.grab_from_url(url)
        self.cutout_size = imsize
        return  self.cutout.copy(), 
    
    def get_image(self,imsize=30*u.arcsec,filt="i",timeout=120):
        """
        Grab a fits image from Pan-STARRS in a
        specific band.

        Args:
            imsize (Quantity): Angular size of the image desired
            filt (str): One of 'g','r','i','z','y' (default: 'i')
            timeout (int): Number of seconds to timout the query (default: 120 s)
        Returns:
            hdu: fits header data unit for the downloaded image
        """
        assert len(filt)==1 and filt in "grizy", "Filter name must be one of 'g','r','i','z','y'"
        url = _get_url(self.coord,imsize=imsize,filt=filt,imgformat='fits')[0]
        imagedat = fits.open(astroutils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))[0]
        return imagedat



def _get_url(coord,imsize=30*u.arcsec,filt="i",output_size=None,imgformat="fits",color=False):
    """
    Returns the url corresponding to the requested image cutout
    Args:
        coord (astropy SkyCoord): Center of the search area.
        imsize (astropy Angle): Length and breadth of the search area.
        filt (str): 'g','r','i','z','y'
        output_size (int): display image size (length) in pixels
        imgformat (str): "fits","png" or "jpg"
    """
    assert imgformat in ['jpg','png','fits'], "Image file can be only in the formats 'jpg', 'png' and 'fits'."
    if color:
        assert len(filt)==3,"Three filters are necessary for a color image"
        assert imgformat in ['jpg','png'], "Color image not available in fits format"
    
    pixsize = int(imsize.to(u.arcsec).value/0.25) #0.25 arcsec per pixel
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    filetaburl = ("{:s}?ra={:f}&dec={:f}&size={:d}&format=fits"
           "&filters={:s}").format(service,coord.ra.value,
                                        coord.dec.value, pixsize,filt)
    file_extensions = Table.read(filetaburl, format='ascii')['filename']

    url = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra={:f}&dec={:f}&size={:d}&format={:s}".format(coord.ra.value,coord.dec.value,
                                                                                                        pixsize,imgformat)
    if output_size:
        url += "&output_size={}".format(output_size)
    if color:
        cols = ['red','green','blue']
        for col,extension in zip(cols,file_extensions):
            url += "&{}={}".format(col,extension)
    else:
        urlbase = url + "&red="
        url = []
        for extensions in file_extensions:
            url.append(urlbase+extensions)
    return url
 
def _check_columns(columns,table,release):
    """
    Checks if the requested columns are present in the
    table from which data is to be pulled. Raises an error
    if those columns aren't found.
    Args:
        columns (list of str): column names to retrieve
        table (str): "mean","stack" or "detection"
        release (str): "dr1" or "dr2"
    """
    dcols = {}
    for col in _ps1metadata(table,release)['name']:
        dcols[col.lower()] = 1
    badcols = []
    for col in columns:
        if col.lower().strip() not in dcols:
            badcols.append(col)
    if badcols:
        raise ValueError('Some columns not found in table: {}'.format(', '.join(badcols)))

def _check_legal(table,release):
    """
    Checks if this combination of table and release is acceptable
    Raises a VelueError exception if there is problem.
    Taken from http://ps1images.stsci.edu/ps1_dr2_api.html
    Args:
        table (str): "mean","stack" or "detection"
        release (str): "dr1" or "dr2"
    """
    
    releaselist = ("dr1", "dr2")
    if release not in releaselist:
        raise ValueError("Bad value for release (must be one of {})".format(', '.join(releaselist)))
    if release=="dr1":
        tablelist = ("mean", "stack")
    else:
        tablelist = ("mean", "stack", "detection")
    if table not in tablelist:
        raise ValueError("Bad value for table (for {} must be one of {})".format(release, ", ".join(tablelist)))

def _ps1metadata(table="stack",release="dr2",
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"):
    """Return metadata for the specified catalog and table
    
    Args:
        table (string): mean, stack, or detection
        release (string): dr1 or dr2
        baseurl: base URL for the request
    
    Returns an astropy table with columns name, type, description
    """
    
    _check_legal(table,release)
    url = "{baseurl}/{release}/{table}/metadata".format(**locals())
    r = requests.get(url)
    r.raise_for_status()
    v = r.json()
    # convert to astropy table
    tab = Table(rows=[(x['name'],x['type'],x['description']) for x in v],
               names=('name','type','description'))
    return tab
