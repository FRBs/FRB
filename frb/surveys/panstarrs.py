"""
Slurp data from Pan-STARRS catalog
"""

import numpy as np

from astropy import units as u,utils as astroutils
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table
from ..galaxies.defs import PanSTARRS_bands

from .images import grab_from_url

import warnings

try:
    from astroquery.vizier import Vizier
except ImportError:
    warnings.warn("Warning: You need to install astroquery to use the survey tools...")

from frb.surveys import surveycoord,catalog_utils,images

#TODO: It's potentially viable to use the same code for other
#catalogs in the VizieR database. Maybe a generalization wouldn't
#be too bad in the future.

# Define the data model for Pan-STARRS data
photom = {}
photom['Pan-STARRS'] = {}
for band in PanSTARRS_bands:
    photom["Pan-STARRS"]["Pan-STARRS"+'_{:s}'.format(band)] = '{:s}mag'.format(band.lower())
    photom["Pan-STARRS"]["Pan-STARRS"+'_{:s}_err'.format(band)] = 'e_{:s}mag'.format(band.lower())
    photom["Pan-STARRS"]["Pan-STARRS_ID"] = 'objID'
photom["Pan-STARRS"]['ra'] = 'RAJ2000'
photom["Pan-STARRS"]['dec'] = 'DEJ2000'
photom["Pan-STARRS"]["Pan-STARRS_field"] = 'field'

class Pan_STARRS_Survey(surveycoord.SurveyCoord):
    """
    A class to access all the catalogs hosted on the
    Vizier database. Inherits from SurveyCoord. This
    is a super class not meant for use by itself and
    instead meant to instantiate specific children
    classes like PAN-STARRS_Survey
    """
    def __init__(self,coord,radius,**kwargs):
        surveycoord.SurveyCoord.__init__(self,coord,radius,**kwargs)

        self.Survey = "Pan_STARRS"
    
    def get_catalog(self,query_fields=None,timeout=120,print_query=False):
        """
        Query a catalog in the VizieR database for
        photometry.

        Args:
            query_fields: list, optional
                A list of query fields to
                get.
            timeout: float, optional
                Query timeout in sec. Exits with an error
                if the query time exceeds this value.
                Default value of 120 s.
            print_query: bool, optional
                If true, prints the SQL query used
                on screen.
        
        Returns:
            catalog: astropy.table.Table
                Contains all query results
        """
        if query_fields is None:
            query_fields = ['objID','RAJ2000','DEJ2000','f_objID','Qual']
            query_fields +=['{:s}mag'.format(band) for band in PanSTARRS_bands]
            query_fields +=['e_{:s}mag'.format(band) for band in PanSTARRS_bands]
        vclient = Vizier(columns=query_fields,timeout=timeout,row_limit=-1)
        tablelist = vclient.query_region(self.coord,radius=self.radius,catalog="Pan-STARRS")
        if len(tablelist)==0:
            self.catalog = Table()
            self.catalog.meta['radius'] = self.radius
            self.catalog.meta['survey'] = self.survey
            # Validate
            self.validate_catalog()
            return self.catalog.copy()
        photom_catalog = tablelist[0]
        pdict = photom['Pan-STARRS']
        photom_catalog = catalog_utils.clean_cat(photom_catalog,pdict)
        #
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
                                to the oiginal cutout size.
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
        url = get_url(self.coord,imsize=imsize,filt=newfilt,output_size=output_size,color=True,imgformat='png')
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
        url = get_url(self.coord,imsize=imsize,filt=filt,imgformat='fits')[0]
        imagedat = fits.open(astroutils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))[0]
        return imagedat



def get_url(coord,imsize=30*u.arcsec,filt="i",output_size=None,imgformat="fits",color=False):

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

