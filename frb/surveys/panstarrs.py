"""
Slurp data from Pan-STARRS catalog using the MAST API.
A lot of this code has been directly taken from
http://ps1images.stsci.edu/ps1_dr2_api.html

"""

import numpy as np

from astropy import units as u,utils as astroutils
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.table import Table, join
from ..galaxies.defs import PanSTARRS_bands
import importlib_resources

import warnings
import requests

from frb.surveys import surveycoord,catalog_utils,images

import os

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

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['Pan-STARRS'] = {'Pan-STARRS_ID': int}

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
    A class to access the Pan-STARRS catalogs hosted on the
    MAST database. Inherits from SurveyCoord.

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """
    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        surveycoord.SurveyCoord.__init__(self,coord,radius,**kwargs)

        self.Survey = "Pan_STARRS"
    
    def get_catalog(self, query_fields: list[str] | None = None,
                    release: str = "dr2", table: str = "stack",
                    print_query: bool = False, use_psf: bool = False,
                    photoz: bool = True) -> Table:
        """
        Query a catalog in the MAST Pan-STARRS database for
        photometry.

        Args:
            query_fields (list of str, optional): A list of query fields to
                get in addition to the default fields.
            release (str, optional): "dr1" or "dr2" (default: "dr2").
                Data release version.
            table (str, optional): "mean","stack" or "detection"
                (default: "stack"). The data table to
                search within.
            print_query (bool, optional): Print the URL of the query
            use_psf (bool, optional): If True, use PSFmag instead of KronMag
            photoz (bool, optional): If True, also download photometric redshifts
                using the Mast CasJobs API.

        Returns:
            astropy.table.Table: Contains all query results

        Raises:
            ValueError: If the ``release``/``table`` combination is not
                allowed, or a requested field is not in the table.
            IOError: If ``photoz`` is True and the MAST_CASJOBS_USER and
                MAST_CASJOBS_PWD environment variables are not set.

        """
        #assert self.radius <= 0.5*u.deg, "Cone serches have a maximum radius"
        #Validate table and release input
        _check_legal(table,release)
        url = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/{:s}/{:s}".format(release,table)
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
        data['format'] = 'csv'
        if print_query:
            print(url)
        ret = requests.get(url,params=data)
        ret.raise_for_status()
        if len(ret.text)==0:
            self.catalog = catalog_utils.ensure_empty_schema(
                Table(), list(photom['Pan-STARRS'].keys()),
                dtypes=schema_dtypes['Pan-STARRS']
            )
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
        
        photom_catalog = catalog_utils.clean_cat(photom_catalog, pdict, mask_photometry=True)

        #Remove bad positions because Pan-STARRS apparently decided
        #to flag some positions with large negative numbers. Why even keep
        #them?
        bad_ra = (photom_catalog['ra']<0)+(photom_catalog['ra']>360)
        bad_dec = (photom_catalog['dec']<-90)+(photom_catalog['dec']>90)
        bad_pos = bad_ra+bad_dec # bad_ra OR bad_dec
        photom_catalog = photom_catalog[~bad_pos]

        # Download photometric redshifts if requested.
        if photoz:
            try:
                import mastcasjobs as mcj

                # Query
                photoz_query = f"""SELECT m.objID, m.z_phot, m.z_photErr, m.class
                                  FROM fGetNearbyObjEq({self.coord.ra.value}, {self.coord.dec.value}, {self.radius.to('arcmin').value}) nb
                                  INNER JOIN catalogRecordRowStore m on m.objid=nb.objid"""
                # Execute
                user = os.getenv('MAST_CASJOBS_USER')
                pwd = os.getenv('MAST_CASJOBS_PWD')
                if user is None or pwd is None:
                    raise IOError("You need to set the MAST_CASJOBS_USER and MAST_CASJOBS_PWD environment variables. Create an account at https://mastweb.stsci.edu/mcasjobs/CreateAccount.aspx to get your credentials. Or set photoz=False in get_catalog.")
                job = mcj.MastCasJobs(context="HLSP_PS1_STRM", username=user, password=pwd)
                photoz_tab = job.quick(photoz_query, task_name="Photo-z cone search")
                photoz_tab.rename_column('objID', 'Pan-STARRS_ID')

                # Merge to the main tab
                photom_catalog = join(photom_catalog, photoz_tab, keys='Pan-STARRS_ID', join_type='left')

                # Now join the 
            except ImportError:
                warnings.warn("mastcasjobs not installed. Cannot download photometric redshifts.")
        
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

    def get_cutout(self, imsize: u.Quantity = 30*u.arcsec, band: str = "irg",
                   output_size: int | None = None, **kwargs) -> tuple:
        """
        Grab a color cutout (PNG) from Pan-STARRS

        Args:
            imsize (astropy.units.Quantity):  Angular size of image desired
            band (str, optional): A string with the three filters to be used
            output_size (int, optional): Output image size in pixels. Defaults
                to the original cutout size.
            **kwargs: Only the deprecated ``filt`` (use ``band``) is accepted.

        Returns:
            tuple: A 1-tuple holding the RGB image (numpy.ndarray).

        Raises:
            TypeError: If an unexpected keyword argument is given, or both
                ``band`` and the deprecated ``filt`` are specified.

        """
        if 'filt' in kwargs:
            warnings.warn(
                "'filt' is deprecated; use 'band' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if band != "irg":
                raise TypeError("Specify only one of 'band' or deprecated 'filt'.")
            band = kwargs.pop('filt')
        if kwargs:
            raise TypeError(f"Unexpected keyword arguments: {list(kwargs.keys())}")

        assert len(band) == 3, "Need three filters for a cutout."
        #Sort filters from red to blue
        band = band.lower() #Just in case the user is cheeky about the filter case.
        reffilt = "yzirg"
        idx = np.argsort([reffilt.find(f) for f in band])
        newband = ""
        for i in idx:
            newband += band[i]
        #Get image url
        url = _get_url(self.coord, imsize=imsize, band=newband, output_size=output_size,
                       color=True, imgformat='png')
        self.cutout = images.grab_from_url(url)
        self.cutout_size = imsize
        return  self.cutout.copy(), 
    
    def get_image(self, imsize: u.Quantity = 30*u.arcsec, band: str = "i",
                  timeout: int | float = 120, **kwargs) -> fits.PrimaryHDU:
        """
        Grab a fits image from Pan-STARRS in a
        specific band.

        Args:
            imsize (astropy.units.Quantity): Angular size of the image desired
            band (str, optional): One of 'g','r','i','z','y' (default: 'i')
            timeout (int or float, optional): Number of seconds to timout the query (default: 120 s)
            **kwargs: Only the deprecated ``filt`` (use ``band``) is accepted.

        Returns:
            astropy.io.fits.PrimaryHDU: fits header data unit for the downloaded image

        Raises:
            TypeError: If an unexpected keyword argument is given, or both
                ``band`` and the deprecated ``filt`` are specified.

        """
        if 'filt' in kwargs:
            warnings.warn(
                "'filt' is deprecated; use 'band' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if band != "i":
                raise TypeError("Specify only one of 'band' or deprecated 'filt'.")
            band = kwargs.pop('filt')
        if kwargs:
            raise TypeError(f"Unexpected keyword arguments: {list(kwargs.keys())}")

        assert len(band) == 1 and band in "grizy", "Filter name must be one of 'g','r','i','z','y'"
        url = _get_url(self.coord, imsize=imsize, band=band, imgformat='fits')[0]
        imagedat = fits.open(astroutils.data.download_file(url,cache=True,show_progress=False,timeout=timeout))[0]
        return imagedat



def _get_url(coord: SkyCoord, imsize: u.Quantity = 30*u.arcsec, band: str = "i",
             output_size: int | None = None, imgformat: str = "fits",
             color: bool = False, **kwargs) -> str | list[str]:
    """
    Returns the url corresponding to the requested image cutout

    Args:
        coord (astropy.coordinates.SkyCoord): Center of the search area.
        imsize (astropy.units.Quantity): Length and breadth of the search area.
        band (str, optional): 'g','r','i','z','y'; three of them for a color image
        output_size (int, optional): display image size (length) in pixels
        imgformat (str, optional): "fits","png" or "jpg"
        color (bool, optional): Request a color image (needs three filters in ``band``
            and a "png" or "jpg" ``imgformat``).
        **kwargs: Only the deprecated ``filt`` (use ``band``) is accepted.

    Returns:
        str or list of str: The URL of the color image if ``color`` is True,
        otherwise a list of URLs, one for each image of the requested band.

    """
    assert imgformat in ['jpg','png','fits'], "Image file can be only in the formats 'jpg', 'png' and 'fits'."
    if 'filt' in kwargs:
        warnings.warn(
            "'filt' is deprecated; use 'band' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if band != "i":
            raise TypeError("Specify only one of 'band' or deprecated 'filt'.")
        band = kwargs.pop('filt')
    if kwargs:
        raise TypeError(f"Unexpected keyword arguments: {list(kwargs.keys())}")

    if color:
        assert len(band) == 3,"Three filters are necessary for a color image"
        assert imgformat in ['jpg','png'], "Color image not available in fits format"
    
    pixsize = int(imsize.to(u.arcsec).value/0.25) #0.25 arcsec per pixel
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    filetaburl = ("{:s}?ra={:f}&dec={:f}&size={:d}&format=fits"
            "&filters={:s}").format(service,coord.ra.value,
                             coord.dec.value, pixsize,band)
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
 
def _check_columns(columns: list[str], table: str, release: str) -> None:
    """
    Checks if the requested columns are present in the
    table from which data is to be pulled. Raises an error
    if those columns aren't found.

    Args:
        columns (list of str): column names to retrieve
        table (str): "mean","stack" or "detection"
        release (str): "dr1" or "dr2"

    Raises:
        ValueError: If any of the columns is not in the table.

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

def _check_legal(table: str, release: str) -> None:
    """
    Checks if this combination of table and release is acceptable
    Raises a ValueError exception if there is problem.
    Taken from http://ps1images.stsci.edu/ps1_dr2_api.html

    Args:
        table (str): "mean","stack" or "detection"
        release (str): "dr1" or "dr2"

    Raises:
        ValueError: If ``release`` is not allowed, or ``table`` is not
            available for it.

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

def _ps1metadata(table: str = "stack", release: str = "dr2",
                 baseurl: str = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"
                 ) -> Table:
    """
    Return metadata for the specified catalog and table

    Args:
        table (str, optional): mean, stack, or detection
        release (str, optional): dr1 or dr2
        baseurl (str, optional): base URL for the request

    Returns:
        astropy.table.Table: Table with columns name, datatype, description

    Raises:
        IOError: If the server has no metadata and there is no local copy.

    """
    
    _check_legal(table,release)
    url = f"{baseurl}/{release}/{table}/metadata"
    r = requests.get(url)
    r.raise_for_status()
    v = r.json()
    # convert to astropy table
    local_metadata_path = importlib_resources.files('frb').joinpath('data','Public', 'Pan-STARRS','ps1_{}_{}_metadata.csv'.format(release,table))
    try:
        tab = Table(rows=[(x['name'],x['datatype'],x['description']) for x in v],
                    names=('name','datatype','description'))
        # Cache locally
        # Create directory if it doesn't exist
        if not os.path.isfile(local_metadata_path):
            os.makedirs(os.path.dirname(local_metadata_path), exist_ok=True)
        tab.write(local_metadata_path, overwrite=True)

    # The following catches the case when there is a server issue
    # and we have a local copy of the metadata.
    except KeyError:
        if not os.path.isfile(local_metadata_path):
            raise IOError("PS1 metadata not available from server and no local copy found.")
        else:
            tab = Table.read(local_metadata_path)

    return tab

