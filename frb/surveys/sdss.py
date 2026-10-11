""" Methods related to SDSS/BOSS queries """

import warnings
import numpy as np

from astropy import units
from astropy.coordinates import Angle, SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.io.fits import PrimaryHDU
from astropy.table import Table

try:
    from astroquery.sdss import SDSS
except ImportError:
    print("Warning: You need to install astroquery to use the survey tools...")

from frb.surveys import surveycoord
from frb.surveys import catalog_utils
from frb.surveys import images
from frb.surveys.skyview import SkyView_Survey

# Define the data model for SDSS data
photom = {}
photom['SDSS'] = {}
SDSS_bands = ['u', 'g', 'r', 'i', 'z']
for band in SDSS_bands:
    photom['SDSS']['SDSS_{:s}'.format(band)] = 'modelMag_{:s}'.format(band.lower())
    photom['SDSS']['SDSS_{:s}_err'.format(band)] = 'modelMagErr_{:s}'.format(band.lower())
photom['SDSS']['SDSS_ID'] = 'objid'
photom['SDSS']['ra'] = 'ra'
photom['SDSS']['dec'] = 'dec'
photom['SDSS']['SDSS_field'] = 'field'

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['SDSS'] = {'SDSS_ID': 'uint64', 'SDSS_field': int}

class SDSS_Survey(SkyView_Survey):
    """
    Class to handle queries on the SDSS database

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.skyview.SkyView_Survey`
            (e.g. ``verbose``)

    """
    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        SkyView_Survey.__init__(self, coord, radius, 'sdss', **kwargs)
        #
        self.survey = 'SDSS'

    def get_image(self, imsize: units.Quantity, band: str = 'r') -> PrimaryHDU | None:
        """
        Retrieve a SkyView FITS image for SDSS.

        Args:
            imsize (astropy.units.Quantity): Angular size of desired image.
            band (str, optional): One of ``u``, ``g``, ``r``, ``i``, ``z``.

        Returns:
            astropy.io.fits.PrimaryHDU or None: FITS image product.
        """
        return SkyView_Survey.get_image(self, imsize=imsize, band=band)

    def get_catalog(self, photoobj_fields: list[str] | None = None,
                    timeout: int | float = 120, print_query: bool = False) -> Table:
        """
        Query SDSS for all objects within ``self.radius``
        of ``self.coord``.

        Merges photometry with photo-z

        TODO -- Expand to include spectroscopy
        TODO -- Consider grabbing all of the photometry fields

        Args:
            photoobj_fields (list of str, optional): Fields for querying.
                If None, the position, IDs, type, model magnitudes (and
                errors) and extinctions are queried.
            timeout (int or float, optional): Time to wait in seconds
                before timing out the queries. Default value - 120 s.
            print_query (bool, optional): Print the SQL query for the photo-z values

        Returns:
            astropy.table.Table: Contains all measurements retrieved.
            Empty if SDSS has no sources in the cone.
            *WARNING* :: The SDSS photometry table frequently has multiple entries for a given
            source, with unique objid values

        Raises:
            RuntimeError: If the SDSS photometry query fails.
            ValueError: If the SDSS spectroscopic sources cannot be matched
                to the photometric ones within 1.5 arcsec.

        """
        if photoobj_fields is None:
            photoobj_fs = ['ra', 'dec', 'objid', 'run', 'rerun', 'camcol', 'field','type']
            mags = ['modelMag_'+band for band in SDSS_bands]
            magsErr = ['modelMagErr_'+band for band in SDSS_bands]
            extinct = ["extinction_"+band for band in SDSS_bands]
            photoobj_fields = photoobj_fs+mags+magsErr+extinct

        # Call
        photom_catalog = SDSS.query_region(self.coord, radius=self.radius, timeout=timeout,
                                           photoobj_fields=photoobj_fields)
        if photom_catalog is None:
            self.catalog = catalog_utils.ensure_empty_schema(
                Table(), list(photom['SDSS'].keys()), dtypes=schema_dtypes['SDSS']
            )
            self.catalog.meta['radius'] = self.radius
            self.catalog.meta['survey'] = self.survey
            # Validate
            self.validate_catalog()
            return self.catalog.copy()
        elif '<html>' in photom_catalog.colnames[0]:
            raise RuntimeError("SDSS photometry query appears to have failed. Error message: {}".format(photom_catalog.colnames[0]+photom_catalog[0][0]))

        # Now query for photo-z
        query = "SELECT GN.distance, "
        query += "p.objid, "

        query += "pz.z as redshift, pz.zErr as redshift_error\n"
        query += "FROM PhotoObj as p\n"
        query += "JOIN dbo.fGetNearbyObjEq({:f},{:f},{:f}) AS GN\nON GN.objID=p.objID\n".format(
            self.coord.ra.value,self.coord.dec.value,self.radius.to('arcmin').value)
        query += "JOIN Photoz AS pz ON pz.objID=p.objID\n"
        query += "ORDER BY distance"

        if print_query:
            print(query)

        # SQL command
        photz_cat = SDSS.query_sql(query,timeout=timeout)

        # Match em up
        if photz_cat is not None:
            # Was there an error in the photo-z query?
            if (len(photz_cat.colnames) == 1) &('<html>' in photz_cat.colnames[0]):
                warnings.warn("Photo-z query appears to have failed; no photo-z will be included",category=RuntimeWarning)
                photz_cat = None
                matches = np.full(len(photom_catalog), fill_value=-1, dtype=int)
            else:
                matches = catalog_utils.match_ids(photz_cat['objid'], photom_catalog['objid'], require_in_match=False)
        else:
            matches = np.full(len(photom_catalog), fill_value=-1, dtype=int)
        gdz = matches > 0
        # Init
        photom_catalog['photo_z'] = -9999.
        photom_catalog['photo_zerr'] = -9999.
        # Fill
        if np.any(gdz):
            photom_catalog['photo_z'][matches[gdz]] = photz_cat['redshift'][np.where(gdz)]
            photom_catalog['photo_zerr'][matches[gdz]] = photz_cat['redshift_error'][np.where(gdz)]

        # Trim down catalog
        trim_catalog = trim_down_catalog(photom_catalog, keep_photoz=True)

        # Clean up
        trim_catalog = catalog_utils.clean_cat(trim_catalog, photom['SDSS'], mask_photometry=True)

        # Spectral info
        spec_fields = ['ra', 'dec', 'z', 'run2d', 'plate', 'fiberID', 'mjd', 'instrument']
        spec_catalog = SDSS.query_region(self.coord,spectro=True, radius=self.radius,
                                         timeout=timeout, specobj_fields=spec_fields) # Duplicates may exist
        
        # Make sure the returned spec_catalog isn't bad
        if spec_catalog == None:
            bad_spec = True
        elif len(spec_catalog) == 0:
            bad_spec = True
        elif len(spec_catalog.colnames) == 1 and '<html>' in spec_catalog.colnames[0]:
            bad_spec = True
        else:
            bad_spec = False

        if not bad_spec:
            trim_spec_catalog = trim_down_catalog(spec_catalog)
            # Match
            spec_coords = SkyCoord(ra=trim_spec_catalog['ra'], dec=trim_spec_catalog['dec'], unit='deg')
            phot_coords = SkyCoord(ra=trim_catalog['ra'], dec=trim_catalog['dec'], unit='deg')
            idx, d2d, d3d = match_coordinates_sky(spec_coords, phot_coords, nthneighbor=1)
            # Check
            if np.max(d2d).to('arcsec').value > 1.5:
                raise ValueError("Bad match in SDSS")
            # Fill me
            zs = -1 * np.ones_like(trim_catalog['ra'].data)
            zs[idx] = trim_spec_catalog['z']
            trim_catalog['z_spec'] = zs
        else:
            trim_catalog['z_spec'] = -1.

        # Sort by offset
        catalog = trim_catalog.copy()
        self.catalog = catalog_utils.sort_by_separation(catalog, self.coord, radec=('ra','dec'), add_sep=True)

        # Meta
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey

        # Validate
        self.validate_catalog()

        # Return
        return self.catalog.copy()

    def get_cutout(self, imsize: units.Quantity, scale: float = 0.396127) -> tuple:
        """
        Grab a cutout from SDSS

        Args:
            imsize (astropy.units.Quantity):  Size of image desired
            scale (float, optional): Pixel scale in arcsec/pixel

        Returns:
            tuple: ``(cutout, None)``: the image (numpy.ndarray, from the JPEG
            served by SDSS), stored in ``self.cutout``, and a None to match the
            image header (not provided by SDSS)

        """
        # URL
        sdss_url = get_url(self.coord, imsize=imsize.to('arcsec').value,
                           scale=scale)
        # Image
        self.cutout = images.grab_from_url(sdss_url)
        self.cutout_size = imsize

        # Return
        return self.cutout, None


def get_url(coord: SkyCoord, imsize: float = 30., scale: float = 0.396127,
            grid: bool = False, label: bool = False, invert: bool = False) -> str:
    """
    Generate the SDSS URL for an image retrieval

    Args:
        coord (astropy.coordinates.SkyCoord): Center of image
        imsize (float, optional): Image size (rectangular) in arcsec and without units
        scale (float, optional): Pixel scale in arcsec/pixel
        grid (bool, optional): Overlay a grid on the image
        label (bool, optional): Label the image
        invert (bool, optional): Invert the colors of the image

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


def trim_down_catalog(catalog: Table, keep_photoz: bool = False,
                      cut_within: Angle | units.Quantity = 1.5*units.arcsec) -> Table:
    """
    Cut down a catalog to keep only 1 source within cut_within

    Args:
        catalog (astropy.table.Table):  Input source catalog
        keep_photoz (bool, optional): Prefer sources with a photo-z (which
            must be in the 'photo_z' column) when choosing which to keep
        cut_within (astropy.coordinates.Angle or astropy.units.Quantity, optional):  Cut radius

    Returns:
        astropy.table.Table:  Catalog trimmed down

    """
    if len(catalog) == 1:
        return catalog

    # All good
    keep = np.ones_like(catalog, dtype=bool)

    coords = SkyCoord(ra=catalog['ra'], dec=catalog['dec'], unit='deg')
    # Search on closest next neighbor
    idx, d2d, d3d = match_coordinates_sky(coords, coords, nthneighbor=2)
    too_close = np.where(d2d < cut_within)[0]
    for idx in too_close:
        # Already purged?
        if not keep[idx]:
            continue
        # Find the matches
        seps = coords[idx].separation(coords)
        orig_matches = np.where((seps < cut_within) & keep)[0]
        final_matches = orig_matches.copy()
        # Keep photo-z?
        if keep_photoz:
            good_pz = catalog['photo_z'][final_matches] > -9000.
            if np.any(good_pz):
                final_matches = final_matches[good_pz]
        # Take the first one -- Any reason to do otherwise?
        final_match = final_matches[0]
        # Zero out the rest
        zero_me = orig_matches != final_match
        keep[orig_matches[zero_me]] = False
    # Finish
    return catalog[keep]


