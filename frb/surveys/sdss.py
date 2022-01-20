""" Methods related to SDSS/BOSS queries """

import numpy as np
from IPython import embed

from astropy import units
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.table import Table

try:
    from astroquery.sdss import SDSS
except ImportError:
    print("Warning: You need to install astroquery to use the survey tools...")

from frb.surveys import surveycoord
from frb.surveys import catalog_utils
from frb.surveys import images

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

class SDSS_Survey(surveycoord.SurveyCoord):
    """
    Class to handle queries on the SDSS database

    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate

    """
    def __init__(self, coord, radius, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        #
        self.survey = 'SDSS'

    def get_catalog(self, photoobj_fields=None, timeout=120, print_query=False):
        """
        Query SDSS for all objects within a given
        radius of the input coordinates.

        Merges photometry with photo-z

        TODO -- Expand to include spectroscopy
        TODO -- Consider grabbing all of the photometry fields

        Args:
            coord: astropy.coordiantes.SkyCoord
            radius: Angle, optional
              Search radius
            photoobj_fields: list
              Fields for querying
            timeout: float, optional
              Default value - 120 s.
            print_query: bool, optional
              Print the SQL query for the photo-z values

        Returns:
            catalog: astropy.table.Table
              Contains all measurements retieved
              *WARNING* :: The SDSS photometry table frequently has multiple entries for a given
              source, with unique objid values

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
            self.catalog = Table()
            self.catalog.meta['radius'] = self.radius
            self.catalog.meta['survey'] = self.survey
            # Validate
            self.validate_catalog()
            return

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
            matches = catalog_utils.match_ids(photz_cat['objid'], photom_catalog['objid'], require_in_match=False)
        else:
            matches = -1 * np.ones(len(photom_catalog), dtype=int)
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
        trim_catalog = catalog_utils.clean_cat(trim_catalog, photom['SDSS'])

        # Spectral info
        spec_fields = ['ra', 'dec', 'z', 'run2d', 'plate', 'fiberID', 'mjd', 'instrument']
        spec_catalog = SDSS.query_region(self.coord,spectro=True, radius=self.radius,
                                         timeout=timeout, specobj_fields=spec_fields) # Duplicates may exist
        if spec_catalog is not None:
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

    def get_cutout(self, imsize, scale=0.396127):
        """
        Grab a cutout from SDSS

        Args:
            imsize (Quantity):  Size of image desired

        Returns:
            PNG image, None: self.cutout and a None to match the image header (not provided by SDSS)


        """
        # URL
        sdss_url = get_url(self.coord, imsize=imsize.to('arcsec').value,
                           scale=scale)
        # Image
        self.cutout = images.grab_from_url(sdss_url)
        self.cutout_size = imsize

        # Return
        return self.cutout, None


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


def trim_down_catalog(catalog, keep_photoz=False, cut_within=1.5*units.arcsec):
    """
    Cut down a catalog to keep only 1 source within cut_within

    Args:
        catalog (astropy.table.Table):  Input source catalog
        keep_photoz (bool, optional):
        cut_within (Angle or Quantity):  Cut radius

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


