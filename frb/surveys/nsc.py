"""NOIRLab source catalog"""

import re
import warnings

import numpy as np

from astropy import units, utils, wcs
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.table import Table

from frb.surveys import dlsurvey, defs
from frb.surveys import catalog_utils

# Dependencies
try:
    from pyvo.dal import sia
except ImportError:
    print("Warning:  You need to install pyvo to retrieve DES images")
    _svc = None
    DALFormatError = Exception
else:
    from pyvo.dal import DALFormatError
    _svc = sia.SIAService(defs.NOIR_DEF_ACCESS_URL+'nsa')

# The image service holds the individual DECam CCD images (there are no stacks).
# The resampled images are north-up, which the cutout service needs to
# return a cutout of the requested size.
_IMG_PROCTYPE = 'Resampled'
# Half of the short side of a DECam CCD (~9 x 18 arcmin). A cutout is sure
# to fit on a CCD whose center is closer to the target than this, less the
# half-diagonal of the cutout.
_CCD_HALF_WIDTH = 4.4*units.arcmin
# The cutout service rounds the size of a cutout by a few pixels. A cutout
# covering at least this fraction of the requested size is complete.
_COMPLETE_COVERAGE = 0.98

# Define the data model for DES data
photom = {}
photom['NSC'] = {}
photom['NSC']['NSC_ID'] = 'id'
photom['NSC']['ra'] = 'ra'
photom['NSC']['dec'] = 'dec'
photom['NSC']['class_star'] = 'class_star'
NSC_bands = ['u','g', 'r', 'i', 'z', 'Y', 'VR']
for band in NSC_bands:
    photom['NSC']['NSC_{:s}'.format(band)] = '{:s}mag'.format(band.lower())
    # Uncertainty of the mean magnitude (not the scatter between detections, <band>rms)
    photom['NSC']['NSC_{:s}_err'.format(band)] = '{:s}err'.format(band.lower())

# Columns of the catalog that are not floats; used for the schema of empty catalogs
schema_dtypes = {}
schema_dtypes['NSC'] = {'NSC_ID': str}

class NSC_Survey(dlsurvey.DL_Survey):
    """
    Class to handle queries on the NSC survey

    Child of DL_Survey which uses datalab to access NOAO

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.dlsurvey.DL_Survey`
            (e.g. ``verbose``)

    """

    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        dlsurvey.DL_Survey.__init__(self, coord, radius, **kwargs)
        self.survey = 'NSC'
        self.bands = NSC_bands
        self.svc = _svc
        self.qc_profile = "default"
        self.database = "nsc_dr2.object"
        self.default_query_fields = list(photom['NSC'].values())

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
        main_cat = super(NSC_Survey, self).get_catalog(query=query,
                                                       query_fields=query_fields,
                                                       print_query=print_query,**kwargs)
        main_cat = catalog_utils.clean_cat(main_cat, photom['NSC'], mask_photometry=True)
        # Empty catalogs get the standard columns (no-op otherwise)
        main_cat = catalog_utils.ensure_empty_schema(main_cat, list(photom['NSC'].keys()),
                                                     dtypes=schema_dtypes['NSC'])
        
        # Finish
        self.catalog = main_cat
        self.validate_catalog()
        return self.catalog

    @staticmethod
    def _cutout_url(access_url: str, imsize: units.Quantity) -> str:
        """
        Set the size of the cutout requested by an image access URL.

        Args:
            access_url (str): Access URL of an image, as returned by the
                image service. It asks for a cutout of the size of the search.
            imsize (astropy.units.Quantity): Angular size of the cutout wanted.

        Returns:
            str: URL of the cutout of size ``imsize``.

        """
        size = imsize.to(units.deg).value
        url, nsub = re.subn(r'SIZE=[^&]*', f'SIZE={size:.6f},{size:.6f}', access_url)
        if nsub == 0:
            warnings.warn("Could not set the cutout size in the NSC image URL; "
                          "the image may not have the requested size.", RuntimeWarning)
        return url

    @staticmethod
    def _zero_points(column) -> np.ndarray:
        """
        Read the zero points of the images, which measure their depth.

        Args:
            column (astropy.table.Column): The 'magzero' column of the image
                table. It may hold strings and masked entries.

        Returns:
            numpy.ndarray: The zero points as floats; -inf where there is none.

        """
        values = []
        for item in column:
            try:
                value = float(item)
            except (TypeError, ValueError):
                value = np.nan
            values.append(value if np.isfinite(value) else -np.inf)
        return np.array(values, dtype=float)

    def _download_cutout(self, url: str, timeout: int | float = 120) -> fits.PrimaryHDU:
        """
        Download a cutout image.

        Args:
            url (str): URL of the cutout.
            timeout (int or float, optional): Time to wait in seconds before timing out.

        Returns:
            astropy.io.fits.PrimaryHDU: The cutout.

        """
        filename = utils.data.download_file(url, cache=True, show_progress=False, timeout=timeout)
        with fits.open(filename) as hdul:
            return fits.PrimaryHDU(data=np.array(hdul[0].data), header=hdul[0].header.copy())

    @staticmethod
    def _cutout_coverage(hdu: fits.PrimaryHDU, imsize: units.Quantity) -> float:
        """
        Fraction of the requested size that a cutout covers along its
        worst axis. It is below 1 when the target is near the edge of the CCD.

        Args:
            hdu (astropy.io.fits.PrimaryHDU): The cutout.
            imsize (astropy.units.Quantity): Angular size that was requested.

        Returns:
            float: Covered fraction, between 0 and 1.

        """
        if hdu.data is None or hdu.data.ndim != 2:
            return 0.
        pixscale = np.mean(wcs.utils.proj_plane_pixel_scales(wcs.WCS(hdu.header)))*units.deg
        npix = (imsize/pixscale).decompose().value
        # Allow for the rounding of the number of pixels
        return float(min(1., (min(hdu.data.shape) + 1)/npix))

    def get_image(self, imsize: units.Quantity, band: str | None = None,
                  timeout: int | float = 120, verbose: bool = False,
                  search_size: units.Quantity = 0.2*units.deg,
                  max_tries: int = 5) -> fits.PrimaryHDU | None:
        """
        Get a FITS cutout of any size from the DECam images behind the NSC.

        The image service holds the individual (resampled) CCD images. They
        are searched for in a region of at least ``search_size``, because
        small search regions miss images that cover the target, and the cutout
        of size ``imsize`` is then requested from the deepest image that is
        sure to contain it. If the cutout comes back truncated, because the
        target is near the edge of the CCD, the next image is tried.

        Args:
            imsize (astropy.units.Quantity): Angular size of the image.
                Cutouts larger than a DECam CCD (~9 arcmin) cannot be complete.
            band (str, optional): Band of the image, one of ``self.bands``
                (case-insensitive). If None, 'r' is used.
            timeout (int or float, optional): Time to wait in seconds before timing out
            verbose (bool, optional): Print status
            search_size (astropy.units.Quantity, optional): Minimum size of
                the region searched for images. Increase it if no image is found.
            max_tries (int, optional): Maximum number of images to download
                when looking for one that fully covers the cutout.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The cutout; None if the
            image service fails or has no image in this band. If no image
            covers the full cutout, the one that covers most of it is returned
            with a warning.

        Raises:
            RuntimeError: If the image service (``self.svc``) is not set,
                e.g. because pyvo is not installed.
            TypeError: If ``band`` is not one of ``self.bands``.

        """
        if self.svc is None:
            raise RuntimeError("svc attribute cannot be None. Have you installed pyvo?")
        if band is None:
            band = 'r'
        allowed = [item.lower() for item in self.bands]
        if band.lower() not in allowed:
            raise TypeError("Allowed filters (case-insensitive) for {:s} photometric bands are {}".format(
                self.survey, self.bands))

        # Find the CCD images around the target
        search = max(imsize.to(units.deg), search_size.to(units.deg))
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                imgTable = self.svc.search(self.coord, search, verbosity=2).to_table()
        except DALFormatError:
            warnings.warn(f"Image cannot be retrieved. Invalid base URL?: {self.svc._baseurl}.",
                          RuntimeWarning)
            return None
        if verbose:
            print("The full image list contains", len(imgTable), "entries")
        if len(imgTable) > 0:
            selection = np.char.lower(np.asarray(imgTable['obs_bandpass']).astype(str)) == band.lower()
            selection &= np.asarray(imgTable['proctype']).astype(str) == _IMG_PROCTYPE
            selection &= np.asarray(imgTable['prodtype']).astype(str) == 'image'
            imgTable = imgTable[selection]
        if len(imgTable) == 0:
            print('No image available')
            return None

        # Order the images: first those sure to contain the cutout, deepest
        # first, and then the others, starting with the nearest.
        centers = SkyCoord(imgTable['s_ra'], imgTable['s_dec'], unit='deg')
        offset = self.coord.separation(centers).to(units.arcmin)
        depth = self._zero_points(imgTable['magzero'])
        safe = offset + imsize/np.sqrt(2) <= _CCD_HALF_WIDTH
        order = np.lexsort((np.where(safe, -depth, offset.value), ~safe))

        # Download, until a complete cutout is found
        best_hdu, best_coverage = None, -1.
        for idx in order[:max_tries]:
            url = self._cutout_url(str(imgTable['access_url'][idx]), imsize)
            try:
                hdu = self._download_cutout(url, timeout=timeout)
            except Exception as exc:
                if verbose:
                    print(f"Could not download {url}: {exc}")
                continue
            coverage = self._cutout_coverage(hdu, imsize)
            if verbose:
                print(f"Image {offset[idx]:.2f} from the target covers {100*coverage:.0f}% of the cutout")
            if coverage > best_coverage:
                best_hdu, best_coverage = hdu, coverage
            if coverage >= _COMPLETE_COVERAGE:
                break
        if best_hdu is None:
            print('No image available')
        elif best_coverage < _COMPLETE_COVERAGE:
            warnings.warn(f"The NSC image covers only {100*best_coverage:.0f}% of the requested "
                          f"{imsize} along one axis; the target is near the edge of the CCD "
                          "(or the image is larger than a CCD).", RuntimeWarning)
        return best_hdu
