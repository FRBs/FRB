"""SkyView-backed survey image retrieval helpers.

Native pixel scales for SkyView image products
-----------------------------------------------
By default ``_skyview_fetch`` computes the ``pixels`` output dimension from
``imsize`` and the survey's native SkyView pixel scale below, so the returned
image is at full resolution.  Pass an explicit integer ``pixels`` value to
request a coarser (downsampled) grid.

+---------------------------+------------------+----------------------------------------------------+
| Survey                    | Pixel scale      | Source                                             |
+===========================+==================+====================================================+
| VLA FIRST (1.4 GHz)       |  1.8 arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
| NVSS                      | 15   arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
| WENSS                     | 21   arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
| TGSS ADR1                 |  6.2 arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
| GLEAM 72-103 MHz          | 56   arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
| GLEAM 103-134 MHz         | 44   arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
| GLEAM 139-170 MHz         | 34   arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
| GLEAM 170-231 MHz         | 28   arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
| SDSS u/g/r/i/z            |  0.4 arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
|   (camera native)         |  0.396 arcsec/pix| https://www.sdss4.org/instruments/camera/          |
| GALEX NUV / FUV           |  1.5 arcsec/pix  | https://skyview.gsfc.nasa.gov/current/cgi/survey.pl|
|   (mission native)        | ~1.5 arcsec/pix  | https://www.galex.caltech.edu/researcher/techdoc-ch5.html|
| 2MASS J/H/K (Atlas image) |  1.0 arcsec/pix  | https://irsa.ipac.caltech.edu/data/2MASS/docs/releases/allsky/doc/sec2_4.html|
|   (detector sampling)     | ~2.0 arcsec/pix  | https://irsa.ipac.caltech.edu/data/2MASS/docs/releases/allsky/doc/sec3_1b.html|
+---------------------------+------------------+----------------------------------------------------+
"""

import warnings

import numpy as np
from astropy import units as u
from astropy import wcs
from astropy.coordinates import Angle, SkyCoord
from astropy.io.fits import PrimaryHDU

try:
    from astroquery.skyview import SkyView
except ImportError:
    print("Warning:  You need astroquery installed to use SkyView survey tools")

from frb.surveys import surveycoord


class SkyView_Survey(surveycoord.SurveyCoord):
    """
    Class to handle queries to the SkyView service of `astroquery`.

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around.
        radius (astropy.coordinates.Angle): Search radius around the coordinate.
        mission (str): Mission served by SkyView for image searches.
            One of 'first', 'nvss', 'wenss', 'gleam', 'tgss', 'sdss',
            'galex' or '2mass' (case-insensitive).
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """

    SDSS_SURVEYS = {'u': 'SDSSu', 'g': 'SDSSg', 'r': 'SDSSr', 'i': 'SDSSi', 'z': 'SDSSz'}
    GALEX_SURVEYS = {'NUV': 'GALEX Near UV', 'FUV': 'GALEX Far UV'}
    TWOMASS_SURVEYS = {'J': '2MASS-J', 'H': '2MASS-H', 'K': '2MASS-K'}

    # Native SkyView pixel scales in arcsec/pixel; used to compute the ``pixels``
    # output dimension so images are returned at full (native) resolution by default.
    SKYVIEW_PIXEL_SCALES = {
        'VLA FIRST (1.4 GHz)': 1.8,
        'NVSS':                15.0,
        'WENSS':               21.0,
        'TGSS ADR1':            6.2,
        'GLEAM 72-103 MHz':    56.0,
        'GLEAM 103-134 MHz':   44.0,
        'GLEAM 139-170 MHz':   34.0,
        'GLEAM 170-231 MHz':   28.0,
        'SDSSu':                0.396,
        'SDSSg':                0.396,
        'SDSSr':                0.396,
        'SDSSi':                0.396,
        'SDSSz':                0.396,
        'GALEX Near UV':        1.5,
        'GALEX Far UV':         1.5,
        '2MASS-J':              1.0,
        '2MASS-H':              1.0,
        '2MASS-K':              1.0,
    }

    def __init__(self, coord: SkyCoord, radius: Angle, mission: str, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        self.survey = None
        self.mission = mission
        self.skyview = SkyView()

    @staticmethod
    def _coerce_imsize(imsize: u.Quantity | None = None,
                       radius: u.Quantity | None = None) -> u.Quantity:
        """
        Resolve the image size from the ``imsize`` and deprecated ``radius`` arguments.

        Args:
            imsize (astropy.units.Quantity, optional): Angular size (full side
                length) of the image.
            radius (astropy.units.Quantity, optional): Deprecated. Radius of
                the image; used (as ``2*radius``) only if ``imsize`` is None.

        Returns:
            astropy.units.Quantity: Angular size of the image.

        Raises:
            TypeError: If neither ``imsize`` nor ``radius`` is given.

        """
        if imsize is None and radius is None:
            raise TypeError("get_image() requires imsize")
        if radius is not None:
            warnings.warn(
                "radius is deprecated for SkyView-backed image retrieval; use imsize instead.",
                DeprecationWarning,
                stacklevel=3,
            )
            if imsize is None:
                imsize = 2 * radius
        return imsize

    def _skyview_fetch(self, skyview_name: str, imsize: Angle,
                       pixels: int | None = None) -> PrimaryHDU | None:
        """
        Fetch a FITS image from SkyView.

        Args:
            skyview_name (str): SkyView survey identifier.
            imsize (astropy.coordinates.Angle): Angular size of the image (full side length).
            pixels (int, optional): Output image side length in pixels.  If
                ``None`` (default), the side length is computed from ``imsize``
                and the survey's native SkyView pixel scale so the image is
                returned at full resolution.  Provide an integer smaller than
                the native value to request a coarser, downsampled grid.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None (with a
            warning) if SkyView returned no image.

        """
        radius = imsize / 2
        if pixels is None:
            pixel_scale = self.SKYVIEW_PIXEL_SCALES.get(skyview_name)
            if pixel_scale is not None:
                # At least one pixel, however small the requested image is
                pixels = max(1, int(round(imsize.to(u.arcsec).value / pixel_scale)))
        images = SkyView.get_images(
            position=self.coord, survey=skyview_name, radius=radius,
            pixels=str(pixels) if pixels is not None else None,
        )
        if not images or not images[0]:
            warnings.warn(f"SkyView returned no image for {skyview_name}.")
            return None
        return images[0][0]

    def get_image(self, imsize: u.Quantity | None = None, band: str | None = None,
                  radius: u.Quantity | None = None,
                  pixels: int | None = None) -> PrimaryHDU | None:
        """
        Retrieve a FITS image from SkyView for the mission of this survey.

        The image data and header are also stored in ``self.cutout`` and
        ``self.cutout_hdr``, and ``imsize`` in ``self.cutout_size``.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of
                the image. Required unless the deprecated ``radius`` is given.
            band (str, optional): Filter/band for missions that have several
                (SDSS: u/g/r/i/z, GALEX: NUV/FUV, 2MASS: J/H/K, or the
                frequency range of GLEAM, e.g. '170-231 MHz').
                Ignored by the other missions. If None, the mission default is used.
            radius (astropy.units.Quantity, optional): Deprecated; use ``imsize``.
            pixels (int, optional): Output image side length in pixels. If None,
                the native resolution is used. See :meth:`_skyview_fetch`.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None
            if SkyView returned none.

        Raises:
            NotImplementedError: If the mission is not supported.
            TypeError: If ``band`` is not valid for the mission, or no image
                size is given.

        """
        imsize = self._coerce_imsize(imsize=imsize, radius=radius)
        self.cutout_size = imsize

        mission = self.mission.lower()
        if mission == 'first':
            img_hdu = self.get_first(imsize, pixels=pixels)
        elif mission == 'nvss':
            img_hdu = self.get_nvss(imsize, pixels=pixels)
        elif mission == 'wenss':
            img_hdu = self.get_wenss(imsize, pixels=pixels)
        elif mission == 'gleam':
            if band is None:
                img_hdu = self.get_gleam(imsize, pixels=pixels)
            else:
                img_hdu = self.get_gleam(imsize, band=band, pixels=pixels)
        elif mission == 'tgss':
            img_hdu = self.get_tgss(imsize, pixels=pixels)
        elif mission == 'sdss':
            img_hdu = self.get_sdss(imsize, band=band, pixels=pixels)
        elif mission == 'galex':
            img_hdu = self.get_galex(imsize, band=band, pixels=pixels)
        elif mission == '2mass':
            img_hdu = self.get_twomass(imsize, band=band, pixels=pixels)
        else:
            raise NotImplementedError(f"SkyView mission '{self.mission}' is not supported")

        if img_hdu is None:
            self.cutout = None
            self.cutout_hdr = None
            return None

        self.cutout = img_hdu.data
        self.cutout_hdr = img_hdu.header

        mywcs = wcs.WCS(self.cutout_hdr)
        ypix, xpix = self.cutout.shape
        (ra0, dec0), (ra1, dec1), = mywcs.wcs_pix2world([[0, 0], [xpix, ypix]], 0)
        print("Got image spanning (RA, Dec) = ({0} - {1}, {2} - {3})".format(ra0, ra1, dec0, dec1))

        return img_hdu

    def get_cutout(self, imsize: u.Quantity | None = None, band: str | None = None,
                   radius: u.Quantity | None = None,
                   pixels: int | None = None) -> np.ndarray | None:
        """
        Deprecated alias of :meth:`get_image` that returns only the image array.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of
                the image.
            band (str, optional): Filter/band; see :meth:`get_image`.
            radius (astropy.units.Quantity, optional): Deprecated; use ``imsize``.
            pixels (int, optional): Output image side length in pixels.

        Returns:
            numpy.ndarray or None: The image data, or None if no image was
            retrieved. The header is stored in ``self.cutout_hdr``.

        """
        warnings.warn(
            "get_cutout() returns FITS products for this survey and is deprecated; use get_image() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        img_hdu = self.get_image(imsize=imsize, band=band, radius=radius, pixels=pixels)
        if img_hdu is None:
            self.cutout = None
            self.cutout_hdr = None
            return None
        self.cutout = img_hdu.data
        self.cutout_hdr = img_hdu.header
        return self.cutout

    def get_first(self, imsize: u.Quantity, pixels: int | None = None) -> PrimaryHDU | None:
        """
        Retrieve a VLA FIRST (1.4 GHz) FITS image from SkyView.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            pixels (int, optional): Output image side length in pixels.
                If None, the native resolution is used.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None if SkyView returned none.

        """
        return self._skyview_fetch('VLA FIRST (1.4 GHz)', imsize, pixels=pixels)

    def get_nvss(self, imsize: u.Quantity, pixels: int | None = None) -> PrimaryHDU | None:
        """
        Retrieve a NVSS FITS image from SkyView.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            pixels (int, optional): Output image side length in pixels.
                If None, the native resolution is used.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None if SkyView returned none.

        """
        return self._skyview_fetch('NVSS', imsize, pixels=pixels)

    def get_wenss(self, imsize: u.Quantity, pixels: int | None = None) -> PrimaryHDU | None:
        """
        Retrieve a WENSS FITS image from SkyView.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            pixels (int, optional): Output image side length in pixels.
                If None, the native resolution is used.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None if SkyView returned none.

        """
        return self._skyview_fetch('WENSS', imsize, pixels=pixels)

    def get_gleam(self, imsize: u.Quantity, band: str = '170-231 MHz',
                  pixels: int | None = None) -> PrimaryHDU | None:
        """
        Retrieve a GLEAM FITS image from SkyView.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            band (str, optional): Frequency range, one of '72-103 MHz',
                '103-134 MHz', '139-170 MHz' or '170-231 MHz'.
            pixels (int, optional): Output image side length in pixels.
                If None, the native resolution is used.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None if SkyView returned none.

        """
        return self._skyview_fetch(f'GLEAM {band}', imsize, pixels=pixels)

    def get_tgss(self, imsize: u.Quantity, pixels: int | None = None) -> PrimaryHDU | None:
        """
        Retrieve a TGSS ADR1 FITS image from SkyView.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            pixels (int, optional): Output image side length in pixels.
                If None, the native resolution is used.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None if SkyView returned none.

        """
        return self._skyview_fetch('TGSS ADR1', imsize, pixels=pixels)

    def get_sdss(self, imsize: u.Quantity, band: str | None = None,
             pixels: int | None = None) -> PrimaryHDU | None:
        """
        Retrieve a SDSS FITS image from SkyView.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            band (str, optional): Band, u, g, r, i or z (default r).
            pixels (int, optional): Output image side length in pixels.
                If None, the native resolution is used.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None if SkyView returned none.

        Raises:
            TypeError: If ``band`` is not allowed for SDSS.

        """
        band = 'r' if band is None else band.lower()
        if band not in self.SDSS_SURVEYS:
            raise TypeError(f"Allowed filters for SDSS are {list(self.SDSS_SURVEYS)}")
        return self._skyview_fetch(self.SDSS_SURVEYS[band], imsize, pixels=pixels)

    def get_galex(self, imsize: u.Quantity, band: str | None = None,
              pixels: int | None = None) -> PrimaryHDU | None:
        """
        Retrieve a GALEX FITS image from SkyView.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            band (str, optional): Band, NUV or FUV (default NUV).
            pixels (int, optional): Output image side length in pixels.
                If None, the native resolution is used.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None if SkyView returned none.

        Raises:
            TypeError: If ``band`` is not allowed for GALEX.

        """
        band = 'NUV' if band is None else band.upper()
        if band not in self.GALEX_SURVEYS:
            raise TypeError(f"Allowed filters for GALEX are {list(self.GALEX_SURVEYS)}")
        return self._skyview_fetch(self.GALEX_SURVEYS[band], imsize, pixels=pixels)

    def get_twomass(self, imsize: u.Quantity, band: str | None = None,
                pixels: int | None = None) -> PrimaryHDU | None:
        """
        Retrieve a 2MASS FITS image from SkyView.

        Args:
            imsize (astropy.units.Quantity): Angular size (full side length) of the image.
            band (str, optional): Band, J, H or K (default J).
            pixels (int, optional): Output image side length in pixels.
                If None, the native resolution is used.

        Returns:
            astropy.io.fits.PrimaryHDU or None: The image, or None if SkyView returned none.

        Raises:
            TypeError: If ``band`` is not allowed for 2MASS.

        """
        band = 'J' if band is None else band.upper()
        if band not in self.TWOMASS_SURVEYS:
            raise TypeError(f"Allowed filters for 2MASS are {list(self.TWOMASS_SURVEYS)}")
        return self._skyview_fetch(self.TWOMASS_SURVEYS[band], imsize, pixels=pixels)
