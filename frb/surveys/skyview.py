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

from astropy import units as u
from astropy import wcs

try:
    from astroquery.skyview import SkyView
except ImportError:
    print("Warning:  You need astroquery installed to use SkyView survey tools")

from frb.surveys import surveycoord


class SkyView_Survey(surveycoord.SurveyCoord):
    """
    Class to handle queries to the SkyView service of `astroquery`.

    Args:
        coord (SkyCoord): Coordinate for surveying around.
        radius (Angle): Search radius around the coordinate.
        mission (str): Mission served by SkyView for image searches.
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

    def __init__(self, coord, radius, mission, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        self.survey = None
        self.mission = mission
        self.skyview = SkyView()

    @staticmethod
    def _coerce_imsize(imsize=None, radius=None):
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

    def _skyview_fetch(self, skyview_name, imsize, pixels=None):
        """Fetch a FITS image from SkyView.

        Args:
            skyview_name (str): SkyView survey identifier.
            imsize (Angle): Angular size of the image (full side length).
            pixels (int, optional): Output image side length in pixels.  If
                ``None`` (default), the side length is computed from ``imsize``
                and the survey's native SkyView pixel scale so the image is
                returned at full resolution.  Provide an integer smaller than
                the native value to request a coarser, downsampled grid.
        """
        radius = imsize / 2
        if pixels is None:
            pixel_scale = self.SKYVIEW_PIXEL_SCALES.get(skyview_name)
            if pixel_scale is not None:
                pixels = int(round(imsize.to(u.arcsec).value / pixel_scale))
        images = SkyView.get_images(
            position=self.coord, survey=skyview_name, radius=radius,
            pixels=str(pixels) if pixels is not None else None,
        )
        if not images or not images[0]:
            warnings.warn(f"SkyView returned no image for {skyview_name}.")
            return None
        return images[0][0]

    def get_image(self, imsize=None, band=None, radius=None, pixels=None):
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
            img_hdu = self.get_gleam(imsize, pixels=pixels)
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

    def get_cutout(self, imsize=None, band=None, radius=None, pixels=None):
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

    def get_first(self, imsize, pixels=None):
        return self._skyview_fetch('VLA FIRST (1.4 GHz)', imsize, pixels=pixels)

    def get_nvss(self, imsize, pixels=None):
        return self._skyview_fetch('NVSS', imsize, pixels=pixels)

    def get_wenss(self, imsize, pixels=None):
        return self._skyview_fetch('WENSS', imsize, pixels=pixels)

    def get_gleam(self, imsize, band='170-231 MHz', pixels=None):
        return self._skyview_fetch(f'GLEAM {band}', imsize, pixels=pixels)

    def get_tgss(self, imsize, pixels=None):
        return self._skyview_fetch('TGSS ADR1', imsize, pixels=pixels)

    def get_sdss(self, imsize, band=None, pixels=None):
        band = 'r' if band is None else band.lower()
        if band not in self.SDSS_SURVEYS:
            raise TypeError(f"Allowed filters for SDSS are {list(self.SDSS_SURVEYS)}")
        return self._skyview_fetch(self.SDSS_SURVEYS[band], imsize, pixels=pixels)

    def get_galex(self, imsize, band=None, pixels=None):
        band = 'NUV' if band is None else band.upper()
        if band not in self.GALEX_SURVEYS:
            raise TypeError(f"Allowed filters for GALEX are {list(self.GALEX_SURVEYS)}")
        return self._skyview_fetch(self.GALEX_SURVEYS[band], imsize, pixels=pixels)

    def get_twomass(self, imsize, band=None, pixels=None):
        band = 'J' if band is None else band.upper()
        if band not in self.TWOMASS_SURVEYS:
            raise TypeError(f"Allowed filters for 2MASS are {list(self.TWOMASS_SURVEYS)}")
        return self._skyview_fetch(self.TWOMASS_SURVEYS[band], imsize, pixels=pixels)
