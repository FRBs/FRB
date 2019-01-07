""" Surveys to be accessed through the HEASARC interface (via astroquery"""

import pdb

from astropy.table import Table
from astropy import units, wcs

try:
    from astroquery.heasarc import Heasarc
    from astroquery.skyview import SkyView
except ImportError:
    print("Warning:  You need astroquery installed to use the surveys from HEASARC and SkyView")

from frb.surveys import surveycoord
from frb.surveys import catalog_utils


class HEASARC_Survey(surveycoord.SurveyCoord):
    def __init__(self, coord, radius, mission, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        #
        self.survey = None
        self.mission = mission
        # Instantiate astroquery object
        self.heasarc = Heasarc()

    def get_catalog(self):
        try:
            catalog = self.heasarc.query_region(self.coord,
                                                mission=self.mission,
                                                radius=self.radius)
        except (ValueError, TypeError):  # No table found
            self.catalog = Table()
        else:
            # Clean
            catalog.rename_column("RA", "ra")
            catalog.rename_column("DEC", "dec")
            for key in ['ra', 'dec']:
                catalog[key].unit = units.deg
            # Sort
            self.catalog = catalog_utils.sort_by_separation(catalog,
                                                            self.coord,
                                                            radec=('ra', 'dec'))
        # Add meta, etc.
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey
        # Validate
        self.validate_catalog()
        # Return
        return self.catalog


class SkyView_Survey(surveycoord.SurveyCoord):
    def __init__(self, coord, radius, mission, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        #
        self.survey = None
        self.mission = mission
        # Instantiate astroquery object
        self.skyview = SkyView()

    def get_cutout(self):
        if self.mission.lower() == 'first':
            img_hdu = self.get_first()
        elif self.mission.lower() == 'nvss':
            img_hdu = self.get_nvss()
        elif self.mission.lower() == 'wenss':
            img_hdu = self.get_wenss()
        elif self.mission.lower() == 'gleam':
            img_hdu = self.get_gleam()
        elif self.mission.lower() == 'tgss':
            img_hdu = self.get_tgss()

        self.cutout = img_hdu.data
        self.cutout_hdr = img_hdu.header

        mywcs = wcs.WCS(self.cutout_hdr)
        ypix, xpix = self.cutout.shape
        (ra0, dec0), (ra1, dec1), = mywcs.wcs_pix2world([[0, 0], [xpix, ypix]],
                                                        0)
        print("Got image spanning (RA, Dec) = ({0} - {1}, {2} - {3})"
              .format(ra0, ra1, dec0, dec1))

        return self.cutout

    def get_first(self):
        return SkyView.get_images(position=self.coord,
                                  survey='VLA FIRST (1.4 GHz)',
                                  radius=self.radius)[0][0]

    def get_nvss(self):
        return SkyView.get_images(position=self.coord, survey='NVSS',
                                  radius=self.radius)[0][0]

    def get_wenss(self):
        return SkyView.get_images(position=self.coord, survey='WENSS',
                                  radius=self.radius)[0][0]

    def get_gleam(self, band="170-231 MHz"):
        return SkyView.get_images(position=self.coord,
                                  survey='GLEAM {0}'.format(band),
                                  radius=self.radius)[0][0]

    def get_tgss(self):
        return SkyView.get_images(position=self.coord,
                                  survey='TGSS ADR1',
                                  radius=self.radius)[0][0]


class NVSS_Survey(HEASARC_Survey, SkyView_Survey):
    def __init__(self, coord, radius, **kwargs):
        HEASARC_Survey.__init__(self, coord, radius, 'nvss', **kwargs)
        SkyView_Survey.__init__(self, coord, radius, 'nvss', **kwargs)
        self.survey = 'NVSS'


class FIRST_Survey(HEASARC_Survey, SkyView_Survey):
    def __init__(self, coord, radius, **kwargs):
        HEASARC_Survey.__init__(self, coord, radius, 'first', **kwargs)
        SkyView_Survey.__init__(self, coord, radius, 'first', **kwargs)
        self.survey = 'FIRST'


class WENSS_Survey(HEASARC_Survey, SkyView_Survey):
    def __init__(self, coord, radius, **kwargs):
        HEASARC_Survey.__init__(self, coord, radius, 'wenss', **kwargs)
        SkyView_Survey.__init__(self, coord, radius, 'wenss', **kwargs)
        self.survey = 'WENSS'
