""" PSRCat survey """

import numpy as np

from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord
from astropy import units

try:
    from pulsars import io as pio
except ImportError:
    print("Warning:  You need FRB/pulsars installed to use PSRCat")
    pio = None

from frb.surveys import surveycoord
from frb.surveys import catalog_utils

    
class PSRCAT_Survey(surveycoord.SurveyCoord):
    """
    Class to handle queries on the PSRCAT catalog

    Args:
        coord (astropy.coordinates.SkyCoord): Coordinate for surveying around
        radius (astropy.coordinates.Angle): Search radius around the coordinate
        **kwargs: Passed to :class:`frb.surveys.surveycoord.SurveyCoord`
            (e.g. ``verbose``)

    """
    def __init__(self, coord: SkyCoord, radius: Angle, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        #
        self.survey = 'PSRCAT'

    def get_catalog(self) -> Table:
        """
        Grab the catalog of pulsars around the input coordinate to the search radius

        Returns:
            astropy.table.Table:  Catalog of sources returned.
            Empty (with 'ra' and 'dec' columns) if there are no pulsars in the cone.

        Raises:
            ImportError: If the FRB/pulsars package is not installed.

        """
        if pio is None:
            raise ImportError("You need FRB/pulsars installed to use PSRCat")
        # Load em
        pulsars = pio.load_pulsars()

        # Coords
        pcoord = SkyCoord(pulsars['RAJ'], pulsars['DECJ'], unit=(units.hourangle, units.deg))

        # Query
        gdp = pcoord.separation(self.coord) <= self.radius

        if not np.any(gdp):
            self.catalog = catalog_utils.ensure_empty_schema(Table(), ['ra', 'dec'])
        else:
            catalog = pulsars[gdp]

            # Clean
            catalog['ra'] = pcoord[gdp].ra.value
            catalog['dec'] = pcoord[gdp].dec.value
            for key in ['ra', 'dec']:
                catalog[key].unit = units.deg
            # Sort
            self.catalog = catalog_utils.sort_by_separation(catalog, self.coord,
                                                            radec=('ra', 'dec'))
        # Add meta, etc.
        self.catalog.meta['radius'] = self.radius
        self.catalog.meta['survey'] = self.survey
        # Validate
        self.validate_catalog()
        # Return
        return self.catalog

