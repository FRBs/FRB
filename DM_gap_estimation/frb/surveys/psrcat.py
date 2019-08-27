""" PSRCat survey """

import pdb

import numpy as np

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units

try:
    from pulsars import io as pio
except ImportError:
    print("Warning:  You need FRB/pulsars installed to use PSRCat")

from frb.surveys import surveycoord
from frb.surveys import catalog_utils

    
class PSRCAT_Survey(surveycoord.SurveyCoord):
    """
    Class to handle queries on the PSRCAT catalog
    
    Args:
        coord (SkyCoord): Coordiante for surveying around
        radius (Angle): Search radius around the coordinate
        
    """
    def __init__(self, coord, radius, **kwargs):
        surveycoord.SurveyCoord.__init__(self, coord, radius, **kwargs)
        #
        self.survey = 'PSRCAT'

    def get_catalog(self):
        """
        Grab the catalog of pulsars around the input coordinate to the search radius
        
        Returns:
            astropy.table.Table:  Catalog of sources returned

        """
        # Load em
        pulsars = pio.load_pulsars()

        # Coords
        pcoord = SkyCoord(pulsars['RAJ'], pulsars['DECJ'], unit=(units.hourangle, units.deg))

        # Query
        gdp = pcoord.separation(self.coord) <= self.radius

        if not np.any(gdp):
            self.catalog = Table()
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

