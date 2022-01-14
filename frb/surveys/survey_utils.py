""" utils related to SurveyCoord objects"""

from frb.surveys.sdss import SDSS_Survey
from frb.surveys.des import DES_Survey
from frb.surveys.wise import WISE_Survey
from frb.surveys.decals import DECaL_Survey
from frb.surveys.psrcat import PSRCAT_Survey
from frb.surveys import heasarc
from frb.surveys.panstarrs import Pan_STARRS_Survey
from frb.surveys.nsc import NSC_Survey
from frb.surveys.vista import VISTA_Survey

from astropy.coordinates import SkyCoord
from astropy import units as u

import numpy as np
import warnings

allowed_surveys = ['SDSS', 'DES', 'NVSS', 'FIRST', 'WENSS', 'DECaL', 'PSRCAT', 'WISE', 'Pan-STARRS', 'NSC', 'VISTA']


def load_survey_by_name(name, coord, radius, **kwargs):
    """
    Load up a Survey class object for the named survey
    allowed_surveys = ['SDSS', 'DES', 'NVSS', 'FIRST', 'WENSS', 'DECaL', 'PSRCAT', 'WISE', 'Pan-STARRS']

    Args:
        name (str): Name of the survey 
        coord (astropy.coordiantes.SkyCoord): Coordinate to define survey around 
        radius (astropy.units.Quanity): Outer radius of the survey
        **kwargs: Passed the Survey object

    Returns:
        frb.surveys.SurveyCoord: Child of this parent given by input survey name

    """

    if name not in allowed_surveys:
        raise IOError("Not ready for input survey.\n These are allowed: {}".format(allowed_surveys))

    # Do it
    if name == 'SDSS':
        survey = SDSS_Survey(coord, radius, **kwargs)
    elif name == 'DES':
        survey = DES_Survey(coord, radius, **kwargs)
    elif name == 'NVSS':
        survey = heasarc.NVSS_Survey(coord, radius, **kwargs)
    elif name == 'WENSS':
        survey = heasarc.WENSS_Survey(coord, radius, **kwargs)
    elif name == 'FIRST':
        survey = heasarc.FIRST_Survey(coord, radius, **kwargs)
    elif name == 'DECaL':
        survey = DECaL_Survey(coord, radius, **kwargs)
    elif name == 'WISE':
        survey = WISE_Survey(coord, radius, **kwargs)
    elif name == 'PSRCAT':
        survey = PSRCAT_Survey(coord, radius, **kwargs)
    elif name == 'Pan-STARRS':
        survey = Pan_STARRS_Survey(coord, radius,**kwargs)
    elif name == 'NSC':
        survey = NSC_Survey(coord, radius, **kwargs)
    elif name == 'VISTA':
        survey = VISTA_Survey(coord, radius, **kwargs)

    # Return
    return survey

def is_inside(surveyname:str, coord:SkyCoord)->bool:
    """
    Tests if a coordinate is within a survey footprint.
    Args:
        surveyname (str): Name of the survey
        coord (astropy.coordiantes.SkyCoord): Coordinate to check
    Returns:
        inside (bool): True if coord is within the footprint.
    """

    # Instantiate survey and run a cone search with 1 arcmin radius
    survey = load_survey_by_name(surveyname, coord, 1*u.arcmin)
    cat = survey.get_catalog()

    # Are there any objects in the returned catalog?
    if cat is None or len(cat) == 0:
        return False
    # is the coordinate in a region that's close to the edge but on the exterior?
    else:
        # Compute PAs of all returned objects
        cat_coords = SkyCoord(cat['ra'], cat['dec'], unit="deg")
        pa_list = coord.position_angle(cat_coords)

        # If the object is in the interior, then the max PA
        # difference is likely greater than 180 degrees. 
        max_pa_diff = np.max(pa_list) - np.min(pa_list)
        if np.abs(max_pa_diff) <= 180*u.deg:
            # Are there too few objects? Warn the user if so.
            if len(cat)<=5:
                # Make sure there is no source at the search location.
                # if that's not the case, then no need to warn.
                separations = coord.separation(cat_coords)
                if np.min(separations)>2*u.arcsec:
                    warnings.warn("Only {} objects found in {} within 1'. Check location manually.".format(len(cat), surveyname),
                                  RuntimeWarning, stacklevel=2)
            return False
        else:
            return True

def in_which_survey(coord:SkyCoord)->dict:
    """
    Check if a particular coord is inside any
    survey that can be currently queried from
    `frb.surveys` module.
    Args:
        coord (astropy.coordiantes.SkyCoord): Coordinate to check
    Returns:
        inside (dict): A dict which tells which surveys the coordinate
            is inside.
    """
    # Loop through known surveys and check them one by one.
    inside = {}
    for surveyname in allowed_surveys:
        # Skip PSRCAT
        if surveyname == "PSRCAT":
            continue
        inside[surveyname] = is_inside(surveyname, coord)
    
    return inside
