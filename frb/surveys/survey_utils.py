""" utils related to SurveyCoord objects"""

from urllib.error import HTTPError
from frb.surveys.sdss import SDSS_Survey
from frb.surveys.des import DES_Survey
from frb.surveys.wise import WISE_Survey
from frb.surveys.decals import DECaL_Survey
from frb.surveys.psrcat import PSRCAT_Survey
from frb.surveys import heasarc
from frb.surveys.panstarrs import Pan_STARRS_Survey
from frb.surveys.nsc import NSC_Survey
from frb.surveys.vista import VISTA_Survey
from frb.surveys.catalog_utils import xmatch_and_merge_cats

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

import numpy as np
import warnings

optical_surveys = ['Pan-STARRS', 'WISE', 'SDSS', 'DES',  'DECaL', 'VISTA', 'NSC']
radio_surveys = ['NVSS', 'FIRST', 'WENSS', 'PSRCAT']
allowed_surveys = optical_surveys+radio_surveys


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

def in_which_survey(coord:SkyCoord, optical_only:bool=True)->dict:
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
    if optical_only:
        all_surveys = optical_surveys
    else:
        all_surveys = allowed_surveys
    for surveyname in all_surveys:
        # Skip PSRCAT
        if surveyname == "PSRCAT":
            continue
        inside[surveyname] = is_inside(surveyname, coord)
    
    return inside


def search_all_surveys(coord:SkyCoord, radius:u.Quantity, include_radio:bool=False,
                       seed_cat:Table=None):
    """
    A method to query all allowed surveys and combine
    the results into one table.
    Args:
        coord (SkyCoord): Central coordinates of cone search.
        radius (Quantity): Search radius in angular units.
        include_radio (bool, optional): Want to include results from the HEASARC surveys?
            Include at your own risk. Untested. Might break in unexpected ways.
        seed_cat (Table, optional): If you'd like to merge the survey results
            with another photometry table that you already have.

    Returns:
        combined_cat (Table): Table of merged query results.
    """

    # Start with the seed table
    if seed_cat is not None:
        assert isinstance(seed_cat, Table), "The seed_cat must be an astropy table"
        assert np.all(np.isin(['ra','dec'],seed_cat.colnames)), "The seed catalog doesn't have 'ra' and/or 'dec' columns." 
        combined_cat = seed_cat
    else:
        # Start with an empty table
        combined_cat = Table()
    # Select surveys
    if include_radio: # Careful! NOT TESTED!
        surveys = allowed_surveys 
    else:
        surveys = optical_surveys
    # Loop over them
    for surveyname in surveys:
        if surveyname=='Pan-STARRS':
            if radius>0.5*u.deg:
                warnings.warn("Pan=STARRS doesn't allow cone searches wider than 0.5 deg. Skipping.", RuntimeWarning)
                continue
        survey = load_survey_by_name(name=surveyname, coord=coord, radius=radius)
        try:
            survey.get_catalog()
        except ConnectionError:
            warnings.warn("Couldn't connect to {:s}. Skipping this for now.".format(surveyname), RuntimeWarning)
        except HTTPError:
            warnings.warn("Couldn't connect to {:s}. Skipping this for now.".format(surveyname), RuntimeWarning)
        # Did the survey return something?
        if (survey.catalog is not None) & (len(survey.catalog)>0):
            print("Found {:d} objects in {:s}".format(len(survey.catalog), surveyname))
            if len(combined_cat)==0:
                # First time
                combined_cat = survey.catalog
            else:
                # Combine otherwise
                # TODO: Need to deal with duplicate column names more elegantly.
                combined_cat = xmatch_and_merge_cats(combined_cat, survey.catalog,)
        # No objects found?
        elif len(survey.catalog)==0:
            print("Empty table in "+surveyname)
    
    return combined_cat

            
