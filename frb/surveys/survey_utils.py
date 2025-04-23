""" utils related to SurveyCoord objects"""

from urllib.error import HTTPError
from frb.surveys.nedlvs import NEDLVS
from frb.surveys.sdss import SDSS_Survey
from frb.surveys.des import DES_Survey
from frb.surveys.wise import WISE_Survey
from frb.surveys.decals import DECaL_Survey
from frb.surveys.psrcat import PSRCAT_Survey
from frb.surveys import heasarc
from frb.surveys.panstarrs import Pan_STARRS_Survey
from frb.surveys.nsc import NSC_Survey
from frb.surveys.delve import DELVE_Survey
from frb.surveys.vista import VISTA_Survey
from frb.surveys.cluster_search import TullyGroupCat
from frb.surveys.galex import GALEX_Survey
from frb.surveys.hsc import HSC_Survey, QueryError
from frb.surveys.catalog_utils import xmatch_and_merge_cats, remove_duplicates

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table, join
from pyvo.dal import DALServiceError
from requests import ReadTimeout, HTTPError

import numpy as np
import warnings

optical_surveys = ['Pan-STARRS', 'GALEX','SDSS', 'DES', 'DELVE',  'DECaL', 'VISTA', 'HSC', 'NEDLVS', 'WISE'] # 'NSC',
group_catalogs = ['TullyGroupCat']
radio_surveys = ['NVSS', 'FIRST', 'WENSS', 'PSRCAT']
allowed_surveys = optical_surveys+radio_surveys+group_catalogs


def load_survey_by_name(name, coord, radius, **kwargs):
    """
    Load up a Survey class object for the named survey
    allowed_surveys = ['SDSS', 'DES', 'NVSS', 'FIRST', 'WENSS', 'DECaL', 
    'PSRCAT', 'WISE', 'Pan-STARRS', 'NEDLVS']

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
    elif name == 'DELVE':
        survey = DELVE_Survey(coord, radius, **kwargs)
    elif name == 'VISTA':
        survey = VISTA_Survey(coord, radius, **kwargs)
    elif name == 'HSC':
        survey = HSC_Survey(coord, radius, **kwargs)
    elif name == 'NEDLVS':
        survey = NEDLVS(coord, radius, **kwargs)
    elif name == 'TullyGroupCat':
        survey = TullyGroupCat(coord, radius, **kwargs)
    elif name == 'GALEX':
        survey = GALEX_Survey(coord, radius, **kwargs)

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
    if surveyname == "NEDLVS":
        radius = 1*u.deg # Not as deep as the rest and hence a larger radius.
    else:
        radius = 1*u.arcmin
    survey = load_survey_by_name(surveyname, coord, radius)
    try:
        cat = survey.get_catalog()
    except DALServiceError:
        warnings.warn("Couldn't reach IPAC.", RuntimeWarning)
        cat = None
    except ReadTimeout:
        warnings.warn("Couldn't reach NOIRLAB.", RuntimeWarning)
        cat = None
    except QueryError:
        warnings.warn("Do not have credentials to search HSC.", RuntimeWarning)
        cat = None
    except HTTPError:
        warnings.warn("Couldn't reach MAST for PS1.", RuntimeWarning)
        cat = None
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
        all_surveys = optical_surveys+radio_surveys
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
        surveys = optical_surveys+radio_surveys
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
        except (ConnectionError, HTTPError, QueryError):
            warnings.warn("Couldn't connect to {:s}. Skipping this for now.".format(surveyname), RuntimeWarning)

        # Did the survey return something?
        if (survey.catalog is not None):
            if len(survey.catalog)>0:
                print("Found {:d} objects in {:s}".format(len(survey.catalog), surveyname))
                if len(combined_cat)==0:
                    # First time
                    combined_cat = survey.catalog
                else:
                    # Combine otherwise
                    # Deal with duplicate columns first
                    # remove separations
                    if 'separation' in survey.catalog.colnames:
                        survey.catalog.remove_column('separation')
                    # Now other duplicates
                    duplicate_colnames = np.array(survey.catalog.colnames)[np.isin(survey.catalog.colnames, combined_cat.colnames)]
                    # Ignore ra, dec
                    duplicate_colnames = duplicate_colnames[~np.isin(duplicate_colnames, ['ra', 'dec'])]
                    # Rename them
                    if len(duplicate_colnames)>0:
                        renamed_duplicates = [colname+"_"+surveyname for colname in duplicate_colnames]
                        survey.catalog.rename_columns(duplicate_colnames.tolist(), renamed_duplicates)

                    # Now merge
                    combined_cat = xmatch_and_merge_cats(combined_cat, survey.catalog)
            # No objects found?
            elif len(survey.catalog)==0:
                print("Empty table in "+surveyname)
    if len(combined_cat)>0:
        # Fill in any empty separations and sort them.
        combined_cat['separation'] = coord.separation(SkyCoord(combined_cat['ra'], combined_cat['dec'], unit='deg')).to(u.arcmin)
        combined_cat.sort('separation')

        # Make the ra, dec, separation the first columns
        colnames = combined_cat.colnames
        other_cols = np.setdiff1d(colnames, ['ra', 'dec', 'separation'])
        combined_cat = combined_cat[['ra', 'dec', 'separation']+other_cols.tolist()]
    
    return combined_cat
           
def PS1_tile(coord:SkyCoord, side:u.Quantity=1*u.deg, **kwargs)->Table:
    """
    Tile multiple 20' cone searches of 
    the Pan-STARRS catalog to cover a larger
    roughly square region in the sky. Doing
    this manually because MAST doesn't
    allow cone searches of radius greater than
    30'.
    Args:
        coord (SkyCoord): Center of search region.
        side (astropy Quantity): Angular size
            of the square region to be searched.
        kwargs: Additional keyword arguments
            to be passed onto the Pan-STARRS_Survey
            get_catalog method.
    Returns:
        combined_tab (Table): A PS1 catalog.
    """
    assert side>=30*u.arcmin, "Use a regular Pan-STARRS search for this radius."
    RA0, DEC0 = coord.ra, coord.dec

    radius = 20*u.arcmin

    # Make a grid of RA, Dec to search
    n = 2*int((side/2)//(radius*np.sqrt(2)))+1
    ra_vals = RA0 + radius*np.sqrt(2)*np.linspace(-1,1,n)
    dec_vals = DEC0 + radius*np.sqrt(2)*np.linspace(-1,1,n)

    # Loop over grid points, search and collate catalogs
    print("Looping through {:d} pointings.".format(n*n))
    combined_table = Table()
    for i, ra in enumerate(ra_vals):
        for j, dec in enumerate(dec_vals):
            print("Getting pointing #",i*len(ra_vals)+j+1)
            coord = SkyCoord(ra,dec)
            survey = load_survey_by_name("Pan-STARRS", coord, radius=radius)
            cat = survey.get_catalog(**kwargs)
            # Remove columns so that
            # the join function further down the 
            # line has no trouble merging duplicates.
            cat.remove_column('separation')
            if len(cat)>0:
                if len(combined_table)>0:
                    combined_table = join(combined_table, cat, join_type='outer')
                else:
                    combined_table = cat
            else:
                continue
    combined_table = remove_duplicates(combined_table, 'Pan-STARRS_ID')
    return combined_table