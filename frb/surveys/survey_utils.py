""" utils related to SurveyCoord objects"""

from frb.surveys.sdss import SDSS_Survey
from frb.surveys.des import DES_Survey


def load_survey_by_name(name, coord, radius, **kwargs):
    allowed_surveys = ['SDSS', 'DES']

    if name not in allowed_surveys:
        raise IOError("Not ready for input survey.\n These are allowed: {}".format(allowed_surveys))

    # Do it
    if name == 'SDSS':
        survey = SDSS_Survey(coord, radius, **kwargs)
    elif name == 'DES':
        survey = DES_Survey(coord, radius, **kwargs)

    # Return
    return survey


