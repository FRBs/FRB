""" Module for EM calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import os

from pkg_resources import resource_filename

from astropy import units
from astropy import constants

def em_from_halpha(sb_obs, z, T=1e4*units.K):
    """
    Estimate EM from the observed Halpha surface brightness

    Follows the Reynolds 1977 formalism

    Args:

    Returns:
        Quantity: EM

    """
    # Correct for surface brightness dimming
    sb_corr = sb_obs * (1+z)**4

    # Halpha energy
    Ha = 6564.613 * units.Angstrom
    E_Ha_photon = constants.c * constants.h / Ha

    # Intensity in Rayleighs
    I_R = (sb_corr * units.ph / E_Ha_photon).to('rayleigh')

    # EM
    EM_Ha = 2.75 * units.pc / units.cm**6 * (T.value/1e4)**0.9 * I_R.to('rayleigh').value

    # Return
    return EM_Ha

