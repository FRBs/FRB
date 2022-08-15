""" Module for EM calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

from astropy import units
from astropy import constants


def em_from_halpha(sb_obs, z, T=1e4*units.K):
    """
    Estimate EM from the observed Halpha surface brightness

    Follows the Reynolds 1977 formalism

    Args:
        sb_obs (Quantity):
            Observed surface brightness
        z (float):
            Redshift of the galaxy
        T (Quantity, optional):
            Temperature for the analysis

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
    EM_Ha = 2.75 * units.pc / units.cm**6 * (T.to('K').value/1e4)**0.9 * I_R.to('rayleigh').value

    # Return
    return EM_Ha


def dm_from_em(EM, L, ff=1., eps=1., cloudcloud=2.):
    """
    This follows the formalism presented in Tendulkar+2017
    which follows Reynolds 1977 and Cordes+2016

    Args:
        EM (Quantity):
          Emission measure
        L (Quantity):
          Linear size of the source
        ff (float, optional):
          Filling factor
        eps (float, optional):
             fractional variation inside discrete clouds 
             due to turbulent-like density variations
        cloudcloud (float, optional):
            cloud-to-cloud density variations in the ionized region of depth 
            L in kpc.

    Returns:
        Quantity: DM at the source;  correct for (1+z)^-1 at your liking
    """
    # DM at the source
    DM_s = 387 * units.pc / units.cm**3 * np.sqrt(L.to('kpc').value) * np.sqrt(
        ff/(cloudcloud * (1+eps**2)/4)) * np.sqrt(EM.to('pc/cm**6').value/600)
    # Return
    return DM_s
