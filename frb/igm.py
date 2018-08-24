""" Module for IGM calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from scipy.interpolate import interp1d

from astropy import units
from astropy import constants
from astropy.utils import isiterable
from astropy.cosmology import Planck15

from frb.io import load_dla_fits
from frb.turb_scattering import Turbulence


def average_DM(z, cosmo=None):
    """
    Calculate the average DM 'expected' based on our empirical
    knowledge of baryon distributions and their ionization state.

    Parameters
    ----------
    z: float or ndarray

    Returns
    -------
    DM: float or ndarray
      mirrors z type

    """
    # Cosmology
    if cosmo is None:
        cosmo = Planck15
    #

def avg_rhoMstar(z):
    """
    Kludging for now

    Parameters
    ----------
    z

    Returns
    -------
    rho_Mstar: float or ndarray

    """
    if not isiterable(z):
        z = np.array([z])
        flg_z = 0
    else:
        flg_z = 1
    # Eye-ball of Madau & Dickinson (2014)
    zval = [0., 0.5, 1., 1.5, 2., 3., 4.]
    logrho_val = [8.8, 8.7, 8.55, 8.35, 8.15, 7.7, 7.3]
    # Output
    logrho_Mstar_unitless = np.zeros_like(z)

    # Extrema
    highz = z > zval[-1]
    logrho_Mstar_unitless[highz] = logrho_val[-1]

    # Interpolate
    fint = interp1d(zval, logrho_val, kind='cubic')
    logrho_Mstar_unitless[~highz] = fint(z[~highz])

    # Finish
    rho_Mstar = 10**logrho_Mstar_unitless * units.Msun / units.Mpc**3

    # Return
    if flg_z:
        return rho_Mstar
    else:
        return rho_Mstar[0]


def avg_rhoSFR(z):
    """
    Average SFR density

    Based on Madau & Dickinson (2014)

    Parameters
    ----------
    z: float or ndarray

    Returns
    -------
    SFR: float or ndarray
      Units of Msun/yr/Mpc^3

    """
    rho_SFR_unitless = 0.015 * (1+z)**2.7 / (1 + ((1+z)/2.9)**5.6)
    rho_SFR = rho_SFR_unitless * units.Msun / units.yr / units.Mpc**3

    # Return
    return rho_SFR
