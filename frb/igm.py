""" Module for IGM calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from pkg_resources import resource_filename

from scipy.interpolate import interp1d

from astropy import units
from astropy.table import Table
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
    Return mass density in stars as calculated by
    Madau & Dickinson (2014)

    Parameters
    ----------
    z

    Returns
    -------
    rho_Mstar: float or ndarray

    """
    # Load
    stellar_mass_file = resource_filename('frb', 'data/IGM/stellarmass.dat')
    rho_mstar_tbl = Table.read(stellar_mass_file, format='ascii')
    # float or ndarray?
    if not isiterable(z):
        z = np.array([z])
        flg_z = 0
    else:
        flg_z = 1
    # Output
    rho_Mstar_unitless = np.zeros_like(z)

    # Extrema
    highz = z > rho_mstar_tbl['z'][-1]
    rho_Mstar_unitless[highz] = rho_mstar_tbl['rho_Mstar'][-1]

    # Interpolate
    fint = interp1d(rho_mstar_tbl['z'], rho_mstar_tbl['rho_Mstar'], kind='cubic')
    rho_Mstar_unitless[~highz] = fint(z[~highz])

    # Finish
    rho_Mstar = rho_Mstar_unitless * units.Msun / units.Mpc**3
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
