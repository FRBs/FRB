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
from astropy import constants

from frb.io import load_dla_fits
from frb.turb_scattering import Turbulence

# Fukugita 2004 (Table 1)
def fukugita04_dict():
    f04_dict = {}
    f04_dict['M_sphere'] = 0.0015
    f04_dict['M_disk'] = 0.00055
    f04_dict['M_HI'] = 0.00062
    f04_dict['M_H2'] = 0.00016
    f04_dict['M_WD'] = 0.00036
    f04_dict['M_NS'] = 0.00005
    f04_dict['M_BH'] = 0.00007
    f04_dict['M_BD'] = 0.00014
    # Return
    return f04_dict

def average_fHI(z, z_reion=7.):
    """
    Average HI fraction


    Parameters
    ----------
    z
    z_reion

    Returns
    -------

    """
    z, flg_z = z_to_array(z)
    fHI = np.zeros_like(z)
    #
    zion = z > z_reion
    fHI[zion] = 1.
    # Return
    if flg_z:
        return fHI
    else:
        return fHI[0]

def average_DM(z, cosmo=None, cumul=False, neval=10000,
               mu=1.3):
    """
    Calculate the average DM 'expected' based on our empirical
    knowledge of baryon distributions and their ionization state.

    Parameters
    ----------
    z: float
    mu: float
      Reduced mass correction for He when calculating n_H
    cumul: bool, optional
      Return the DM as a function of z

    Returns
    -------
    DM: Quantity
      One value if cumul=False
      else evaluated at a series of z (neval)
    zeval: ndarray
      Evaluation redshifts if cumul=True

    """
    # Cosmology
    if cosmo is None:
        cosmo = Planck15
    # Init
    zeval = np.linspace(0., z, neval)

    # Get baryon mass density
    rho_b = cosmo.Ob(zeval) * cosmo.critical_density0.to('Msun/Mpc**3')

    # Dense components
    rho_Mstar = avg_rhoMstar(zeval, remnants=True)
    rho_ISM = avg_rhoISM(zeval)

    # Diffuse
    rho_diffuse = rho_b - (rho_Mstar+rho_ISM)

    # Here we go
    n_H = (rho_diffuse/constants.m_p/mu)
    n_He = n_H / 12.  # 25% He mass fraction

    pdb.set_trace()



def avg_rhoISM(z):
    """
    Co-moving Mass density of the ISM

    Assumes z=0 properties for z<1
    and otherwise M_ISM = M* for z>1

    Parameters
    ----------
    z

    Returns
    -------
    rhoISM : Quantity
      Units of Msun/Mpc^3

    """
    # Init
    z, flg_z = z_to_array(z)
    rhoISM_unitless = np.zeros_like(z)
    # Mstar
    rhoMstar = avg_rhoMstar(z, remnants=False)
    # z<1
    f04_dict = fukugita04_dict()
    M_ISM = f04_dict['M_HI'] + f04_dict['M_H2']
    f_ISM = M_ISM/(f04_dict['M_sphere']+f04_dict['M_disk'])
    lowz = z<1
    rhoISM_unitless[lowz] = f_ISM * rhoMstar[lowz].value
    # z>1
    rhoISM_unitless[~lowz] = rhoMstar[~lowz].value
    # Finish
    rhoISM = rhoISM_unitless * units.Msun / units.Mpc**3
    #
    return rhoISM



def avg_rhoMstar(z, remnants=True):
    """
    Return mass density in stars as calculated by
    Madau & Dickinson (2014)

    Parameters
    ----------
    z: float or ndarray
    remnants: bool, optional
      Include remnants in the calculation?

    Returns
    -------
    rho_Mstar: Quantity
      Units of Msun/Mpc^3

    """
    # Init
    z, flg_z = z_to_array(z)
    # Load
    stellar_mass_file = resource_filename('frb', 'data/IGM/stellarmass.dat')
    rho_mstar_tbl = Table.read(stellar_mass_file, format='ascii')
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

    # Remnants
    if remnants:
        # Fukugita 2004 (Table 1)
        f04_dict = fukugita04_dict()
        f_remnants = (f04_dict['M_WD'] + f04_dict['M_NS'] + f04_dict['M_BH'] + f04_dict['M_BD']) / (
                f04_dict['M_sphere'] + f04_dict['M_disk'])
        # Apply
        rho_Mstar *= (1+f_remnants)

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
    SFR: Quantity
      Units of Msun/yr/Mpc^3

    """
    rho_SFR_unitless = 0.015 * (1+z)**2.7 / (1 + ((1+z)/2.9)**5.6)
    rho_SFR = rho_SFR_unitless * units.Msun / units.yr / units.Mpc**3

    # Return
    return rho_SFR

def z_to_array(z):
    # float or ndarray?
    if not isiterable(z):
        z = np.array([z])
        flg_z = 0
    else:
        flg_z = 1
    # Return
    return z, flg_z
