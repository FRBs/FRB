""" Module for IGM calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
from IPython import embed

from pkg_resources import resource_filename

from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from astropy import units
from astropy.table import Table
from astropy.utils import isiterable
from astropy.cosmology import Planck15
from astropy import constants

from frb import halos
from frb import mw

def fukugita04_dict():
    """
    Data from Fukugita 2004, Table 1

    Returns
    -------
    f04_dict: dict

    """
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
    1 = neutral
    0 = fully ionized


    Parameters
    ----------
    z
    z_reion

    Returns
    -------
    fHI: float or ndarray
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


def average_He_nume(z, z_HIreion=7.):
    """
    Average number of electrons contributed by He as a function of redshift
    per He nucleus

    Following Kulkarni, Worseck & Hennawi 2018
        https://arxiv.org/abs/1807.09774


    Parameters
    ----------
    z: float or ndarray
      Redshift

    Returns
    -------
    neHe: ndarray
      Number of free electrons per Helium nucelus

    """
    z, flg_z = z_to_array(z)
    # Load Kulkarni Table
    He_file = resource_filename('frb', os.path.join('data','IGM','qheIII.txt'))
    qHeIII = Table.read(He_file, format='ascii')
    # Fully re-ionized
    first_ionized = np.where(qHeIII['Q_HeIII_18'] >= 1.)[0][0]
    z_HeIIreion = qHeIII['z'][first_ionized]
    #
    fHeI = np.zeros_like(z)
    fHeII = np.zeros_like(z)
    # HeI ionized at HI reionization
    zion = z > z_HIreion
    fHeI[zion] = 1.
    # HeII ionized at HeII reionization
    zion2 = (z > z_HeIIreion) & (z < z_HIreion)
    fi_HeIII = interp1d(qHeIII['z'], qHeIII['Q_HeIII_18'])
    fHeII[zion2] = 1. - fi_HeIII(z[zion2])
    # Combine
    neHe = (1.-fHeI) + (1.-fHeII)  #  No 2 on the second term as the first one gives you the first electron
    # Return
    if flg_z:
        return neHe
    else:
        return neHe[0]


def z_from_DM(DM, cosmo=None, coord=None, corr_nuisance=True):
    """
    Report back an estimated redshift from an input IGM DM
    Any contributions from the Galaxy and/or host need to have been 'removed'

    Args:
        DM (Quantity):
        cosmo (astropy.cosmology, optional):
        coord (astropy.coordinate.SkyCoord, optional):
           If provided, use it to remove the ISM
        corr_nuisance (bool, optional):
            If True, correct for the MW Halo and the host with
            100 DM units

    Returns:
        float: Redshift

    """
    if coord is not None:
        DM_ISM = mw.ismDM(coord)
        DM_use = DM - DM_ISM
    else:
        DM_use = DM

    # Correct
    if corr_nuisance:
        DM_use -= 100 * units.pc/units.cm**3

    # Calculate DMs
    all_DM, zeval = average_DM(20., cosmo=cosmo, neval=20000, cumul=True)
    # Interpolate
    fint = interp1d(all_DM.value, zeval)
    # Evaluate
    z = fint(DM_use.to('pc/cm**3').value)
    # Return
    return z


def average_DM(z, cosmo=None, cumul=False, neval=10000, mu=1.3):
    """
    Calculate the average cosmic DM 'expected' based on our empirical
    knowledge of baryon distributions and their ionization state.

    This includes both the IGM and galactic halos, i.e. any and all diffuse gas

    Args:
        z: float
          Redshift
        mu: float
          Reduced mass correction for He when calculating n_H
        cumul: bool, optional
          Return the DM as a function of z

    Returns:
        Quantity (cumul=False) or Quantity, ndarray (cumul=True): DM, zeval
    """
    # Cosmology
    if cosmo is None:
        cosmo = Planck15
    # Init
    zeval = np.linspace(0., z, neval)[1:]

    # Get baryon mass density
    rho_b = cosmo.Ob0 * cosmo.critical_density0.to('Msun/Mpc**3') * (1+zeval)**3

    # Dense components
    rho_Mstar = avg_rhoMstar(zeval, remnants=True)
    rho_ISM = avg_rhoISM(zeval)

    # Diffuse
    rho_diffuse = rho_b - (rho_Mstar+rho_ISM)*(1+zeval)**3

    # Here we go
    n_H = (rho_diffuse/constants.m_p/mu).to('cm**-3')
    n_He = n_H / 12.  # 25% He mass fraction

    n_e = n_H * (1.-average_fHI(zeval)) + n_He*(average_He_nume(zeval))

    # Cosmology -- 2nd term is the (1+z) factor for DM
    denom = np.sqrt((1+zeval)**3 * cosmo.Om0 + cosmo.Ode0) * (1+zeval) * (1+zeval)

    # Time to Sum
    dz = zeval[1] - zeval[0]
    DM_cum = ((constants.c/cosmo.H0) * np.cumsum(n_e * dz / denom)).to('pc/cm**3')

    # Return
    if cumul:
        return DM_cum, zeval
    else:
        return DM_cum[-1]


def avg_DMhalos(z, logMmin=10.3, f_diffuse=0.75, cumul=False, rmax=1.):
    """
    Average DM_halos term from halos along the sightline to an FRB

    Args:
        z: float
          Redshift of the FRB
        logMmin: float, optional
          Lowest mass halos to consider
          Cannot be much below 10.3 or the Halo code barfs
          The code deals with h^-1 factors, i.e. do not impose it yourself
        f_diffuse: float, optional
          Fraction of the cosmic baryon fraction in diffuse gas
        cumul: bool, optional
          Return a cumulative evaluation?
        rmax: float, optional
          Size of a halo in units of r200

    Returns:
        DM: Quantity
          One value if cumul=False
          else evaluated at a series of z
        zeval: ndarray, optional
          Evaluation redshifts if cumul=True
    """

    # Cosmic DM
    DM_cosmic, zeval = average_DM(z, cumul=True)
    # Halo mass fraction
    zvals = np.linspace(0., z, 20)
    fhalos = halos.frac_in_halos(zvals, 10**logMmin, 1e16, rmax = rmax)
    fhalos_interp = IUS(zvals, fhalos)

    # DM halos
    dDM = DM_cosmic.value - np.roll(DM_cosmic.value,1)
    dDM[0] = dDM[1]
    DM_halos = np.cumsum(dDM*fhalos_interp(zeval)*f_diffuse)*units.pc/units.cm**3

    # Return
    if cumul:
        return DM_halos, zeval
    else:
        return DM_halos[-1]


def average_DMIGM(z, fb=1., rmax=1., neval=10000):
    """
    Estimate DM_IGM in a cumulative fashion

    Args:
        z (float):
            Redshift for the evaluation
        fb (float, optional):
            Fraction of the cosmic baryons in a given halo
            1 means halos have all of their alloted baryons to r200
        rmax (float, optional):
        neval (int, optional):
            Number of redshift evaluations

    Returns:
        ndarray, Quantity

    """
    # Redshifts
    zvals = np.linspace(0., z, 50)  # Keep it cheap for fhalos
    dz_vals = zvals[1] - zvals[0]
    # Fraction of DM mass in halos
    fhalos = halos.frac_in_halos(zvals, 10**(10.3), 1e16, rmax=rmax)

    # Fraction of mass in IGM
    fIGM = 1. - fhalos*fb

    # DM_cosmic
    DM_cosmic_cumul, zeval = average_DM(z, cumul=True, neval=neval)
    dzeval = zeval[1] - zeval[0]
    dDM_cosmic = DM_cosmic_cumul - np.roll(DM_cosmic_cumul, 1)
    dDM_cosmic[0] = dDM_cosmic[1]
    DM_interp = IUS(zeval, dDM_cosmic)
    dDM_IGM = DM_interp(zvals) * fIGM * dz_vals / dzeval
    dDM_cosm = DM_interp(zvals) * dz_vals / dzeval

    # Interpolate
    sub_DM_IGM = np.cumsum(dDM_IGM)
    f_IGM = IUS(zvals, sub_DM_IGM)

    # Evaluate
    DM_IGM = f_IGM(zeval)

    # Return
    return DM_IGM * units.pc / units.cm**3, zeval

def avg_rhoISM(z):
    """
    Co-moving Mass density of the ISM

    Assumes z=0 properties for z<1
    and otherwise M_ISM = M* for z>1

    Parameters
    ----------
    z: float or ndarray
      Redshift

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
      Redshift
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
    stellar_mass_file = resource_filename('frb', os.path.join('data','IGM','stellarmass.dat'))
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
      Redshift

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
    """
    Convert input scalar or array to an array

    Parameters
    ----------
    z: float or ndarray
      Redshift

    Returns
    -------
    z: ndarray
    flg_z: int
      0 -- Input was a scalar
      1 -- Input was an array

    """
    # float or ndarray?
    if not isiterable(z):
        z = np.array([z])
        flg_z = 0
    else:
        flg_z = 1
    # Return
    return z, flg_z
