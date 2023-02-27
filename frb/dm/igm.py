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
from astropy import constants

from frb.halos import hmf as frb_hmf
from frb import mw
from frb import defs

def fukugita04_dict():
    """
    Data from Fukugita 2004, Table 1

    Returns:
      f04_dict (dict): data dict.

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

    Args:
      z (float or ndarray): redshift
      z_reion (float, optional): redshift
        of reionization.

    Returns:
      fHI (float or ndarray): float or ndarray
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

    Args:
      z (float or ndarray): Redshift
      z_HIreion (float, optional): Helium reionization
        redshift.

    Returns:
      neHe (float or ndarray): Number of free electrons
        per Helium nucelus.

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


def z_from_DM(DM, cosmo=defs.frb_cosmo, coord=None, corr_nuisance=True):
    """
    Report back an estimated redshift from an input IGM DM
    Any contributions from the Galaxy and/or host need to have been 'removed'

    Args:
      DM (Quantity): Dispersion measure.
      cosmo (Cosmology, optional): Cosmology
        of the universe. LambdaCDM with the Repo cosmology 
        used by default.
      coord (SkyCoord, optional): If provided, use it to remove the ISM
      corr_nuisance (bool, optional): If True, correct for the MW Halo
        and the host with 100 DM units.
    Returns:
        z (float): Redshift

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
    all_DM, zeval = average_DM(5., cosmo=cosmo, neval=20000, cumul=True)
    # Interpolate
    fint = interp1d(all_DM.value, zeval)
    # Evaluate
    z = fint(DM_use.to('pc/cm**3').value)
    # Return
    return z

def f_diffuse(z, cosmo=defs.frb_cosmo, return_rho = False,
              perturb_Mstar=None):
  """
  Calculate the cosmic fraction of baryons
  in diffuse gas phase based on our empirical
  knowledge of baryon distributions and their 
  ionization state.

  Args:
    z (float or ndarray): Redshift
    cosmo (Cosmology, optional): Cosmology of
      the universe.
    return_rho (bool, optional): If true, 
      the diffuse gas density
      is returned too.
    perturb_Mstar (float, optional):
      If provided, scale rho_Mstar by this value.
      Useful for exploring the uncertainty in f_diffuse

  Returns:
    f_diffuse (float, ndarray): Diffuse gas baryon fraction.
    rho_diffuse (Quantity): Physical diffuse gas density.
      Returned if return_rho is set to true.
  """
  # Get comoving baryon mass density
  rho_b = cosmo.Ob0 * cosmo.critical_density0.to('Msun/Mpc**3')

  # Dense components
  rho_Mstar = avg_rhoMstar(z, remnants=True)
  if perturb_Mstar is not None:
      rho_Mstar *= perturb_Mstar
      
  rho_ISM = avg_rhoISM(z, cosmo=cosmo)

  # Diffuse gas fraction
  f_diffuse = 1 - ((rho_Mstar+rho_ISM)/rho_b).value
  if not return_rho:
    return f_diffuse
  else:
    return f_diffuse, rho_b*f_diffuse*(1+z)**3


def ne_cosmic(z, cosmo = defs.frb_cosmo, mu = 4./3):
  """
  Calculate the average cosmic electron
  number density as a function of redshift.
  Args:
    z (float or ndarray): Redshift
    cosmo (Cosmology, optional): Cosmology in 
      which the calculations are to be performed.
    mu (float): Reduced mass
  Returns:
    ne_cosmic (Quantity): Average physical number
      density of electrons in the unverse in cm^-3.
  """
  # Get diffuse gas density
  _, rho_diffuse = f_diffuse(z, cosmo=cosmo, return_rho=True)

  # Number densities of H and He
  n_H = (rho_diffuse/constants.m_p/mu).to('cm**-3')
  n_He = n_H / 12.  # 25% He mass fraction

  # Compute electron number density
  ne_cosmic = n_H * (1.-average_fHI(z)) + n_He*(average_He_nume(z))
  return ne_cosmic

def average_DM(z, cosmo = defs.frb_cosmo, cumul=False, neval=10000, mu=4/3):
    """
    Macquart Relation

    Calculate the average cosmic DM 'expected' based on our empirical
    knowledge of baryon distributions and their ionization state.

    This includes both the IGM and galactic halos, i.e. any and all diffuse gas

    Args:
        z (float): Redshift
        mu (float): Reduced mass correction for He when calculating n_H
        cumul (bool, optional): Return the DM as a function of z

    Returns:
        DM (Quantity or Quantity array): DM values evaluated at
          the required redshifts. An array is returned only if
          cumul is True.
        zeval (ndarray): evaluation redshifts. Only returned if
          cumul is True.
    """
    # Init
    zeval, dz = np.linspace(0., z, neval,retstep=True)

    # Get n_e as a function of z
    n_e = ne_cosmic(zeval, cosmo=cosmo)

    # Cosmology -- 2nd term is the (1+z) factor for DM
    denom = cosmo.H(zeval) * (1+zeval) * (1+zeval)

    # Time to Sum
    DM_cum = (constants.c * np.cumsum(n_e * dz / denom)).to('pc/cm**3')

    # Return
    if cumul:
        return DM_cum, zeval
    else:
        return DM_cum[-1]


def average_DMhalos(z, cosmo = defs.frb_cosmo, f_hot = 0.75, rmax=1., logMmin=10.3, logMmax=16., neval = 10000, cumul=False):
    """
    Average DM_halos term from halos along the sightline to an FRB

    Args:
        z (float): Redshift of the FRB
        cosmo (Cosmology): Cosmology in which the calculations
          are to be performed.
        f_hot (float, optional): Fraction of the halo baryons in diffuse phase.
        rmax (float, optional): Size of a halo in units of r200
        logMmin (float, optional): Lowest mass halos to consider
          Cannot be much below 10.3 or the Halo code barfs
          The code deals with h^-1 factors, i.e. do not impose it yourself
        logMmax (float, optional): Highest halo mass. Default to 10^16 Msun
        neval (int, optional): Number of redshift values between
          0 and z the function is evaluated at.
        cumul (bool, optional): Return a cumulative evaluation?

    Returns:
        DM_halos (Quantity or Quantity array): One value if cumul=False
          else evaluated at a series of z
        zeval (ndarray): Evaluation redshifts if cumul=True
    """

    zeval, dz = np.linspace(0, z, neval, retstep = True)

    # Electron number density in the universe
    ne_tot = ne_cosmic(zeval, cosmo = cosmo)

    # Diffuse gas mass fraction
    f_diff = f_diffuse(zeval, cosmo = cosmo)

    # Fraction of total mass in halos
    zvals = np.linspace(0, z, 20)
    fhalos = frb_hmf.frac_in_halos(zvals, Mlow = 10**logMmin, Mhigh = 10**logMmax, rmax = rmax)
    fhalos_interp = IUS(zvals, fhalos)(zeval)

    # Electron number density in halos only
    ne_halos = ne_tot*fhalos_interp*f_hot/f_diff

    # Cosmology -- 2nd term is the (1+z) factor for DM
    denom = cosmo.H(zeval) * (1+zeval) * (1+zeval)

    # DM halos
    DM_halos = (constants.c * np.cumsum(ne_halos * dz / denom)).to('pc/cm**3')

    # Return
    if cumul:
        return DM_halos, zeval
    else:
        return DM_halos[-1]
    
def average_DMIGM(z, cosmo = defs.frb_cosmo, f_hot = 0.75, rmax=1., logMmin=10.3, neval = 10000, cumul=False):
    """
    Estimate DM_IGM in a cumulative fashion

    Args:
        z (float): Redshift of the FRB
        cosmo (Cosmology, optional): Cosmology in which 
          the calculations are to be performed. LambdaCDM
          with the Repo cosmology assumed by default.
        f_hot (float, optional): Fraction of the halo
          baryons in diffuse phase.
        rmax (float, optional):
          Size of a halo in units of r200
        logMmin (float, optional):
          Lowest mass halos to consider. Cannot be much below
          10.3 or the Halo code barfs. The code deals with
           h^-1 factors, i.e. do not impose it yourself
        neval (int, optional): Number of redshift values between
          0 and z the function is evaluated at.
        cumul (bool, optional):
          Return a cumulative evaluation?
    Returns:
        DM (Quantity or Quantity array): One value if cumul=False
          else evaluated at a series of z
        zeval (ndarray, optional): Evaluation redshifts if cumul=True
    """
    # DM cosmic
    DM_cosmic, zeval = average_DM(z, cosmo = cosmo, cumul=True, neval=neval)

    # DM_halos
    DM_halos, _ = average_DMhalos(z,cosmo = cosmo, logMmin = logMmin,
                              f_hot=f_hot, cumul = True, rmax = rmax, neval = neval)
    
    # Subtract the two
    DM_IGM = DM_cosmic - DM_halos

    # Return
    if cumul:
        return DM_IGM, zeval
    else:
        return DM_IGM[-1]

def avg_rhoISM(z, cosmo=defs.frb_cosmo):
    """
    Co-moving Mass density of the ISM

    Interpolates from z=0 values to z=1 where
    we assume M_ISM = M* and also for z>1

    Args:
        z (float or ndarray): Redshift
        cosmo (Cosmology, optional): Cosmology in which
          the calculations are to be performed. LambdaCDM
          with defs.frb_cosmo parameters assumed by default.

    Returns:
      rhoISM (Quantity): Units of Msun/Mpc^3
    """
    # Init
    z, flg_z = z_to_array(z)

    # Mstar
    rhoMstar = avg_rhoMstar(z, remnants=False)

    # z=0 (Fukugita+ 2004)
    f04_dict = fukugita04_dict()
    M_ISM = f04_dict['M_HI'] + f04_dict['M_H2']
    f_ISM_0 = M_ISM/(f04_dict['M_sphere']+f04_dict['M_disk'])

    # Assume M_ISM = M* at z=1
    f_ISM_1 = 1.

    # Ages
    t0 = cosmo.age(0.).to('Gyr').value
    t1 = cosmo.age(1.).to('Gyr').value
    t1_2 = (t0+t1)/2.
    tval = cosmo.age(z).to('Gyr').value

    # Interpolate
    f_ISM = interp1d([t0, t1_2, t1], [f_ISM_0, 0.58, f_ISM_1], kind='quadratic',
                     bounds_error=False, fill_value=1.)

    # Calculate
    rhoISM_unitless = f_ISM(tval) * rhoMstar.value

    # Finish
    rhoISM = rhoISM_unitless * units.Msun / units.Mpc**3
    #
    return rhoISM


def avg_rhoMstar(z, remnants=True):
    """
    Return mass density in stars as calculated by
    Madau & Dickinson (2014)

    Args:
      z (float or ndarray): Redshift
      remnants (bool, optional): Include remnants in the calculation?

    Returns:
      rho_Mstar (Quantity): Units of Msun/Mpc^3
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
