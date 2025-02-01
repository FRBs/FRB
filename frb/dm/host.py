""" Methods related to DM in FRB Host Galaxies """
import numpy as np

import warnings 

from astropy import units
from astropy.units import Quantity
from astropy.cosmology import Planck18

try:
    import extinction
except ImportError:
    warnings.warn("extinction package not loaded.  Extinction corrections will fail")

from frb.galaxies import nebular
from frb import em
from frb.halos import models

def dm_host_from_Halpha(z:float, Halpha:Quantity, reff_ang:Quantity,
                             AV:float=None):
    """Estimate DM_Host from Halpha and angular size (and redshift)

    Args:
        z (float): Redshift
        Halpha (Quantity): Total Halpha flux of the galaxy
        reff_ang (Quantity): Galaxy angular size
        AV (float, optional): Correct for dust if provided

    Returns:
        Quantity: DM_host as observed (i.e. includes the 1/1+z term)
    """
    # Dust correction?
    if AV is not None:
        wave = 6564. 
        al = extinction.fm07(np.atleast_1d(wave), AV)[0]
    else:
        al = 0.

    # H-alpha surface brightness
    halpha_SB = Halpha * 10**(al/2.5) / (np.pi * reff_ang**2)  

    # EM & DM
    em_R1 = em.em_from_halpha(halpha_SB, z)

    # Return
    return em.dm_from_em(em_R1, 1*units.kpc) / (1+z)

def dm_host_from_ssfr(z:float, ssfr:Quantity):
    """Estimate DM_host from the surface density of SFR 

    Args:
        z (float): Redshift
        ssfr (Quantity): Surface density of SFR (with units)

    Returns:
        Quantity: DM_host as observed (i.e. includes the 1/1+z term)

    """
    # ----------------------------
    # estimate DM_host from SSFR
    # grab SSFR
    #ssfr = Host.SSFR * units.Msun/units.yr/units.kpc**2

    halpha_kpc2 = ssfr / nebular.Ha_conversion * units.erg / units.s  # convert to Halpha surface density
    kpc_arcmin = Planck18.kpc_proper_per_arcmin(z)

    halpha_sqarcsec = halpha_kpc2 * kpc_arcmin**2   # convert to per square arcsec
    halpha_sqarcsec.to('erg/s/arcsec**2')

    D_L = Planck18.luminosity_distance(z)  # luminosity distance

    halpha_SB2 = halpha_sqarcsec / (4*np.pi*D_L**2)  # proper surface brightness
    halpha_SB2.to('erg/cm**2/s/arcsec**2')
    em_burst = em.em_from_halpha(halpha_SB2, z)

    dm_host_ssfr = em.dm_from_em(em_burst, 100*units.pc) / (1+z)  # get DM

    # Return
    return dm_host_ssfr

def dm_host_halo(R:units.Quantity, log10_Mstar:float, 
                 z:float, HMR:str='Moster',
                 mNFW:models.ModifiedNFW=None):

    # Halo mass
    if HMR == 'Moster':
        log10_Mhalo = models.halomass_from_stellarmass(log10_Mstar, z=z)
    elif HMR == 'Kravstov':
        log10_Mhalo = models.halomass_from_stellarmass_kravtsov(log10_Mstar)
    else:
        raise IOError(f"Not ready for this HMR: {HMR}")

    # Halo
    if mNFW is None:
        mNFW = models.ModifiedNFW(log_Mhalo=log10_Mhalo, 
                                  z=z, 
                                  alpha=2., 
                                  y0=2., 
                                  f_hot=0.55)


    # Find the DM
    DM_host_halo = mNFW.Ne_Rperp(R) / 2

    # Return
    return DM_host_halo
