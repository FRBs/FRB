from pkg_resources import resource_filename
import os
import numpy as np

from frb.galaxies import hosts

from IPython import embed

'''
# Globals -- to speed up calculations
r_dat, mag_uniq, _ = hosts.read_r_mags(
    resource_filename('frb', os.path.join('data', 'Galaxies', 'driver2016_t3data.fits')))
eb17_spl = hosts.interpolate.UnivariateSpline(x=mag_uniq,
                                   y=np.log10(r_dat),
                                   bbox=[-100, 100],
                                   k=3)
'''

# Spline parameters(globals) are for rmag vs sigma
driver_tck = (np.array([15., 15., 15., 15., 30., 30., 30., 30.]),
       np.array([-6.41580144, -3.53188049, -1.68500105, -0.63090954, 0., 0., 0., 0.]), 3)
driver_spl = hosts.interpolate.UnivariateSpline._from_tck(driver_tck)



def driver_sigma(rmag):
    """
<<<<<<< HEAD
    Estimated incidence of galaxies per sq arcsec with r > rmag
=======
    Estimated incidence of galaxies per sq arcsec with r > rmag 
>>>>>>> e4248983c1a7885084debb4e5f4a4281f63735b9
    using Driver et al. 2016 number counts.

    Spline parameters (globals) are for rmag vs sigma

    Args:
        rmag (float or np.ndarray):

    Returns:
        float or np.ndarray:  Galaxy number density

    """
    return 10**driver_spl(rmag)


def bloom_sigma(rmag):
    """
    Estimated incidence of galaxies per sq arcsec with r > rmag

    Args:
        rmag (float or np.ndarray):

    Returns:
        float or np.ndarray:  Galaxy density

    """
    # Sigma(m)
    sigma = 1. / (3600. ** 2 * 0.334 * np.log(10)) * 10 ** (0.334 * (rmag - 22.963) + 4.320)
    return sigma


def pchance(rmag, sep, r_half, sigR, scale_rhalf=2., nsigma=2., ndens_eval='driver'):
    """

    Args:
        rmag (np.ndarray):
            Extinction corrected apparent magnitude
        sep (np.ndarray):
            Angular separations; arcsec
        r_half (np.ndarray):
            Half light radii of the galaxies; arcsec
        sigR (float):
            1 sigma error in FRB localization; assumed symmetric; arcsec
        scale_rhalf (float, optional):
            Weighting factor for the half-light radius
        nsigma:
        ndens_eval (str, optinal):
            Number count source used
            'bloom': Hogg et al.
            'driver':  Driver et al. 2016

    Returns:
        np.ndarray:

    """

    # Reff - More conservative than usual
    Rs = np.stack([scale_rhalf * r_half, np.ones_like(r_half)* nsigma * sigR,
                   np.sqrt(sep ** 2 + (scale_rhalf * r_half) ** 2)])
    reff = np.max(Rs, axis=0)

    # Number density
    if ndens_eval =='bloom':
        nden = bloom_sigma(rmag)
    elif ndens_eval =='driver':
        nden = driver_sigma(rmag)
    else:
        raise IOError("Bad ndens evaluation")

    # Nbar
    Nbar = np.pi * reff ** 2 * nden

    # Return Pchance
    return 1. - np.exp(-Nbar)
