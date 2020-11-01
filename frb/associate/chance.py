from pkg_resources import resource_filename
import os
import numpy as np

from frb.galaxies import hosts

from IPython import embed

# Globals -- to speed up calculations
r_dat, mag_uniq, _ = hosts.read_r_mags(
    resource_filename('frb', os.path.join('data', 'Galaxies', 'driver2016_t3data.fits')))
eb17_spl = hosts.interpolate.UnivariateSpline(x=mag_uniq,
                                   y=np.log10(r_dat),
                                   bbox=[-100, 100],
                                   k=3)
def n_gal(m_r):
    return 10 ** eb17_spl(m_r)


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


def pchance(rmag, sep, r_half, sigR, scale_rhalf=3., nsigma=3., ndens_eval='bloom'):
    """

    Args:
        rmag (np.ndarray):
        sep (np.ndarray):
        r_half (np.ndarray):
            Half light radius of the galaxy
        sigR (float):
            1 sigma error in FRB localization; assumed symmetric
        scale_rhalf:
        nsigma:

    Returns:

    """

    # Reff - More conservative than usual
    Rs = np.stack([scale_rhalf * r_half, np.ones_like(r_half)* nsigma * sigR,
                   np.sqrt(sep ** 2 + (scale_rhalf * r_half) ** 2)])
    reff = np.max(Rs, axis=0)

    # Number density
    if ndens_eval =='bloom':
        nden = bloom_sigma(rmag)
    elif ndens_eval =='eb17':
        embed(header='83 of pchance')
        # SPEED UP THE FOLLOWING!  SPLINE IT TOO
        ndens = hosts.quad(n_gal, 0, rmag)[0]
    else:
        raise IOError("Bad ndens evaluation")

    # Nbar
    Nbar = np.pi * reff ** 2 * nden

    # Return Pchance
    return 1. - np.exp(-Nbar)
