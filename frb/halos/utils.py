import numpy as np
from scipy.optimize import fsolve

from astropy.utils import isiterable

def stellarmass_from_halomass(log_Mhalo, z=0):
    """ Stellar mass from Halo Mass from Moster+2013
    https://doi.org/10.1093/mnras/sts261

    Args:
        log_Mhalo (float): log_10 halo mass
            in solar mass units.
        z (float, optional): halo redshift.
            Assumed to be 0 by default.
    Returns:
        log_mstar (float): log_10 galaxy stellar mass
            in solar mass units.
    """

    # Define model parameters from Table 1
    # of the paper.
    N10 = 0.0351
    N11 = -0.0247
    beta10 = 1.376
    beta11 = -0.826
    gamma10 = 0.608
    gamma11 = 0.329
    M10 = 11.59
    M11 = 1.195

    # Get redshift dependent parameters
    # from equations 11-14.
    z_factor = z / (1 + z)
    N = N10 + N11 * z_factor
    beta = beta10 + beta11 * z_factor
    gamma = gamma10 + gamma11 * z_factor
    logM1 = M10 + M11 * z_factor
    M1 = 10 ** logM1

    M_halo = 10 ** log_Mhalo

    # Simple
    log_mstar = log_Mhalo + np.log10(2 * N) - np.log10((M_halo / M1) ** -beta + (M_halo / M1) ** gamma)
    # Done
    return log_mstar


def halomass_from_stellarmass(log_mstar, z=0):
    """ Halo mass from Stellar mass (Moster+2013).
    Inverts the function `stellarmass_from_halomass`
    numerically.

    Args:
        log_mstar (float or numpy.ndarray): log_10 stellar mass
            in solar mass units.
        z (float, optional): galaxy redshift

    Returns:
        log_Mhalo (float): log_10 halo mass
            in solar mass units.
    """
    try:
        log_mstar * z
    except ValueError:
        raise TypeError(
            "log_mstar and z can't be broadcast together for root finding. Use numpy arrays of same length or scalar values.")

    f = lambda x: stellarmass_from_halomass(x, z=z) - log_mstar
    guess = 2 + log_mstar
    if isiterable(log_mstar):
        return fsolve(f, guess)
    else:
        return fsolve(f, guess)[0]