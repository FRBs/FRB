""" Code for calculations of P(DM|z) and P(z|DM)"""
import numpy as np

from scipy.stats import norm, lognorm

from IPython import embed

class P_DM_z(object):
    pass

def prob_DMcosmic_FRB(frb, DM_min=0., DM_max=5000., step=1.,
                      ISMfrac=0.10, DM_MWhalo=50.):
    """
    
    Args:
        frb (:class:frb.frb.FRB):
        DM_min (float, optional):
        DM_max (float, optional):
        step (float, optional):
            Step size of DM array in units of pc/cm**3
        ISMfrac (float, optional):
            Fraction of DM_ISM to adopt as the 1-sigma error
        DM_MWhalo (float, optional):
            Fixed value to use for the MW halo


    Returns:
        tuple:  numpy.ndarray, numpy.ndarray
            DM_cosmic values (units of pc/cm**3), P(DM_cosmic) normalized to unity

    """
    # Init
    DMcosmics = np.arange(DM_min, DM_max+step, step)
    P_DM_cosmic = np.zeros_like(DMcosmics)

    # ISM
    scale = np.pi * ISMfrac * frb.DMISM.value
    p_ISM = norm(loc=frb.DMISM.value, scale=scale)

    # Pre calculate
    DM_ISMs = DMcosmics
    pdf_ISM = p_ISM.pdf(DM_ISMs)

    # Host 
    # TODO Should use the MCMC chains to do this right!
    #  And should fix Omega_b true
    exp_u = 68.2  # Median
    sigma_host = 0.88  # Median
    lognorm_floor=0.
    p_host = lognorm(s=sigma_host, loc=lognorm_floor, scale=exp_u)

    # Loop time
    for kk, DMcosmic in enumerate(DMcosmics):
        DM_host = frb.DM.value - DM_MWhalo - DM_ISMs - DMcosmic
        # Prob time
        Prob = pdf_ISM * p_host.pdf(DM_host*(1+frb.z))
        # Sum
        P_DM_cosmic[kk] = np.sum(Prob)

    # Normalize
    P_DM_cosmic = P_DM_cosmic / np.sum(P_DM_cosmic)

    # Return
    return DMcosmics, P_DM_cosmic




    
