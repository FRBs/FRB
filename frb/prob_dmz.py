""" Code for calculations of P(DM|z) and P(z|DM)"""
import numpy as np
import os

from pkg_resources import resource_filename

from scipy.stats import norm, lognorm
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from astropy.table import Table

from frb import igm

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


def grid_P_DMcosmic_z(beta=3., F=0.31):
    """
    Generate a grid of P(DM_cosmic|z)
    
    Returns:
        tuple: z, DM_cosmic, P(DM_cosmic|z)
    """
    # Check
    if not np.isclose(beta, 3.):
        raise IOError("Not prepared for this beta value (yet)")
    # Load
    # sigma_DM
    spl_sigma = grab_sigma_spline()
    f_C0_3 = grab_C0_spline()

    # Grid
    zvals = np.linspace(0., 2., 200)
    DM_cosmics = np.linspace(1., 5000., 1000)

    PDF_grid = np.zeros((DM_cosmics.size, zvals.size))

    # Loop
    for kk, zval in enumerate(zvals):
        # z=0
        if kk == 0:
            PDF_grid[0,0] = 1.
            continue
        avgDM = igm.average_DM(zval).value
        # Params
        sigma = F / np.sqrt(zval)
        C0 = f_C0_3(sigma)
        #  Delta
        Delta = DM_cosmics / avgDM
        # PDF time
        PDF = DMcosmic_PDF(Delta, C0, sigma)
        # Normalize
        PDF_grid[:,kk] = PDF / np.sum(PDF)

    # Return
    return zvals, DM_cosmics, PDF_grid


def grab_sigma_spline(redo=False):
    """

    Args:
        redo:

    Returns:

    """
    sigma_file = os.path.join(resource_filename('frb', 'data'),
                              'DM', 'sigma_sigma.ascii')
    if redo:
        raise IOError("Not ready to redo")
        npt = 200
        sigma_DMps = np.linspace(0., 0.55, npt)
        sigmas = []
        for sigma_DMp in sigma_DMps:
            res = minimize_scalar(deviate2, args=(f_C0, sigma_DMp))
            sigmas.append(float(res.x))
        # Write
        tbl = Table()
        tbl['sigma_DMp'] = sigma_DMps
        tbl['sigma'] = sigmas
        tbl.write(sigma_file, format='ascii.fixed_width', overwrite=True)
    # Load me up
    tbl = Table.read(sigma_file, format='ascii.fixed_width')
    spl_sigma = IUS(tbl['sigma_DMp'], tbl['sigma'])
    return spl_sigma

def grab_C0_spline(max_log10_sigma=0., npt=100, ret_all=False,
                    redo=False, beta=4., ifile=None):
    """
    Generate a spline of C0 vs sigma values for the
    McQuinn formalism
    """
    if redo:
        raise NotImplementedError('Not ready for this')
        sigmas = 10 ** np.linspace(-2, max_log10_sigma, npt)
        C0s = np.zeros_like(sigmas)
        # Loop
        for kk, sigma in enumerate(sigmas):
            # Minimize
            res = minimize_scalar(deviate1, args=(sigma,beta))
            C0s[kk] = res.x
    else:
        # Load from file
        if ifile is None:
            ifile = os.path.join(resource_filename('frb', 'data'),
                                 'DM', 'sigma_C0_beta3.ascii')
        tbl = Table.read(ifile, format='ascii.fixed_width')
        sigmas = tbl['sigma'].data
        C0s = tbl['C0'].data
    # Spline
    f_C0 = IUS(sigmas, C0s)
    # Return
    if ret_all:
        return f_C0, sigmas, C0s
    else:
        return f_C0
    

def DMcosmic_PDF(Delta, C0, sigma, A=1., alpha=3., beta=3.):
    """
    PDF(Delta) following the McQuinn formalism describing the DM_cosmic PDF

    Args:
        Delta (float or ndarray):
            DM / averageDM values
        C0 (float):
            parameter
        A:
        alpha:
        beta:
        sigma (float):

    Returns:

    """
    PDF = A*np.exp(-(Delta**(-alpha)-C0)**2 / (
            2*alpha**2 * sigma**2))*Delta**(-beta)
    return PDF

