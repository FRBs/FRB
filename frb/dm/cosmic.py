""" Utility pieces for DM calculations"""

import os
from pkg_resources import resource_filename

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.optimize import minimize_scalar
from scipy.special import hyp1f1
from scipy.special import erf
from scipy.special import gamma

from astropy.table import Table


# Globals to speed up analysis
const1_A = np.sqrt(2) * gamma(5/6)
const2_A = 9. * gamma(4/3)
const3_A = 3 * 2**(1/6) * 3**(1/3) * np.sqrt(np.pi)

DM_values = np.linspace(1., 20000., 20000)
Delta_values = np.linspace(1./400., 20., 20000)  # Fiducial

# For analysis of the Macquart Relation
gold_frbs = ['FRB180924', 'FRB181112', 'FRB190102', 'FRB190608', 'FRB190711']

def DMcosmic_PDF(Delta, C0, sigma, A=1., alpha=3., beta=3.):
    """
    PDF(Delta) following the McQuinn formalism describing the DM_cosmic PDF

    See Macquart+2020 for details

    Args:
        Delta (float or np.ndarray):
            DM / averageDM values
        C0 (float):
            parameter
        sigma (float):
        A (float, optional):
        alpha (float, optional):
        beta (float, optional):

    Returns:
        float or np.ndarray:

    """
    PDF = A * np.exp(-(Delta ** (-alpha) - C0) ** 2 / (
            2 * alpha ** 2 * sigma ** 2)) * Delta ** (-beta)
    return PDF


def deviate1(C0, sigma, beta, orig=False):
    """
    Calculate deviate to solve fo C0

    Args:
        C0 (float):
        sigma (float):
        beta (float):
        orig (bool, optional):
            use the original approach.  Not recommended

    Returns:
        float: deviate

    """
    if orig:
        # Calculate <D>
        x = -1 * C0**2 / 18 / sigma**2
        num1 = const1_A * hyp1f1(2/3, 3/2, x)
        num2 = const2_A * sigma * hyp1f1(1/6, 1/2, x)
        cerf = erf(C0/(3*np.sqrt(2)*sigma))  # Used again below
        denom = const3_A * sigma**(4/3) * (cerf + 1)
        avgD = (num1 + num2) / denom
    else:
        PDF = DMcosmic_PDF(Delta_values, C0, sigma=sigma, beta=beta)
        avgD = np.sum(Delta_values*PDF)/np.sum(PDF)
    # Return
    return np.abs(avgD-1)


def build_C0_spline(max_log10_sigma=0., npt=100, ret_all=False, beta=4.):
    """
    Generate a spline of C0 vs sigma values for the
    McQuinn formalism

    Args:
        max_log10_sigma (float, optional):
        npt (int, optional):
        ret_all (bool, optional):
            if True, return more items
        beta (float, optional):

    Returns:
        float or tuple:  If ret_all, return f_C), sigmas, COs else return the spline

    """
    sigmas = 10 ** np.linspace(-2, max_log10_sigma, npt)
    C0s = np.zeros_like(sigmas)
    # Loop
    for kk, sigma in enumerate(sigmas):
        # Minimize
        res = minimize_scalar(deviate1, args=(sigma,beta))
        C0s[kk] = res.x
    # Spline
    f_C0 = IUS(sigmas, C0s)
    # Return
    if ret_all:
        return f_C0, sigmas, C0s
    else:
        return f_C0


def grab_sigma_spline():
    """
    Load up the sigma spline

    Args:
        redo:

    Returns:
        scipy.interpolate.InterpolatedUnivariateSpline:

    """
    sigma_file = os.path.join(resource_filename('frb', 'data'),
                              'DM', 'sigma_sigma.ascii')
    # Load me up
    tbl = Table.read(sigma_file, format='ascii.fixed_width')
    spl_sigma = IUS(tbl['sigma_DMp'], tbl['sigma'])
    return spl_sigma


def grab_C0_spline(max_log10_sigma=0., npt=100, ret_all=False,
                   redo=False, beta=3., ifile=None):
    """
    Load up the C0 spline

    Args:
        max_log10_sigma:
        npt:
        ret_all:
        redo:
        beta:
        ifile:

    Returns:
        scipy.interpolate.InterpolatedUnivariateSpline:

    """
    if redo:
        raise NotImplementedError('Not ready for this')
        sigmas = 10 ** np.linspace(-2, max_log10_sigma, npt)
        C0s = np.zeros_like(sigmas)
        # Loop
        for kk, sigma in enumerate(sigmas):
            # Minimize
            res = minimize_scalar(deviate1, args=(sigma, beta))
            C0s[kk] = res.x
    else:
        # Load from file
        if ifile is None:
            if beta != 3.:
                raise IOError("Not ready for anything other than beta=3.")
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


