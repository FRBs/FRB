""" Module for assessing impact of intervening galaxies
   (DLAs) on FRB measurements
   Based on calclations presented in Prochaska & Neeleman 2017
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from scipy.interpolate import interp1d

from astropy import units as u

from .io import load_dla_fits


def approx_avgDM(zeval, dla_model='atan', verbose=False):
    """ Calculate the average DM from intervening galaxies
    This method is approximate (and fast) and accurate
    to better than 1% in-so-far as the analysis is correct.

    From Prochaska & Neeleman 2017

    Parameters
    ----------
    zeval : float or ndarray
      Redshift(s) for evaluation
    dla_model : str, optional

    Returns
    -------
    avgDM : Quantity (depending on type of input z)
      Units of pc/cm**3

    """
    # Init
    mlz = _model_lz(dla_model)
    if isinstance(zeval, float):
        flg_float = True
        zeval = np.array([zeval])
    else:
        flg_float = False
    # Error on upper bound
    if np.max(zeval) > 5.:
        raise IOError("Calculation is only valid to z=5")
    # Calculate
    zcalc = np.linspace(0., 5., 10000)
    dz = np.median(zcalc-np.roll(zcalc,1))
    # Load DLA fits model
    dla_fits = load_dla_fits()
    # Evaluate l(z)
    lz = mlz['eval'](zcalc)

    # Average NHI
    avgNHI = _avgN_dbl_pow(dla_fits=dla_fits)

    # Take ne/nH
    nenH_p = dla_fits['nenH']['loglog']
    nenH = nenH_p['bp'] + nenH_p['m'] * (avgNHI - 20.3)

    # Integrate lz for n(z)
    cumul = np.cumsum(lz * dz)

    # Average <z>
    avgz = np.cumsum(zcalc * lz * dz) / cumul

    '''
    # <DM> for a single DLA (rest-frame)
    DM_DLA = 10. ** (avgNHI + nenH) / u.cm ** 2
    if verbose:
        print("DM for an average DLA = {} (rest-frame)".format(DM_DLA.to('pc/cm**3')))
    '''

    # Altogether now
    avgDM_values = 10 ** avgNHI * 10 ** nenH * cumul / (1 + avgz) #/ u.cm ** 2

    # Finish up
    DM_values = np.zeros_like(zeval)
    for kk,iz in enumerate(zeval):
        iminz = np.argmin(np.abs(iz - zcalc))
        DM_values[kk] = avgDM_values[iminz]
    # Return
    return (DM_values / u.cm**2).to('pc/cm**3')


def monte_DM(zeval, model='atan', nrand=100, verbose=False):
    """
    Parameters
    ----------
    zeval : ndarray
      Array of redshifts for evaluation
    model
    nrand : int, optional
      Number of samples on NHI
    verbose : bool, optional

    Returns
    -------
    rand_DM : ndarray
      Random DM values
      Reported in cm**-2  (unitless array)

    """
    # Init
    dla_fits = load_dla_fits()
    nenH_param = dla_fits['nenH']['loglog']
    mlz = _model_lz(model)
    lgNmax = np.linspace(20.3, 22., 10000)
    intfN = _int_dbl_pow(dla_fits['fN']['dpow'], lgNmax=lgNmax)

    # Interpolate (cubic is *very* slow)
    interp_fN = interp1d(intfN/intfN[-1], lgNmax)

    # l(z) model
    # Evaluate l(z) in small z intervals
    zmax = np.max(zeval)
    z = np.linspace(0., zmax, 50000)
    dz = np.median(z-np.roll(z,1))
    lz = mlz['eval'](z, param=mlz['param'])
    # Setup for n(z) and drawing zdla
    nzc = np.cumsum(lz*dz)  # Note nzc[0] is not 0
    avgz = np.cumsum(z*lz*dz) / nzc
    interp_avgz = interp1d(z, avgz)
    nzc[0] = 0.
    interp_nz = interp1d(z, nzc)
    interp_revnz = interp1d((nzc-nzc[0])/nzc[-1], z)  # Accurate to ~1%
    #
    rand_DM = np.zeros((nrand, zeval.size))
    nz = interp_nz(zeval)
    for kk,inz in enumerate(nz):
        # Random number of DLAs
        rn = np.random.poisson(inz, size=nrand)
        ndla = np.sum(rn)
        if ndla == 0:
            continue
        # Draw NHI
        rval = np.random.uniform(size=ndla)
        rNHI = interp_fN(rval)
        # nenH
        nenH = nenH_param['bp'] + nenH_param['m'] * (rNHI-20.3)
        # Draw zdla
        rval2 = np.random.uniform(size=ndla)
        zdla = interp_revnz(rval2*inz/nzc[-1])
        # DM values
        DMi = 10.**(rNHI + nenH) / (1+zdla)
        # Generate a dummy array
        DMarr = np.zeros((nrand, max(rn)))
        cnt = 0
        for jj in range(nrand): # Fill
            if rn[jj] > 0:
                DMarr[jj,:rn[jj]] = DMi[cnt:cnt+rn[jj]]
                cnt += rn[jj]
        # Fill em up
        rand_DM[:,kk] = np.sum(DMarr,axis=1)

    # Return
    return rand_DM

def _avgN_dbl_pow(lgNmin=20.3, dla_fits=None):
    """  Calculate <NHI> for the double power-law

    Parameters
    ----------
    lgNmin : float, optional

    Returns
    -------
    avglgN : float
      log10 <NHI>

    """
    if dla_fits is None:
        dla_fits = load_dla_fits()
    # Parameters
    param = dla_fits['fN']['dpow']
    # Calculate
    fterm = 1/(param['a3']+2) - 1./(param['a4']+2)
    sterm = (10**(lgNmin-param['Nd']))**(param['a3']+2) / (param['a3']+2)
    # Numerator
    num = (10**param['Nd'])**2 *(fterm-sterm)
    # Denom
    denom = _int_dbl_pow(param, lgNmin=lgNmin)
    return np.log10(num/denom)

def _atan_lz(zval, param=None):
    """ arctan l(z) model
    Parameters
    ----------
    zval : float or ndarray

    Returns
    -------
    atan_lz : float or ndarray

    """
    if param is None:
        dfits = load_dla_fits()
        param = dfits['lz']['atan']
    lz = param['A'] + param['B'] * np.arctan(zval-param['C'])
    return lz

def _int_dbl_pow(param, lgNmin=20.3, lgNmax=None):
    """ Integrate the double power-law for f(N)
    For normalizing with l(z) and for doing random draws

    Parameters
    ----------
    lgNmin : float, optional
    lgNmax : ndarray, optional
      If None, integrate to Infinity

    Returns
    -------
    val : float or ndarray
      Integral of f(N) dN  [modulo the j(z) term]
      Really just the integral of h(N) dN

    """
    # Calculate
    if lgNmax is None:  # Integrate to Infinity
        fterm = 1/(param['a3']+1) - 1./(param['a4']+1)
    else:  # Indefinite integral
        fterm = np.zeros_like(lgNmax)
        high = lgNmax > param['Nd']
        fterm[high] = 1/(param['a3']+1) - 1./(param['a4']+1)
        fterm[high] += (10**(lgNmax[high]-param['Nd']))**(param['a4']+1) / (param['a4']+1)
        fterm[~high] = (10**(lgNmax[~high]-param['Nd']))**(param['a3']+1) / (param['a3']+1)
    # Nmin term
    sterm = (10**(lgNmin-param['Nd']))**(param['a3']+1) / (param['a3']+1)
    # Finish
    val = 10**param['Nd'] * (fterm-sterm)
    return val

def _model_lz(name):
    """ Return the model for l(z)
    Enables multiple ways to model the DLA observations
    Returns
    -------
    mlz : dict
    """
    # Fit parameters
    dfits = load_dla_fits()
    #
    mlz = dict(name=name)
    if name == 'atan':
        mlz['param'] = dfits['lz']['atan']
        mlz['eval'] = _atan_lz
    else:
        raise IOError("Bad lz model name!!")
    return mlz
