""" Module for assessing impact of intervening galaxies
   (DLAs) on FRB measurements
   Based on calclations presented in Prochaska & Neeleman 2017
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from astropy import units as u

from .io import load_dla_fits

def avgDM(zeval, use_boot=False, verbose=False):
    """ Calculate the average DM from intervening galaxies

    Parameters
    ----------
    zeval : float or ndarray
      Redshift(s) for evaluation
    use_boot

    Returns
    -------
    avgDM : Quantity (depending on type of input z)
      Units of pc/cm**3

    """
    # Init
    if isinstance(zeval, float):
        flg_float = True
        zeval = np.array([zeval])
    else:
        flg_float = False
    #
    if np.max(zeval) > 5.:
        raise IOError("Calculation is only valid to z=5")
    #
    if use_boot:
        pdb.set_trace()
    else:
        zcalc = np.linspace(0., 5., 1000)
        dz = np.median(zcalc-np.roll(zcalc,1))
        # Load DLA fits model
        dla_fits = load_dla_fits()
        # Evaluate l(z)
        param = dla_fits['lz']['atan']
        lz = param['A'] + param['B'] * np.arctan(zcalc-param['C'])

        # Average NHI
        avgNHI = avgN_dbl_pow(dla_fits=dla_fits)

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


def avgN_dbl_pow(lgNmin=20.3, dla_fits=None):
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
    denom = int_dbl_pow(param, lgNmin=lgNmin)
    return np.log10(num/denom)

def int_dbl_pow(param, lgNmin=20.3, lgNmax=None):
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