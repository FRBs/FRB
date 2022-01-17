""" Module for assessing impact of intervening galaxies
   (DLAs) on FRB measurements
   Based on calclations presented in Prochaska & Neeleman 2017
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from scipy.interpolate import interp1d

from astropy import units as u

from frb.io import load_dla_fits
from frb.turb_scattering import Turbulence
from frb import defs


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
    zeval : float or ndarray
      Array of redshifts for evaluation
    model
    nrand : int, optional
      Number of samples on NHI
    verbose : bool, optional

    Returns
    -------
    rand_DM : ndarray
      Random DM values
      Reported in pc/cm**3  (unitless array)

    """
    # Convert to array
    if isinstance(zeval, float):
        zeval = np.array([zeval])
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
    unit_conv = (1/u.cm**2).to('pc/cm**3').value
    return rand_DM * unit_conv


def monte_tau(zeval, nrand=100, nHI=0.1, avg_ne=-2.6,
              sigma_ne=0.5, cosmo=None, lobs=50*u.cm, turb=None):
    """ Generate random draws of tau at a series of redshifts

    Parameters
    ----------
    zeval : ndarray
      Array of redshifts for evaluation
    nrand : int, optional
      Number of samples on NHI
    avg_ne : float, optional
      Average log10 electron density / cm**3
    sigma_ne : float, optional
      Error in log10 ne
    nHI : float, optional
      Fiducial value for n_HI;  used for DL value
    lobs : Quantity
      Wavelength for analysis
    turb : Turbulence object, optional
      Usually defined internally and that is the highly recommended approach
    cosmo : astropy.cosmology, optional
      Defaults to defs.frb_cosmo

    Returns
    -------
    rand_tau : ndarray  (nrand, nz)
      Random tau values reported in ms (but without explicit astropy Units)
    """
    # Init
    ne_param = dict(value=avg_ne, sigma=sigma_ne)  # Neeleman+15
    dla_fits = load_dla_fits()
    if cosmo is None:
        cosmo = defs.frb_cosmo
    # Setup NHI
    lgNmax = np.linspace(20.3, 22., 10000)
    intfN = _int_dbl_pow(dla_fits['fN']['dpow'], lgNmax=lgNmax)
    # Spline
    interp_fN = interp1d(intfN/intfN[-1], lgNmax)#, kind='cubic')

    # Setup z
    zvals = np.linspace(0., 7., 10000)
    nz_s = _dla_nz(zvals)
    nz_s[0] = 0.

    # Turbulence
    if turb is None:
        turb = _init_dla_turb()
    f_ne=turb.ne
    zsource = 2.
    turb.set_rdiff(lobs)
    fiducial_tau = turb.temporal_smearing(lobs, zsource)
    # Take out the cosmology
    f_D_S = cosmo.angular_diameter_distance(zsource)
    f_D_L = cosmo.angular_diameter_distance(turb.zL)
    f_D_LS = cosmo.angular_diameter_distance_z1z2(turb.zL, zsource)
    fiducial_tau = fiducial_tau / f_D_LS / f_D_L * f_D_S * (1+turb.zL)**3  # ms/Mpc
    kpc_cm = (1*u.kpc).to('cm').value

    rand_tau = np.zeros((nrand, zeval.size))
    # Loop on zeval
    for ss,izeval in enumerate(zeval):
        avg_nz = _dla_nz(izeval)
        rn = np.random.poisson(avg_nz, size=nrand)
        ndla = np.sum(rn)
        if ndla == 0:
            continue
        # Get random NHI
        rval = np.random.uniform(size=ndla)
        rNHI = interp_fN(rval)
        DL = 10.**rNHI / nHI / kpc_cm
        # Get random z
        imin = np.argmin(np.abs(zvals-izeval))
        interp_z = interp1d(nz_s[0:imin]/nz_s[imin-1], zvals[0:imin])#, kind='cubic')
        rval = np.random.uniform(size=ndla)
        rz = interp_z(rval)
        # Cosmology
        D_S = cosmo.angular_diameter_distance(izeval)
        D_L = cosmo.angular_diameter_distance(rz)
        D_LS = cosmo.angular_diameter_distance_z1z2(rz, izeval)
        # Get random n_e
        rne = 10.**(ne_param['value'] + ne_param['sigma']*np.random.normal(size=ndla))
        # Calculate (scale)
        rtau = fiducial_tau * (D_LS * D_L / D_S) * (rne/f_ne.to('cm**-3').value)**2 / (1+rz)**3
        # Generate, fill
        taus = np.zeros((nrand, np.max(rn)))
        kk = 0
        for jj,irn in enumerate(rn):
            if irn > 0:
                taus[jj,0:irn] = rtau[kk:kk+irn]
                kk += irn
        # Finish -- add in quadrature
        final_tau = np.sqrt(np.sum(taus**2, axis=1))
        # Save
        rand_tau[:,ss] = final_tau

    # Return
    return rand_tau


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


def _dla_nz(zarr, mlz=None, model='atan'):
    """ Calculate the number of DLAs intersected on average
    to a given redshift
    Parameters
    ----------
    zarr : ndarray
    mlz : model, optional
    model : str, optional

    Returns
    -------
    nz : ndarray

    """
    # Load model
    if mlz is None:
        mlz = _model_lz(model)

    z = np.linspace(0., 10., 10000)
    dz = np.median(z-np.roll(z,1))
    lz = mlz['eval'](z, param=mlz['param'])

    # Sum
    nz = np.cumsum(lz*dz)

    # Interpolate onto input redshifts
    interp_nz = interp1d(z, nz)

    # Return
    return interp_nz(zarr)


def _init_dla_turb(ne=4e-3/u.cm**3, zL=1.):
    """ Initialize a Turbulence object for a fiducial DLA
    Parameters
    ----------
    ne : Quantity
      Electron density
      Default is based on Neeleman+15
    zL : float
      Redshift of the DLA

    Returns
    -------
    turb : Turbulence object
    """
    # Sizes
    l0 = 1 * u.AU
    L0 = 0.001 * u.pc
    DL = 1 * u.kpc
    # Init
    turb = Turbulence(ne, l0, L0, zL)
    turb.set_SM_obj(DL)
    # Return
    return turb

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
