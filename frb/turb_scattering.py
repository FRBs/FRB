""" Codes for Scattering by ISM
  Based on formalism derived by Macquart & Koay 2013, ApJ, 776, 125
  Modified by Prochaska & Neeleman 2017
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

from scipy.special import gamma

from astropy import constants
from astropy import units

from frb import defs

from IPython import embed

# Constants
const_re = constants.alpha**2 * constants.a0


def theta_mist(n_e, nu_obs, L=50*units.kpc, R=1*units.pc, fV=1.):
    """
    Estimate the scattering angle for a mist, following the
    calculations by M. McQuinn presented in Prochaska+2019

    Args:
        n_e (Quantity):
            Electron density
        nu_obs (Quantity):
            Frequency of the radiation, observed
        L (Quantity, optional):
            Size of the region of the mist
        R (Quantity, optional):
            Size of the clouds
        fV (float, optional):
            Filling factor of the bubbles

    Returns:
        Quantity:  Angle in radians

    """

    # Constants
    me = 9.10938e-28 #(*electron mass: gr *)
    estat = 4.8027e-10 # (*charge of electron in statcoulombs *)

    phi0 = 2.2  # Numerical estimation

    theta = 2 * np.pi * estat**2 * n_e.to('cm**-3').value * phi0 / (
            me*(2*np.pi)**2 * 1e18 * nu_obs.to('GHz').value**2) * np.sqrt(2*L*fV/R) * units.radian

    # Return
    return theta


def tau_mist(n_e, nu_obs, z_FRB, zL, L=50*units.kpc, R=1*units.pc, fV=1., 
             cosmo=defs.frb_cosmo):
    """
    Temporal broadening for a mist of spherical clouds following the
    calculations by M. McQuinn presented in Prochaska+2019

    Args:
        n_e (Quantity):
            Electron density
        nu_obs (Quantity):
            Frequency of the radiation, observed
        z_FRB (float):
            Redshift of the FRB
        zL (float):
            Redshift of the intervening lens
        L (Quantity, optional):
            Size of the region of the mist
        R (Quantity, optional):
            Size of the clouds
        fV (float, optional):
            Filling factor of the bubbles
        cosmo (Cosmology, optional):

    Returns:
        Quantity: temporal broadening in seconds

    """
    D_S = cosmo.angular_diameter_distance(z_FRB)
    D_L = cosmo.angular_diameter_distance(zL)
    D_LS = cosmo.angular_diameter_distance_z1z2(zL, z_FRB)

    tau = theta_mist(n_e, nu_obs, L=L, R=R, fV=fV).to('radian').value**2 * (1+zL) * (D_LS*D_L/D_S/2/constants.c)

    return tau.to('s')


def ne_from_tau_mist(tau_scatt, z_FRB, zL, nu_obs, L=50*units.kpc, R=1*units.pc,
                     fV=1., cosmo=defs.frb_cosmo, verbose=False):
    """
    n_e from temporal broadening for a mist of spherical clouds following the
    calculations by M. McQuinn presented in Prochaska+2019

    Args:
        tau_scatt (Quantity):
            Observed width of the pulse
        z_FRB (float):
            Redshift of the FRB
        zL (float):
            Redshift of the intervening lens
        nu_obs (Quantity):
            Observed freqency
        L (Quantity, optional):
            Size of the region of the mist
        R (Quantity, optional):
            Size of the clouds
        fV (float, optional):
            Filling factor of the bubbles
        cosmo (Cosmology, optional):
        verbose (bool, optional):

    Returns:
        Quantity: density in cm**-3

    """

    D_S = cosmo.angular_diameter_distance(z_FRB)
    D_L = cosmo.angular_diameter_distance(zL)
    D_LS = cosmo.angular_diameter_distance_z1z2(zL, z_FRB)

    # Scale -- It is more correct to use the distances we used in our paper as reference (although result is very similar)
    D_term = (1.05*.262)/1.23*1000*units.Mpc

    nu_181112 = 1.3 * units.GHz
    nu_scale = (nu_obs / nu_181112)**2
    z_scale = ((1+zL)/(1+0.36738))
    cosmo_scale = ((D_L*D_LS/D_S)/D_term)

    # Branch
    Rlim = 0.011*units.pc * np.sqrt(tau_scatt/(40*units.microsecond)*cosmo_scale)/z_scale
    if R < Rlim: #0.011*units.pc * np.sqrt(tau_scatt/(40*units.microsecond)):
        if verbose:
            print("In R<{} limit".format(Rlim.to('pc')))
        n_e = (0.1*units.cm**-3) * 0.568 * np.sqrt(tau_scatt/(40*units.microsecond)) / np.sqrt(
                (L/50/units.kpc) * (fV/1e-3) * (0.1*units.pc/R))
        # Scalings
        z_scale = z_scale**(3/2)
        cosmo_scale = cosmo_scale**(-1/2)
    else:
        if verbose:
            print("In R>{} limit".format(Rlim.to('pc')))
        n_e = (0.1*units.cm**-3) * 5.03 / np.sqrt((L/50/units.kpc) * (fV/1e-3)) * (
                R/(0.1*units.pc))**(3/2)
        # Scalings
        z_scale = z_scale**2
        cosmo_scale = cosmo_scale**(-1)

    # Return
    return n_e.to('cm**-3') * cosmo_scale * nu_scale * z_scale


def ne_from_tau_kolmogorov(tau_scatt, z_FRB, zL, nu_obs, L=50*units.kpc, L0=1*units.kpc, alpha=1.,
                        cosmo=defs.frb_cosmo, debug=False):
    """
    Estimate n_e based on observed temporal broadening

    Scaled from Equation 1 of Prochaska et al. 2019

    Args:
        tau_scatt (Quantity):
            Observed width of the pulse
        z_FRB (float):
            Redshift of the FRB
        zL (float):
            Redshift of the intervening lens
        nu_obs (Quantity):
            Observed freqency
        L (Quantity, optional):
            Size of the intervening gas
        L0 (Quantity, optional):
            Turbulence length scale
        alpha (float, optional):
            Filling factor and fudge factor term
        cosmo (Cosmology, optional):

    Returns:
        Quantity: <n_e>

    """
    n_e_unscaled = 2e-3 * alpha**(-1) * (L/(50*units.kpc))**(-1/2) * (L0/(1*units.kpc))**(1/3) * (
        tau_scatt/(40*1e-6*units.s))**(5/12)
    # FRB 181112
    D_S_181112 = cosmo.angular_diameter_distance(0.47550)
    D_L_181112 = cosmo.angular_diameter_distance(0.36738)
    D_LS_181112 = cosmo.angular_diameter_distance_z1z2(0.36738, 0.47550)
    D_term = (D_L_181112*D_LS_181112/D_S_181112)**(-5/12)
    zterm = (1+0.36738)**(17/12)

    # Now scale
    D_S = cosmo.angular_diameter_distance(z_FRB)
    D_L = cosmo.angular_diameter_distance(zL)
    D_LS = cosmo.angular_diameter_distance_z1z2(zL, z_FRB)
    cosmo_scale = (D_L*D_LS/D_S)**(-5/12) / D_term

    if debug:
        print("D_S", D_S.to('Gpc'))
        print("D_L", D_L.to('Gpc'))
        print("D_LS", D_LS.to('Gpc'))
        embed(header='211 of turb_scattering')

    # Redshift
    z_scale = (1+zL)**(17/12) / zterm

    # Frequency
    nu_181112 = 1.3 * units.GHz
    nu_scale = (nu_obs / nu_181112)**(22/12)

    # Scale
    n_e = n_e_unscaled * z_scale * cosmo_scale * nu_scale

    return n_e / units.cm**3

class Turbulence(object):
    """ Class for turbulence calculations in a plasma
    Primarily used for scattering calculations
    """
    def __init__(self, ne, l0, L0, zL, beta=11./3, SM=None, verbose=True, **kwargs):
        """
        Parameters
        ----------
        ne : Quantity
          Electron density 
        l0 : Quantity
          Inner scale
        L0 : Quantity
          Outer scale
        SM : Quantity, optional
          Generally calculated but can be input
        zL : float
          Redshift of scattering medium
        beta : float, optional
          Exponent of turbulence.  Default is for Kolmogorov
        **kwargs :
          Passed to init methods
          e.g. sets SM if DL is provided
        """
        # Init
        self.beta = beta
        self.zL = zL
        self.ne = ne
        self.l0 = l0.to('pc')
        self.L0 = L0.to('kpc')
        self.verbose = verbose

        # Set SM?
        self.SM = None
        if SM is not None:
            self.SM = SM
        else:
            if 'DL' in kwargs.keys():
                self.set_SM_obj(kwargs['DL'])


        # Might check for units here (SM and beta)

        # Set rdiff based on its 'regime'
        self.regime = 0 # Undefined
        if self.SM is not None:
            if 'lobs' in kwargs.keys():
                self.set_rdiff(kwargs['lobs'])

    @property
    def CN2_gal(self):
        """ Amplitude of the turbulence per unit length
        Equation 29 from Macquarty & Koay 2013
        Assumes Kolmogorov
        
        Returns
        -------
        CN2 : Quantity
        """
        # Simple expression
        CN2 = 1.8e-3 * (self.ne/(1e-2*units.cm**(-3)))**2 * (self.L0/(0.001*units.pc))**(-2/3)
        return (CN2 * units.m**(-20/3.)).si

    @property
    def SMeff(self):
        """ Effective SM 
        Returns
        -------
        SMeff : Quantity

        """
        if self.SM is None:
            return None
        else:
            return self.SM / (1+self.zL)**2


    def set_SM_obj(self, DL):
        """ Specify SM for a discrete object (e.g. galaxy)
        Equation 31 from Macquart & Koay 2013
        Assumes Kolmogorov
        
        Parameters 
        ----------
        DL : Quantity
          Thickness of the object

        Returns
        -------

        """
        self.SM = (self.CN2_gal * DL).decompose()
        if self.verbose:
            print("Set SM={}".format(self.SM.decompose()))

    def set_cloudlet_rdiff(self, lobs, fa):
        """
        Taken from JP notes

        Args:
            lobs:
            fa (int): Number of clouds intersected

        Returns:

        """
        # ASSUMING rdiff > l0 for now
        self.rdiff = self.L0**(-1/5) * (
                2*const_re**2 * (lobs/(1+self.zL))**2 * self.ne**2 * fa)**(-3/5)

    def set_rdiff(self, lobs):
        """ Calculate rdiff in the two regimes and adopt the right one
        Requires that SM was set first

        Parameters
        ----------
        lobs : Quantity
          Observed wavelength
        Returns
        -------
        Nothing;  sets self.rdiff and self.regime
        
        """
        # Check
        if self.SM is None:
            raise IOError("Need to set SM first!")
        # Useful expression
        C = (np.pi*const_re**2 * lobs**2)/(1+self.zL)**2
        # Is rdiff < l0?
        r1 = 1. / np.sqrt(C * self.SM * (self.l0**(self.beta-4.) * (self.beta/4.) *
          gamma(-self.beta/2.)))
        # Is rdiff >> l0?
        r2 = np.power(2**(2-self.beta) * C * self.beta * self.SM * gamma(-self.beta/2.) /
          gamma(self.beta/2.), 1./(2-self.beta))
        # Query
        if r1 < self.l0:
            if self.verbose:
                print('In the regime rdiff < l_0')
            self.rdiff = r1.to('m')
            self.regime = 1  # rdiff < l0
        elif r2 > 10*self.l0:
            if self.verbose:
                print('In the regime rdiff >> l_0')
            self.rdiff = r2.to('m')
            self.regime = 2  # rdiff >> l0
        else: # Undefined
            if self.verbose:
                print('In the regime rdiff >~ l_0.  Be careful here!')
            self.rdiff = r2.to('m')
            self.regime = 2  # rdiff >> l0

    def angular_broadening(self, lobs, zsource, cosmo=defs.frb_cosmo):
        """ Broadening of a point source due to turbulent scattering

        Parameters
        ----------
        lobs : Quantity
          Observed wavelength
        zsource : float
          Redshift of radio source

        Returns
        -------
        theta : Quantity
          Angular broadening.  Radius (half-width at half-max)
        """
        if self.regime == 0:
            raise ValueError("Need to set rdiff and the regime first!")
        # f
        if (self.regime == 1) or np.isclose(self.beta,4.):
            f = 1.18
        elif (self.regime == 2) and np.isclose(self.beta, 11/3.):
            f = 1.01
        # Distances
        D_S = cosmo.angular_diameter_distance(zsource)
        D_LS = cosmo.angular_diameter_distance_z1z2(self.zL, zsource)
        if self.verbose:
            print("D_LS={}, D_S={}".format(D_LS, D_S))
        D_LS_D_S = D_LS/D_S

        # Evaluate
        k = 2*np.pi / (lobs / (1+self.zL))  # Are we sure about this (1+z) factor?!
        self.theta = f * D_LS_D_S / (k * self.rdiff) * units.radian

        return self.theta.to('arcsec')

    def temporal_smearing(self, lobs, zsource, cosmo=defs.frb_cosmo):
        """ Temporal smearing due to turbulent scattering

        Parameters
        ----------
        lobs : Quantity
          Observed wavelength
        zsource : float
        cosmo : astropy.cosmology, optional

        Returns
        -------
        tau : Quantity
          temporal broadening
        """
        D_S = cosmo.angular_diameter_distance(zsource)
        D_L = cosmo.angular_diameter_distance(self.zL)
        D_LS = cosmo.angular_diameter_distance_z1z2(self.zL, zsource)
        # Angular
        theta = self.angular_broadening(lobs, zsource, cosmo=cosmo)
        # Calculate
        tau = D_L*D_S*(theta.to('radian').value)**2 / D_LS / constants.c / (1+self.zL)
        # Return
        return tau.to('ms')


    def __repr__(self):
        txt = '<{:s}'.format(self.__class__.__name__)
        #
        txt = txt + ' ne={},'.format(self.ne.to('cm**-3'))
        txt = txt + ' l0={:.3E},'.format(self.l0.to('pc'))
        txt = txt + ' L0={},'.format(self.L0.to('pc'))
        txt = txt + ' beta={},'.format(self.beta)
        txt = txt + ' zL={}'.format(self.zL)
        #txt = txt + ' SMeff={}'.format(self.SMeff)
        txt = txt + '>'
        return (txt)




