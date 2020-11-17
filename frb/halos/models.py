""" Module for DM Halo calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

import warnings

from pkg_resources import resource_filename

from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.special import hyp2f1
from scipy.interpolate import interp1d

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.cosmology import Planck15 as cosmo
from astropy import constants
from astropy.cosmology import z_at_value
from astropy.table import Table

from IPython import embed

# Speed up calculations
m_p = constants.m_p.cgs.value  # g


class ModifiedNFW(object):
    """ Generate a modified NFW model, e.g. Mathews & Prochaska 2017
    for the hot, virialized gas.

    Parameters:
        log_Mhalo: float, optional
          log10 of the Halo mass (solar masses)
        c: float, optional
          concentration of the halo
        f_hot: float, optional
          Fraction of the baryons in this hot phase
          Will likely use this for all diffuse gas
        alpha: float, optional
          Parameter to modify NFW profile power-law
        y0: float, optional
          Parameter to modify NFW profile position.
        z: float, optional
          Redshift of the halo
        cosmo: astropy cosmology, optional
          Cosmology of the universe. Planck15 by default.

    Attributes:
        H0: Quantity;  Hubble constant
        fb: float; Cosmic fraction of baryons (stars+dust+gas) in the entire halo
           Default to 0.16
        r200: Quantity
           Virial radius
        rho0: Quantity
           Density normalization
        M_b: Quantity
           Mass in baryons of the


    """
    def __init__(self, log_Mhalo=12.2, c=7.67, f_hot=0.75, alpha=0.,
                 y0=1., z=0., cosmo=cosmo, **kwargs):
        # Init
        # Param
        self.log_Mhalo = log_Mhalo
        self.M_halo = 10.**self.log_Mhalo * constants.M_sun.cgs
        self.c = c
        self.alpha = alpha
        self.y0 = y0
        self.z = z
        self.f_hot = f_hot
        self.zero_inner_ne = 0. # kpc
        self.cosmo = cosmo

        # Init more
        self.setup_param(cosmo=self.cosmo)

    def setup_param(self,cosmo):
        """ Setup key parameters of the model
        """
        # Cosmology
        if cosmo is None:
            self.rhoc = 9.2e-30 * units.g / units.cm**3
            self.fb = 0.16       # Baryon fraction
            self.H0 = 70. *units.km/units.s/ units.Mpc
        else:
            self.rhoc = self.cosmo.critical_density(self.z)
            self.fb = cosmo.Ob0/cosmo.Om0
            self.H0 = cosmo.H0
        # Dark Matter
        self.q = self.cosmo.Ode0/(self.cosmo.Ode0+self.cosmo.Om0*(1+self.z)**3) 
        #r200 = (((3*Mlow*constants.M_sun.cgs) / (4*np.pi*200*rhoc))**(1/3)).to('kpc')
        self.rhovir = (18*np.pi**2-82*self.q-39*self.q**2)*self.rhoc
        self.r200 = (((3*self.M_halo) / (4*np.pi*self.rhovir))**(1/3)).to('kpc')
        self.rho0 = self.rhovir/3 * self.c**3 / self.fy_dm(self.c)   # Central density
        # Baryons
        self.M_b = self.M_halo * self.fb
        self.rho0_b = (self.M_b / (4*np.pi) * (self.c/self.r200)**3 / self.fy_b(self.c)).cgs
        # Misc
        self.mu = 1.33   # Reduced mass correction for Helium

    def fy_dm(self, y):
        """ Enclosed mass function for the Dark Matter NFW
        Assumes the NFW profile

        Parameters
        ----------
        y : float or ndarray
          y = c(r/r200)

        Returns
        -------
        f_y : float or ndarray
        """
        f_y = np.log(1+y) - y/(1+y)
        #
        return f_y

    def fy_b(self, y):
        """ Enclosed mass function for the baryons

        Parameters
            y: float or ndarray

        Returns
        -------
            f_y: float or ndarray
              Enclosed mass
        """
        f_y = (y/(self.y0 + y))**(1+self.alpha) * (
                self.y0**(-self.alpha) * (self.y0 + y)**(1+self.alpha) * hyp2f1(
            1+self.alpha, 1+self.alpha, 2+self.alpha, -1*y/self.y0)
                - self.y0) / (1+self.alpha) / self.y0
        return f_y

    def ne(self, xyz):
        """ Calculate n_e from n_H with a correction for Helium
        Assume 25% mass is Helium and both electrons have been stripped

        Parameters
        ----------
        xyz : ndarray (3, npoints)
          Coordinate(s) in kpc

        Returns
        -------
        n_e : float or ndarray
          electron density in cm**-3

        """
        ne = self.nH(xyz) * 1.1667
        if self.zero_inner_ne > 0.:
            rad = np.sum(xyz**2, axis=0)
            inner = rad < self.zero_inner_ne**2
            if np.any(inner):
                if len(xyz.shape) == 1:
                    ne = 0.
                else:
                    ne[inner] = 0.
        # Return
        return ne

    def nH(self, xyz):
        """ Calculate the Hydrogen number density
        Includes a correction for Helium

        Parameters
        ----------
        xyz : ndarray
          Coordinate(s) in kpc

        Returns
        -------
        nH : float or ndarray
          Density in cm**-3

        """
        nH = (self.rho_b(xyz) / self.mu / m_p).cgs.value
        # Return
        return nH

    def rho_b(self, xyz):
        """ Mass density in baryons in the halo; modified

        Parameters
        ----------
        xyz : ndarray
          Position (assumes kpc)

        Returns
        -------
        rho : Quantity
          Density in g / cm**-3

        """
        radius = np.sqrt(rad3d2(xyz))
        y = self.c * (radius/self.r200.to('kpc').value)
        rho = self.rho0_b * self.f_hot / y**(1-self.alpha) / (self.y0+y)**(2+self.alpha)
        # Return
        return rho

    def Ne_Rperp(self, Rperp, step_size=0.1*units.kpc, rmax=1., add_units=True, cumul=False):
        """ Calculate N_e at an input impact parameter Rperp
        Just a simple sum in steps of step_size

        Parameters
        ----------
        Rperp : Quantity
          Impact parameter, typically in kpc
        step_size : Quantity, optional
          Step size used for numerical integration (sum)
        rmax : float, optional
          Maximum radius for integration in units of r200
        add_units : bool, optional
          Speed up calculations by avoiding units
        cumul: bool, optional

        Returns
        -------
        if cumul:
          zval: ndarray (kpc)
             z-values where z=0 is the midplane
          Ne_cumul: ndarray
             Cumulative Ne values (pc cm**-3)
        else:
          Ne: Quantity
             Column density of total electrons
        """
        dz = step_size.to('kpc').value

        # Cut at rmax*rvir
        if Rperp > rmax*self.r200:
            if add_units:
                return 0. / units.cm**2
            else:
                return 0.
        # Generate a sightline to rvir
        zmax = np.sqrt((rmax*self.r200) ** 2 - Rperp ** 2).to('kpc')
        zval = np.arange(-zmax.value, zmax.value+dz, dz)  # kpc
        # Set xyz
        xyz = np.zeros((3,zval.size))
        xyz[0, :] = Rperp.to('kpc').value
        xyz[2, :] = zval

        # Integrate
        ne = self.ne(xyz) # cm**-3
        if cumul:
            Ne_cumul = np.cumsum(ne) * dz * 1000  # pc cm**-3
            return zval, Ne_cumul
        Ne = np.sum(ne) * dz * 1000  # pc cm**-3

        # Return
        if add_units:
            return Ne * units.pc / units.cm**3
        else:
            return Ne

    def RM_Rperp(self, Rperp, Bparallel, step_size=0.1*units.kpc, rmax=1.,
                 add_units=True, cumul=False, zmax=None):
        """ Calculate RM at an input impact parameter Rperp
        Just a simple sum in steps of step_size
        Assumes a constant Magnetic field

        Parameters
        ----------
        Rperp : Quantity
          Impact parameter, typically in kpc
        Bparallel (Quantity):
          Magnetic field
        step_size : Quantity, optional
          Step size used for numerical integration (sum)
        rmax : float, optional
          Maximum radius for integration in units of r200
        add_units : bool, optional
          Speed up calculations by avoiding units
        cumul: bool, optional
        zmax: float, optional
          Maximum distance along the sightline to integrate.
          Default is rmax*rvir

        Returns
        -------
        if cumul:
          zval: ndarray (kpc)
             z-values where z=0 is the midplane
          Ne_cumul: ndarray
             Cumulative Ne values (pc cm**-3)
        else:
          RM: Quantity
             Column density of total electrons
        """
        dz = step_size.to('kpc').value

        # Cut at rmax*rvir
        if Rperp > rmax*self.r200:
            if add_units:
                return 0. / units.cm**2
            else:
                return 0.
        # Generate a sightline to rvir
        if zmax is None:
            zmax = np.sqrt((rmax*self.r200) ** 2 - Rperp ** 2).to('kpc')
        zval = np.arange(-zmax.value, zmax.value+dz, dz)  # kpc
        # Set xyz
        xyz = np.zeros((3,zval.size))
        xyz[0, :] = Rperp.to('kpc').value
        xyz[2, :] = zval

        # Integrate
        ne = self.ne(xyz) # cm**-3
        # Using Akahori & Ryu 2011
        RM = 8.12e5 * Bparallel.to('microGauss').value * \
             np.sum(ne) * dz / 1000  # rad m**-2

        if cumul:
            RM_cumul = 8.12e5 * Bparallel.to('microGauss') * np.cumsum(
                ne) * dz / 1000  # rad m**-2
            return zval, RM_cumul

        # Return
        if add_units:
            return RM * units.rad / units.m**2
        else:
            return RM

    def mass_r(self, r, step_size=0.1*units.kpc):
        """ Calculate baryonic halo mass (not total) to a given radius
        Just a simple sum in steps of step_size

        Parameters
        ----------
        r : Quantity
          Radius, typically in kpc
        step_size : Quantity, optional
          Step size used for numerical integration (sum)

        Returns
        -------
          Mr: Quantity
             Enclosed baryonic mass within r
             Msun units
        """
        dr = step_size.to('kpc').value

        # Generate a sightline to rvir
        rval = np.arange(0., r.to('kpc').value+dr, dr)  # kpc

        # Set xyz
        xyz = np.zeros((3,rval.size))
        xyz[2, :] = rval

        # Integrate
        nH = self.nH(xyz)  # cm**-3
        Mr_number = 4*np.pi*np.sum(nH*rval**2) * dr * self.mu * m_p  # g kpc**3/cm**3
        Mr = Mr_number * units.g * (units.kpc**3)/(units.cm**3)#

        # Return
        return Mr.to('M_sun')

    def __repr__(self):
        txt = '<{:s}: {:s} {:s}, logM={:f}, r200={:g}'.format(
                self.__class__.__name__,
                self.coord.icrs.ra.to_string(unit=units.hour,sep=':',pad=True),
                self.coord.icrs.dec.to_string(sep=':',pad=True,alwayssign=True),
                np.log10(self.M_halo.to('Msun').value),
            self.r200)
        # Finish
        txt = txt + '>'
        return (txt)



class MB04(ModifiedNFW):
    """
    Halo based on the Maller & Bullock (2004) model of
    virialized halo gas.

    Parameters:
        Rc: Quantity
          cooling radius

    """
    def __init__(self, Rc=167*units.kpc, log_Mhalo=12.2, c=7.67, f_hot=0.75, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot, **kwargs)

        # Setup
        self.Rs = self.r200/self.c
        self.Rc = Rc
        self.Cc = (self.Rc/self.Rs).decompose().value
        self.rhoV = 1. * constants.m_p/units.cm**3  # Will be renormalized

        # For development
        self.debug=False

        # Normalize
        self.norm_rhoV()

    def norm_rhoV(self):
        """
        Normalize the density constant from MB04

        Returns:

        """
        # Set rhoV to match expected baryon mass
        r = np.linspace(1., self.r200.to('kpc').value, 1000)  # kpc
        # Set xyz
        xyz = np.zeros((3,r.size))
        xyz[2, :] = r
        #
        dr = r[1] - r[0]
        Mass_unnorm = 4 * np.pi * np.sum(r**2 * self.rho_b(xyz)) * dr * units.kpc**3 # g * kpc**3 / cm**3
        # Ratio
        rtio = (Mass_unnorm/self.M_b).decompose().value
        self.rhoV = self.rhoV.cgs/rtio
        #
        print("rhoV normalized to {} to give M_b={}".format((self.rhoV/constants.m_p).cgs,
                                                            self.M_b.to('Msun')))

    def rho_b(self, xyz):
        """
        Baryonic density profile

        Args:
            xyz: ndarray
              Position array assumed in kpc

        Returns:

        """
        radius = np.sqrt(rad3d2(xyz))
        x = radius/self.Rs.to('kpc').value
        #
        rho = self.rhoV * (1+ (3.7/x)*np.log(1+x) - (3.7/self.Cc) * np.log(1+self.Cc))**(3/2)
        if self.debug:
            pdb.set_trace()
        #
        return rho


class YF17(ModifiedNFW):
    """
    Y. Faerman et al (2017) model of the Milky Way

    For the un-normalized density profile, we adopt the
    average of the warm and hot components in
    """
    def __init__(self, log_Mhalo=12.18, c=7.67, f_hot=0.75, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot, **kwargs)

        # Read
        #faerman_file = resource_filename('pyigm', '/data/CGM/Models/Faerman_2017_ApJ_835_52-density-full.txt')
        faerman_file = resource_filename('frb', '/data/Halos/Faerman_2017_ApJ_835_52-density-full.txt')
        self.yf17 = Table.read(faerman_file, format='ascii.cds')
        self.yf17['nH'] = self.yf17['nHhot'] + self.yf17['nHwarm']

        # For development
        self.debug=False

        # Setup
        self.rhoN = constants.m_p/units.cm**3
        self.setup_yfdensity()

    def setup_yfdensity(self):
        """
        Normalize the density profile from the input mass

        Returns:
            Initializes self.rhoN, the density normalization

        """
        # Setup Interpolation
        self.yf17_interp = interp1d(self.yf17['Radius'], self.yf17['nH'], kind='cubic', bounds_error=False, fill_value=0.)

        # Set rhoN to match expected baryon mass
        r = np.linspace(1., self.r200.to('kpc').value, 1000)  # kpc
        # Set xyz
        xyz = np.zeros((3,r.size))
        xyz[2, :] = r
        #
        dr = r[1] - r[0]
        Mass_unnorm = 4 * np.pi * np.sum(r**2 * self.rho_b(xyz)) * dr * units.kpc**3 # g * kpc**3 / cm**3
        # Ratio
        rtio = (Mass_unnorm/self.M_b).decompose().value
        self.rhoN = self.rhoN.cgs/rtio
        #
        print("rhoN normalized to {} to give M_b={}".format((self.rhoN/constants.m_p).cgs,
                                                            self.M_b.to('Msun')))

    def rho_b(self, xyz):
        """
        Calculate the baryonic density

        Args:
            xyz: ndarray
              Coordinates in kpc

        Returns:
            rho: Quantity array
              Baryonic mass density (g/cm**3)

        """
        radius = np.sqrt(rad3d2(xyz))
        #
        rho = self.rhoN * self.yf17_interp(radius)
        if self.debug:
            pdb.set_trace()
        #
        return rho


class MB15(ModifiedNFW):
    """
    Encodes the Galactic halo profile from
    Miller & Bregman 2015, ApJ, 800, 14
    https://ui.adsabs.harvard.edu/abs/2015ApJ...800...14M/abstract

    The default normalization and beta values are taken from their Table 2, last row.
    The models presented there do not appear to vary too much.

    """
    def __init__(self, log_Mhalo=12.18, c=7.67, f_hot=0.75, **kwargs):
        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot, **kwargs)

        # Best parameters
        self.beta = 0.45
        self.n0_rc3b = 0.79e-2  # Last entry of Table 2; Crazy units

    def nH(self, xyz):
        """
        Calculate the number density of Hydrogen

        Args:
            xyz: ndarray
              Coordinates in kpc

        Returns:
            ndarray: Number density with units of 1/cm**3

        """
        radius = np.sqrt(rad3d2(xyz))
        #  Equation 2 of Miller & Bregman 2015
        nH = self.n0_rc3b / radius**(3*self.beta)
        #
        return nH # / units.cm**3

class MilkyWay(ModifiedNFW):
    """ Fiducial model for the Galaxy

    Halo mass follows latest constraints

    Density profile is similar to Maller & Bullock 2004

    """
    def __init__(self, log_Mhalo=12.18, c=7.67, f_hot=0.75, alpha=2, y0=2, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot,
                             alpha=alpha, y0=y0, **kwargs)


class M31(ModifiedNFW):
    """
    Preferred model for M31

    Taking mass from van der Marel 2012

    """
    def __init__(self, log_Mhalo=12.18, c=7.67, f_hot=0.75, alpha=2, y0=2, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot,
                             alpha=alpha, y0=y0, **kwargs)
        # Position from Sun
        self.distance = 752 * units.kpc # (Riess, A.G., Fliri, J., & Valls - Gabaud, D. 2012, ApJ, 745, 156)
        self.coord = SkyCoord('J004244.3+411609', unit=(units.hourangle, units.deg),
                              distance=self.distance)

    def DM_from_Galactic(self, scoord, **kwargs):
        """
        Calculate DM through M31's halo from the Sun
        given a direction

        Args:
            scoord:  SkyCoord
               Coordinates of the sightline
            **kwargs:
               Passed to Ne_Rperp

        Returns:
            DM: Quantity
              Dispersion measure through M31's halo
        """
        # Setup the geometry
        a=1
        c=0
        x0, y0 = self.distance.to('kpc').value, 0. # kpc
        # Seperation
        sep = self.coord.separation(scoord)
        # More geometry
        atan = np.arctan(sep.radian)
        b = -1 * a / atan
        # Restrct to within 90deg (everything beyond is 0 anyhow)
        if sep > 90.*units.deg:
            return 0 * units.pc / units.cm**3
        # Rperp
        Rperp = np.abs(a*x0 + b*y0 + c) / np.sqrt(a**2 + b**2)  # kpc
        # DM
        DM = self.Ne_Rperp(Rperp*units.kpc, **kwargs).to('pc/cm**3')
        return DM


class LMC(ModifiedNFW):
    """
    Preferred model for LMC

    Taking data from D'Onghia & Fox ARAA 2016

    """
    def __init__(self, log_Mhalo=np.log10(1.7e10), c=12.1, f_hot=0.75, alpha=2, y0=2, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot,
                             alpha=alpha, y0=y0, **kwargs)
        # Position from Sun
        self.distance = 50 * units.kpc
        self.coord = SkyCoord('J052334.6-694522', unit=(units.hourangle, units.deg),
                              distance=self.distance)

class SMC(ModifiedNFW):
    """
    Preferred model for SMC

    Taking data from D'Onghia & Fox ARAA 2016

    """
    def __init__(self, log_Mhalo=np.log10(2.4e9), c=15.0, f_hot=0.75, alpha=2, y0=2, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot,
                             alpha=alpha, y0=y0, **kwargs)
        # Position from Sun
        self.distance = 61 * units.kpc
        self.coord = SkyCoord('J005238.0-724801', unit=(units.hourangle, units.deg),
                              distance=self.distance)

class M33(ModifiedNFW):
    """
    Preferred model for SMC

    Taking data from Corbelli 2006

    """
    def __init__(self, log_Mhalo=np.log10(5e11), c=8.36, f_hot=0.75, alpha=2, y0=2, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot,
                             alpha=alpha, y0=y0, **kwargs)
        # Position from Sun
        self.distance = 840 * units.kpc
        self.coord = SkyCoord(ra=23.4621*units.deg, dec=30.6600*units.deg, distance=self.distance)


class ICM(ModifiedNFW):
    """
    Intracluster medium (ICM) model following the analysis
    of Vikhilnin et al. 2006

    We scale the model to the profile fitted to A907

    """
    def __init__(self, log_Mhalo=np.log10(5e14), c=5, f_hot=0.70, **kwargs):
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot, **kwargs)

    def setup_param(self, cosmo=None):
        super(ICM, self).setup_param(cosmo=cosmo)
        # Scale the profile by r200
        self.scale_profile()

    def scale_profile(self):
        # Using the Vihilnin et al. 2006 values for A907
        self.a907_r200 = 1820 * units.kpc  # Derived in the method below and hard-coded here
        self.a907_c200 = 5.28
        # A907 values
        self.a907_n0 = 6.252e-3 #/ u.cm**3
        self.a907_rc = 136.9 * (self.r200/self.a907_r200).decompose() #* u.kpc
        self.a907_rs = 1887.1 * (self.r200/self.a907_r200).decompose() #* u.kpc
        self.a907_alpha = 1.556
        self.a907_beta = 0.594
        self.a907_epsilon = 4.998
        self.a907_n02 = 0.

        # Scale/set
        self.rc = self.a907_rc * (self.r200/self.a907_r200).decompose() #* u.kpc
        self.rs = self.a907_rs * (self.r200/self.a907_r200).decompose() #* u.kpc
        self.alpha = self.a907_alpha
        self.beta = self.a907_beta
        self.epsilon = self.a907_epsilon
        self.n02 = self.a907_n02
        self.n0 = 6.252e-3 #/ u.cm**3  (temporary)

        # Fixed
        self.gamma = 3

        # Now the hot gas mass for the central density
        Mb_M200 = self.mass_r(self.r200)
        self.n0 *= (self.M_b*self.f_hot/Mb_M200).decompose()

    def a907_nfw(self):
        """
        Code to regenrate the r200 and c200 values for A907
        Now hard-coded
        """
        self.a907_c500 = 3.5
        self.a907_M500 = 5e14 * units.Msun
        self.a907_r500 = (((3*self.a907_M500) / (4*np.pi*500*self.rhoc))**(1/3)).to('kpc')
        self.a907_Rs = self.a907_r500 / self.a907_c500  # Do not confuse with rs
        # Code to re-calculate these
        fy_500 = self.fy_dm(self.a907_r500 / self.a907_Rs)
        yval = np.linspace(3.5, 10, 100)
        rval = self.a907_Rs * yval
        Mval = self.a907_M500 * self.fy_dm(yval) / fy_500
        avg_rho = Mval / (4 * np.pi * rval ** 3 / 3.)
        scaled_rho = (avg_rho / (200 * self.rhoc)).decompose()
        srt = np.argsort(scaled_rho)
        f_Mr = IUS(scaled_rho[srt], rval[srt])
        self.a907_r200 = float(f_Mr(1.))*units.kpc
        self.a907_c200 = (self.a907_r200 / self.a907_Rs).decompose()
        self.a907_M200 = self.a907_M500 * self.fy_dm(self.a907_r200/self.a907_Rs) / fy_500

    def ne(self, xyz):
        """

        Parameters
        ----------
        xyz : ndarray
          Coordinate(s) in kpc

        Returns
        -------
        n_e : float or ndarray
          electron density in cm**-3

        """
        radius = np.sqrt(rad3d2(xyz))

        npne = np.zeros_like(radius)

        # Zero out inner 10kpc
        ok_r = radius > 10.

        # This ignores the n02 term
        npne[ok_r] = self.n0**2 * (radius[ok_r]/self.rc)**(-self.alpha) / (
                (1+(radius[ok_r]/self.rc)**2)**(3*self.beta - self.alpha/2.)) * (1 /
                                                                           (1+(radius[ok_r]/self.rs)**self.gamma)**(self.epsilon/self.gamma))
        if self.n02 > 0:
            pdb.set_trace()  # Not coded yet

        ne = np.sqrt(npne * 1.1667)
        # Return
        return ne

    def nH(self, xyz):
        """
        Scale by He

        Args:
            xyz:

        Returns:

        """
        return self.ne(xyz) / 1.1667


class Virgo(ICM):
    """
    Parameterization of Virgo following the Planck Collaboration
    paper:  A&A 596 A101 (2016)
    """
    def __init__(self, log_Mhalo=np.log10(1.2e14*(cosmo.Om0/cosmo.Ob0)), **kwargs):
        ICM.__init__(self, log_Mhalo=log_Mhalo, **kwargs)

        # Position from Sun
        self.distance = 18 * units.Mpc
        self.coord = SkyCoord('J123049+122328',  # Using M87
                              unit=(units.hourangle, units.deg),
                              distance=self.distance)

    def setup_param(self, cosmo=None):
        """ Setup key parameters of the model
        """
        self.r200 = 1.2 * units.Mpc

    def ne(self, xyz):
        radius = np.sqrt(rad3d2(xyz))

        # Equation 8
        ne = 8.5e-5 / (radius/1e3)**1.2

        # Return
        return ne


def rad3d2(xyz):
    """ Calculate radius to x,y,z inputted
    Assumes the origin is 0,0,0

    Parameters
    ----------
        xyz : Tuple or ndarray

    Returns
    -------
        rad3d : float or ndarray

    """
    return xyz[0]**2 + xyz[1]**2 + xyz[-1]**2

