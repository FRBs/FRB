""" Module for DM Halo calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
from IPython import embed

import warnings

from pkg_resources import resource_filename

from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.special import hyp2f1
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy import constants
from astropy.table import Table

# Speed up calculations
m_p = constants.m_p.cgs.value  # g

def init_hmf():
    """
    Initialize the Aemulus Halo Mass Function

    WARNING: This uses the original version which codes Tinker+2008
    We may refactor to use the more accurate, new version

    Returns:

    """
    # Hidden here to avoid it becoming a dependency
    import hmf_emulator
    # Setup HMF
    # https://github.com/astropy/astropy/blob/master/astropy/cosmology/parameters.py
    #sigma8 = 0.8159
    ns = 0.9667
    Neff = 3.046
    #cosmo_dict = {"om":cosmo.Om0,"ob":cosmo.Ob0,"ol":1.-cosmo.Om0,"ok":0.0,
    #              "h":cosmo.h,"s8":sigma8,"ns":ns,"w0":-1.0,"Neff":Neff} # "wa":0.0 is assumed internally
    cosmo_dict = {"omega_cdm":(cosmo.Om0-cosmo.Ob0)*cosmo.h**2,
                  "omega_b":cosmo.Ob0*cosmo.h**2,"ok":0.0,
                  "ln10As": 3.098, # THIS REPLACES sigma8
                  "H0":cosmo.H0.to('km/(s*Mpc)').value,
                  "n_s":ns,"w0":-1.0,"N_eff":Neff} # "wa":0.0 is assumed internally
    hmfe = hmf_emulator.hmf_emulator()
    hmfe.set_cosmology(cosmo_dict)
    # Return
    return hmfe

# Stroing for use
try:
    import hmf_emulator
except:
    pass
else:
    hmfe = init_hmf()

def frac_in_halos(zvals, Mlow, Mhigh, rmax=1.):
    """
    Calculate the fraction of matter in collapsed halos
     over a mass range and at a given redshift

    Note that the fraction of DM associated with these halos
    will be scaled down by an additional factor of f_diffuse

    Requires Aemulus HMF to be installed

    Args:
        zvals: ndarray
        Mlow: float
          In h^-1 units already so this will be applied for the halo mass function
        Mhigh: float
          In h^-1 units already
        rmax: float
          Extent of the halo in units of rvir
          WARNING: This calculation assumes a single concentration for all halos

    Returns:
        ratios: ndarray
          rho_halo / rho_m
    """

    M = np.logspace(np.log10(Mlow*cosmo.h), np.log10(Mhigh*cosmo.h), num=1000)
    lM = np.log(M)

    ratios = []
    for z in zvals:
        # Setup
        #dndlM = np.array([hmfe.dndlnM(Mi, a)[0] for Mi in M])
        dndlM = M*hmfe.dndM(M, z)
        M_spl = IUS(lM, M * dndlM)

        # Integrate
        rho_tot = M_spl.integral(np.log(Mlow*cosmo.h), np.log(Mhigh*cosmo.h)) * units.M_sun / units.Mpc ** 3
        # Cosmology
        rho_M = cosmo.critical_density(z) * cosmo.Om(z)/(1+z)**3  # Tinker calculations are all mass 
        ratio = (rho_tot*cosmo.h**2 / rho_M).decompose()
        #
        ratios.append(ratio)
    ratios = np.array(ratios)
    # Boost halos if extend beyond rvir (homologous in mass, but constant concentration is an approx)
    if rmax != 1.:
        #from pyigm.cgm.models import ModifiedNFW
        c = 7.7
        nfw = ModifiedNFW(c=c)
        M_ratio = nfw.fy_dm(rmax * nfw.c) / nfw.fy_dm(nfw.c)
        ratios *= M_ratio
    # Return
    return np.array(ratios)


def halo_incidence(Mlow, zFRB, radius=None, hmfe=None, Mhigh=1e16, nsample=20,
                   cumul=False):
    """
    Calculate the (approximate) average number of intersections to halos of a
    given minimum mass to a given zFRB.

    Requires Aemulus HMF to be installed

    Args:
        Mlow: float
          Mass of minimum halo in Solar masses
          The code deals with h^-1 factors so that you do not
          The minimum value is 2e10
        zFRB: float
          Redshift of the FRB
        radius: Quantity, optional
          The calculation will specify this radius as rvir derived from
           Mlow unless this is specified. And this rvir *will* vary with redshift
        hmfe (hmf.hmf_emulator, optional): Halo mass function emulator from Aeumulus
        Mhigh: float, optional
          Mass of maximum halo in Solar masses
        nsammple: int, optional
          Number of samplings in redshift
          20 should be enough
        cumul: bool, optional
          Return the cumulative quantities instead

    Returns:
        If cumul is False
        Navg: float
          Number of average intersections
        elif cumul is True
        zeval: ndarray
        Ncumul: ndarray
    """
    # Mlow limit
    if Mlow < 2e10:
        warnings.warn("Calculations are limited to Mlow > 2e10")
        return
    # HMF
    if hmfe is None:
        hmfe = init_hmf()
    #
    zs = np.linspace(0., zFRB, nsample)
    # Mean density
    ns = []
    for iz in zs:
        ns.append(hmfe.n_in_bins((Mlow * cosmo.h, Mhigh * cosmo.h), iz) * cosmo.h**3)  # * units.Mpc**-3
    # Interpolate
    ns = units.Quantity(ns*units.Mpc**-3)
    # Radii
    if radius is None:
        rhoc = cosmo.critical_density(zs)
        #https://arxiv.org/pdf/1312.4629.pdf eq5
        q = cosmo.Ode0/(cosmo.Ode0+cosmo.Om0*(1+zs)**3)
        rhovir = (18*np.pi**2-82*q-39*q**2)*rhoc
        r200 = (((3*Mlow*constants.M_sun.cgs) / (4*np.pi*rhovir))**(1/3)).to('kpc')
    else:
        r200 = np.ones_like(zs) * radius
    # Ap
    Ap = np.pi * r200**2

    # l(X)
    loX = ((constants.c/cosmo.H0) * ns * Ap).decompose().value

    # dX
    X = cosmo.absorption_distance(zs)
    dX = X - np.roll(X,1)
    dX[0] = 0.

    # Finish
    if cumul:
        Navg = np.cumsum(loX * dX)
        return zs, Navg
    else:
        Navg = np.sum(loX * dX)
        return Navg


def build_grid(z_FRB=1., ntrial=10, seed=12345, Mlow=1e10, r_max=2., outfile=None, dz_box=0.1,
    dz_grid=0.01, f_hot=0.75, verbose=True):
    """
    Generate a universe of dark matter halos with DM measurements
    Mainly an internal function for generating useful output grids.

    Requires the Aemulus Halo Mass function

    Args:
        z_FRB: float, optional
        ntrial: int, optional
        seed: int, optional
        Mlow: float, optional
          h^-1 mass
        r_max: float, optional
          Extent of the halo in units of rvir
        outfile: str, optional
          Write
        dz_box: float, optional
          Size of the slice of the universe for each sub-calculation
        dz_grid: float, optional
          redshift spacing in the DM grid
        f_hot: float
          Fraction of the cosmic fraction of matter in diffuse gas (for DM)

    Returns:
        DM_grid: ndarray (ntrial, nz)
        halo_tbl: Table
          Table of all the halos intersected

    """
    Mhigh = 1e16  # Msun
    # mNFW
    y0 = 2.
    alpha = 2.

    warnings.warn("Ought to do concentration properly someday!")
    cgm = ModifiedNFW(alpha=alpha, y0=y0, f_hot=f_hot)
    icm = ICM()

    # Random numbers
    rstate = np.random.RandomState(seed)

    # Init HMF
    hmfe = init_hmf()

    # Boxes
    nbox = int(z_FRB / dz_box)
    nz = int(z_FRB / dz_grid)
    dX = int(np.sqrt(ntrial))+1
    #
    npad = 6 # Mpc
    base_l = 2*dX + npad
    print('L_base = {} cMpc'.format(base_l))
    warnings.warn("Worry about being big enough given cMpc vs pMpc")

    DM_grid = np.zeros((ntrial,nz))

    # Spline distance to z
    D_max = cosmo.comoving_distance(z_FRB)
    D_val = np.linspace(1e-3,D_max.value,200) # IS THIS FINE ENOUGH?
    z_val = np.array([z_at_value(cosmo.comoving_distance, iz) for iz in D_val*units.Mpc])
    D_to_z = IUS(D_val, z_val)

    # Save halo info
    #halos = [[] for i in range(ntrial)]
    halo_i, M_i, R_i, DM_i, z_i = [], [], [], [], []

    # Loop me
    prev_zbox = 0.
    #for ss in range(nbox):
    #for ss in [0]:
    for ss in [5]:
        zbox = ss*dz_box + dz_box/2.
        print('zbox = {}'.format(zbox))
        a = 1./(1.0 + zbox) # Scale factor
        # Mass function
        M = np.logspace(np.log10(Mlow*cosmo.h), np.log10(Mhigh*cosmo.h), num=1000)
        lM = np.log(M)
        dndlM = np.array([hmf.dndlM(Mi, a) for Mi in M])
        n_spl = IUS(lM, dndlM)
        cum_n = np.array([n_spl.integral(np.log(Mlow*cosmo.h), ilM) for ilM in lM])
        ncum_n = cum_n/cum_n[-1]
        # As z increases, we have numerical issues at the high mass end (they are too rare)
        try:
            mhalo_spl = IUS(ncum_n, lM)
        except ValueError:
            # Kludge me
            print("REDUCING Mhigh by 2x")
            Mhigh /= 2.
            M = np.logspace(np.log10(Mlow*cosmo.h), np.log10(Mhigh*cosmo.h), num=1000)
            lM = np.log(M)
            dndlM = np.array([hmf.dndlM(Mi, a) for Mi in M])
            n_spl = IUS(lM, dndlM)
            cum_n = np.array([n_spl.integral(np.log(Mlow*cosmo.h), ilM) for ilM in lM])
            ncum_n = cum_n/cum_n[-1]
            #
            mhalo_spl = IUS(ncum_n, lM)

        # Volume -- Box with base l = 2Mpc
        D_zn = cosmo.comoving_distance(zbox + dz_box/2.) # Full box
        D_zp = cosmo.comoving_distance(ss*dz_box) # Previous
        D_z = D_zn - D_zp
        V = D_z * (base_l*units.Mpc)**2

        # Average N_halo
        avg_n = hmf.n_bin(Mlow*cosmo.h, Mhigh*cosmo.h, a) * cosmo.h**3 * units.Mpc**-3
        avg_N = (V * avg_n).value

        # Assume Gaussian stats for number of halos
        N_halo = int(np.round(avg_N + np.sqrt(avg_N)*rstate.randn(1)))

        # Random masses
        randM = rstate.random_sample(N_halo)
        rM = np.exp(mhalo_spl(randM)) / cosmo.h

        # r200
        r200 = (((3*rM*units.M_sun.cgs) / (4*np.pi*200*cosmo.critical_density(zbox)))**(1/3)).to('kpc')

        # Random locations (X,Y,Z)
        X_c = rstate.random_sample(N_halo)*base_l # Mpc
        Y_c = rstate.random_sample(N_halo)*base_l # Mpc
        Z_c = (rstate.random_sample(N_halo)*D_z.to('Mpc') + D_zp).value

        # Check mass fraction
        if verbose:
            Mtot = np.log10(np.sum(rM))
            M_m = (cosmo.critical_density(zbox)*cosmo.Om(zbox) * V/(1+zbox)**3).to('M_sun')
            #print("N_halo: {}  avg_N: {}".format(N_halo, avg_N))
            print("z: {}  Mhalo/M_m = {}".format(zbox, 10**Mtot/M_m.value))
            print(frac_in_halos([zbox], Mlow, Mhigh))

        # Redshifts
        z_ran = D_to_z(Z_c)

        # Loop on trials
        all_DMs = []
        all_nhalo = []
        all_r200 = []
        for itrial in range(ntrial):
            # X,Y trial
            X_trial = npad//2 + (2*itrial%dX)  # Step by 2Mpc
            Y_trial = npad//2 + 2*itrial // dX
            # Impact parameters
            try:
                R_com = np.sqrt((X_c-X_trial)**2 + (Y_c-Y_trial)**2)  # Mpc
            except:
                pdb.set_trace()
            R_phys = R_com * 1000. / (1+z_ran) * units.kpc
            # Cut
            intersect = R_phys < r_max*r200
            print("We hit {} halos".format(np.sum(intersect)))
            all_nhalo.append(np.sum(intersect))
            if not np.any(intersect):
                all_DMs.append(0.)
                continue
            # Loop -- FIND A WAY TO SPEED THIS UP!
            DMs = []
            for iobj in np.where(intersect)[0]:
                # Init
                if rM[iobj] > 1e14: # Use ICM model
                    model = icm
                else:
                    model = cgm
                model.log_Mhalo=np.log10(rM[iobj])
                model.M_halo = 10.**model.log_Mhalo * constants.M_sun.cgs
                model.z = zbox # To be consistent with above;  should be close enough
                model.setup_param(cosmo=cosmo)
                # DM
                DM = model.Ne_Rperp(R_phys[iobj], rmax=r_max, add_units=False)/(1+model.z)
                DMs.append(DM)
                # Save halo info
                halo_i.append(itrial)
                M_i.append(model.M_halo.value)
                R_i.append(R_phys[iobj].value)
                DM_i.append(DM)
                z_i.append(z_ran[iobj])
                all_r200.append(cgm.r200.value)
            # Save em
            iz = (z_ran[intersect]/dz_grid).astype(int)
            DM_grid[itrial,iz] += DMs
            all_DMs.append(np.sum(DMs))
            #print(DMs, np.log10(rM[intersect]), R_phys[intersect])
            if (itrial % 100) == 0:
                pdb.set_trace()

    # Table the halos
    halo_tbl = Table()
    halo_tbl['trial'] = halo_i
    halo_tbl['M'] = M_i
    halo_tbl['R'] = R_i
    halo_tbl['DM'] = DM_i
    halo_tbl['z'] = z_i

    # Write
    if outfile is not None:
        print("Writing to {}".format(outfile))
        np.save(outfile, DM_grid, allow_pickle=False)
        halo_tbl.write(outfile+'.fits', overwrite=True)

    return DM_grid, halo_tbl


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


def stellarmass_from_halomass(log_Mhalo,z=0):
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
    z_factor = z/(1+z)
    N = N10 + N11*z_factor
    beta = beta10 + beta11*z_factor
    gamma = gamma10 + gamma11*z_factor
    logM1 = M10 + M11*z_factor
    M1 = 10**logM1

    M_halo = 10**log_Mhalo
    
    # Simple
    log_mstar = log_Mhalo + np.log10(2*N) - np.log10((M_halo/M1)**-beta+(M_halo/M1)**gamma)
    # Done
    return log_mstar


def halomass_from_stellarmass(log_mstar,z=0):
    """ Halo mass from Stellar mass (Moster+2013).
    Inverts the function `stellarmass_from_halomass`
    numerically.

    Args:
        log_mstar (float): log_10 stellar mass
            in solar mass units.
        z (float, optional): galaxy redshift

    Returns:
        log_Mhalo (float): log_10 halo mass
            in solar mass units. 
    """
    try:
        log_mstar*z
    except ValueError:
        raise TypeError("log_mstar and z can't be broadcast together for root finding. Use numpy arrays of same length or scalar values.")

    f = lambda x: stellarmass_from_halomass(x, z = z)-log_mstar
    guess = 2+log_mstar
    return fsolve(f, guess)


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
