"""
Generalized Halo Profile Model: g × M_DM Framework

This module implements the baryon density profile framework derived in 
generalized_halo_profile.ipynb. The key insight is to define the enclosed 
baryon mass as:

    M_B(<y) = g(y) × (Ω_B / (Ω_M - Ω_B)) × M_DM(<y)

where g(y) → 1 as y → ∞, ensuring the enclosed baryon fraction approaches 
the cosmic mean at large radii.

The baryon density profile is then derived via the fundamental theorem of calculus:

    ρ_B(y) = (Ω_B / (Ω_M - Ω_B)) × g(y) × ρ_DM(y) 
           + (Ω_B / (Ω_M - Ω_B)) × (c³ / 4π r_vir³) × (1/y²) × (dg/dy) × M_DM(<y)

This framework is agnostic to the specific dark matter profile and can encode
any simulation-based radial gas profile through the choice of g(y).

References:
    - Ayromlou et al. 2023 (closure radius concept)
    - Bryan & Norman 1998 (virial overdensity)
    - Navarro, Frenk & White 1997 (NFW profile)
"""

import numpy as np 
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from functools import partial

from astropy import units
from astropy import constants
from astropy.coordinates import SkyCoord, Angle

from frb.halos.models import rad3d2
from frb.defs import frb_cosmo as cosmo


class GeneralizedHaloProfile:
    """
    Generalized halo baryon profile using the g × M_DM framework.
    
    This class implements a baryon density profile that automatically 
    respects the cosmic baryon fraction constraint at large radii through 
    the modulation function g(y).
    
    Parameters
    ----------
    log_Mhalo : float
        Log10 of halo mass in solar masses. Default: 12.2
    c : float
        NFW concentration parameter. Default: 7.67
    z : float
        Redshift. Default: 0.
    cosmo : astropy.cosmology
        Cosmology object. Default: frb_cosmo
    k : float
        Scale parameter that sets where baryon fraction approaches cosmic 
        mean. The internal y_out is computed as c × k / arctanh(1 - tol).
        Larger k pushes the transition to larger radii. Default: 2.0
    g_form : str
        Functional form for g(y): 'tanh' or 'exp'. Default: 'tanh'
    zero_inner_ne : float
        If > 0, zero out electron density within this radius (in kpc).
    tol : float
        Tolerance for computing y_out from k. Default: 1e-6
    
    Attributes
    ----------
    rvir : Quantity
        Virial radius computed using Bryan & Norman 1998
    rho0 : Quantity
        NFW scale density
    cosmo_frac : float
        Ω_B / (Ω_M - Ω_B)
    y_out : float
        Computed scale parameter for g(y), derived from k and tol
    r_s : Quantity
        NFW scale radius (r_vir / c)
    Delta_vir : float
        Virial overdensity from Bryan & Norman 1998
    fb : float
        Standard cosmic baryon fraction Ω_B / Ω_M
    """
    
    def __init__(self,
                 log_Mhalo=12.2, 
                 c=7.67, 
                 z=0., 
                 cosmo=cosmo,
                 k=2.0,
                 g_form='tanh',
                 zero_inner_ne=0.,
                    tol=1e-6,
                 ):
        
        self.log_Mhalo = log_Mhalo
        self.M_halo = 10.**self.log_Mhalo * constants.M_sun.cgs
        self.c = c
        self.z = z
        self.cosmo = cosmo
        self.k = k
        self.y_out = self.c * self.k /jnp.arctanh(1-tol) 
        self.g_form = g_form
        self.zero_inner_ne = zero_inner_ne
        
        # Setup derived parameters
        self._setup_cosmology()
        self._setup_normalization()
    
    # =========================================================================
    # Setup Methods
    # =========================================================================
    
    def _setup_cosmology(self):
        """Setup cosmological parameters and virial quantities."""
        
        # Critical density at redshift z
        self.rhoc = self.cosmo.critical_density(self.z)
        self.H0 = self.cosmo.H0
        
        # Bryan & Norman 1998 virial overdensity for flat ΛCDM
        # Δ_vir = 18π² - 82q - 39q²
        # where q = Ω_Λ / [Ω_Λ + Ω_m(1+z)³]
        self.q = self.cosmo.Ode0 / (self.cosmo.Ode0 + self.cosmo.Om0 * (1 + self.z)**3)
        self.Delta_vir = 18 * np.pi**2 - 82 * self.q - 39 * self.q**2
        
        # Virial density and radius
        self.rhovir = self.Delta_vir * self.rhoc
        self.rvir = (((3 * self.M_halo) / (4 * np.pi * self.rhovir))**(1/3)).to('kpc')
        
        # Scale radius
        self.r_s = self.rvir / self.c
        
        # Cosmological baryon fraction
        # This is Ω_B / (Ω_M - Ω_B), the factor that relates M_B to M_DM
        self.cosmo_frac = self.cosmo.Ob0 / (self.cosmo.Om0 - self.cosmo.Ob0)
        
        # Standard baryon fraction for reference
        self.fb = self.cosmo.Ob0 / self.cosmo.Om0
        
        # Mean molecular weight (for n_H calculation)
        self.mu = 1.33
    
    def _setup_normalization(self):
        """Setup NFW normalization constant ρ_0."""
        
        # For NFW, M_DM(<r_vir) = M_halo (by definition, all DM is within r_vir)
        # M_DM(<y) = (4π r_s³ ρ_0) × f_DM(y)
        # At y = c: M_halo = (4π r_s³ ρ_0) × f_DM(c)
        
        f_dm_c = self.f_DM(self.c)
        four_pi_rs3 = 4 * np.pi * self.r_s**3
        
        self.rho0 = (self.M_halo / (four_pi_rs3 * f_dm_c)).to('g/cm^3')
        
        # Also store the constant for convenience
        self.const = four_pi_rs3
    
    # =========================================================================
    # g(y) Function and Derivatives
    # =========================================================================
    
    def g(self, y):
        """
        Modulation function g(y) that controls the baryon fraction profile.
        
        Constraints:
        - g(y) → 1 as y → ∞ (cosmic baryon fraction at large radii)
        - g(0) = 0 for both 'tanh' and 'exp' forms
        
        Parameters
        ----------
        y : float or ndarray
            Dimensionless radius y = c × r / r_vir
        
        Returns
        -------
        g_y : float or ndarray
            Value of g at y
        """
        y = jnp.atleast_1d(jnp.asarray(y))
        
        if self.g_form == 'tanh':
            # g(y) = tanh(y / y_out)
            # Approaches 1 as y → ∞, equals 0 at y = 0
            g_y = jnp.tanh(y / self.y_out)
            
        elif self.g_form == 'exp':
            raise ValueError("The 'exp' form is not implemented in this version.")
            
        else:
            raise ValueError(f"Unknown g_form: {self.g_form}")
        
        return g_y
    
    def _g_scalar(self, y):
        """Scalar version of g for JAX gradient."""
        if self.g_form == 'tanh':
            return jnp.tanh(y / self.y_out)
        elif self.g_form == 'exp':
            return 1 - jnp.exp(-y / self.y_out)
    
    def dg_dy(self, y):
        """
        Derivative of g(y) with respect to y.
        
        Parameters
        ----------
        y : float or ndarray
            Dimensionless radius
        
        Returns
        -------
        dg : float or ndarray
            dg/dy at y
        """
        y = jnp.atleast_1d(jnp.asarray(y))
        
        if self.g_form == 'tanh':
            # d/dy[tanh(y/y_out)] = sech²(y/y_out) / y_out
            dg = (1 / jnp.cosh(y / self.y_out))**2 / self.y_out
            
        elif self.g_form == 'exp':
            # d/dy[1 - exp(-y/y_out)] = exp(-y/y_out) / y_out
            dg = jnp.exp(-y / self.y_out) / self.y_out
            
        else:
            # Use JAX autodiff for custom functions
            grad_g = jax.vmap(jax.grad(self._g_scalar))
            dg = grad_g(y)
        
        return dg
    
    # =========================================================================
    # Dark Matter Profile (NFW)
    # =========================================================================
    
    def f_DM(self, y):
        """
        Dimensionless enclosed mass function for NFW profile.
        
        f_DM(y) = ln(1 + y) - y/(1 + y)
        
        Such that M_DM(<y) = (4π r_s³ ρ_0) × f_DM(y)
        
        Parameters
        ----------
        y : float or ndarray
            Dimensionless radius y = c × r / r_vir
        
        Returns
        -------
        f_y : float or ndarray
        """
        y = jnp.atleast_1d(jnp.asarray(y))
        f_y = jnp.log(1 + y) - y / (1 + y)
        return f_y
    
    def rho_DM(self, y):
        """
        NFW dark matter density profile.
        
        ρ_DM(y) = ρ_0 / [y × (1 + y)²]
        
        Parameters
        ----------
        y : float or ndarray
            Dimensionless radius
        
        Returns
        -------
        rho : Quantity
            Dark matter density in g/cm³
        """
        y = jnp.atleast_1d(jnp.asarray(y))
        # Avoid division by zero at y = 0
        y_safe = jnp.where(y > 0, y, 1e-10)
        rho = self.rho0.value / (y_safe * (1 + y_safe)**2)

        return rho * self.rho0.unit
    

    def rho_dm(self, xyz):
        """
        Dark matter density at position xyz.
        
        Parameters
        ----------
        xyz : ndarray
            Position in kpc, shape (3,) or (3, npoints)
        
        Returns
        -------
        rho : Quantity
            Dark matter density in g/cm³
        """
        radius = np.sqrt(rad3d2(xyz))  # kpc
        y = self.c * (radius / self.rvir.to('kpc').value)
        return self.rho_DM(y)
    
    def M_DM(self, y):
        """
        Enclosed dark matter mass within dimensionless radius y.
        
        M_DM(<y) = (4π r_s³ ρ_0) × f_DM(y)
        
        Parameters
        ----------
        y : float or ndarray
            Dimensionless radius
        
        Returns
        -------
        M : Quantity
            Enclosed dark matter mass in solar masses
        """
        return (self.const * self.rho0 * self.f_DM(y)).to('Msun')
    
    # =========================================================================
    # Baryon Profile (g × M_DM Framework)
    # =========================================================================
    
    def M_B(self, y):
        """
        Enclosed baryon mass within dimensionless radius y.
        
        M_B(<y) = g(y) × (Ω_B / (Ω_M - Ω_B)) × M_DM(<y)
        
        Parameters
        ----------
        y : float or ndarray
            Dimensionless radius
        
        Returns
        -------
        M : Quantity
            Enclosed baryonic mass in solar masses
        """
        return self.cosmo_frac * self.g(y) * self.M_DM(y)
    
    def rho_B(self, y):
        """
        Baryon density profile derived from the g × M_DM framework.
        
        ρ_B(y) = (Ω_B / (Ω_M - Ω_B)) × [g(y) × ρ_DM(y) 
                 + (c³ / 4π r_vir³) × (1/y²) × (dg/dy) × M_DM(<y)]
        
        The first term is a scaled dark matter profile modulated by g(y).
        The second term is a "correction" ensuring proper enclosed mass behavior.
        
        Parameters
        ----------
        y : float or ndarray
            Dimensionless radius
        
        Returns
        -------
        rho : Quantity
            Baryon density in g/cm³
        """
        y = jnp.atleast_1d(jnp.asarray(y))
        
        # Avoid division by zero
        y_safe = jnp.where(y > 0, y, 1e-10)
        
        # Term 1: g(y) × ρ_DM(y)
        term1 = self.g(y) * self.rho_DM(y_safe).value
        
        # Term 2: (1/y²) × (dg/dy) × M_DM(<y) × (c³ / 4π r_vir³)
        # Note: (c³ / 4π r_vir³) = 1 / (4π r_s³) = 1 / const
        dg = self.dg_dy(y)
        M_dm = self.M_DM(y_safe).to('g').value
        term2 = (1 / y_safe**2) * dg * M_dm / self.const.to('cm^3').value
        
        # Combine with cosmological factor
        rho = self.cosmo_frac * (term1 + term2)
        
        return rho * units.g / units.cm**3
    
    def f_B_enclosed(self, y):
        """
        Enclosed baryon fraction f_B(<y) = M_B(<y) / [M_DM(<y) + M_B(<y)]
        
        This should approach Ω_B / Ω_M as y → ∞.
        
        Parameters
        ----------
        y : float or ndarray
            Dimensionless radius
        
        Returns
        -------
        f_B : float or ndarray
            Enclosed baryon fraction
        """
        M_b = self.M_B(y).value
        M_dm = self.M_DM(y).value
        return M_b / (M_dm + M_b)
    
    # =========================================================================
    # Interface Methods (for compatibility with existing code)
    # =========================================================================
    
    def rho_b(self, xyz):
        """
        Baryon density at position xyz (for compatibility with existing code).
        
        Parameters
        ----------
        xyz : ndarray
            Position in kpc, shape (3,) or (3, npoints)
        
        Returns
        -------
        rho : Quantity
            Baryon density in g/cm³
        """
        radius = np.sqrt(rad3d2(xyz))  # kpc
        y = self.c * (radius / self.rvir.to('kpc').value)
        return self.rho_B(y)
    
    def nH(self, xyz):
        """
        Hydrogen number density at position xyz.
        
        Parameters
        ----------
        xyz : ndarray
            Position in kpc, shape (3,) or (3, npoints)
        
        Returns
        -------
        nH : float or ndarray
            Hydrogen number density in cm⁻³
        """
        m_p = constants.m_p.cgs  # g (keep units)
        rho = self.rho_b(xyz)
        nH = (rho / self.mu / m_p).to('cm^-3').value
        return nH
    
    def ne(self, xyz):
        """
        Electron number density at position xyz.
        
        Assumes fully ionized gas with 25% Helium by mass.
        n_e = 1.1667 × n_H
        
        Parameters
        ----------
        xyz : ndarray
            Position in kpc, shape (3,) or (3, npoints)
        
        Returns
        -------
        ne : float or ndarray
            Electron density in cm⁻³
        
        Notes
        -----
        If zero_inner_ne > 0, electron density is set to zero within
        that radius to model feedback-evacuated cores.
        """
        ne = self.nH(xyz) * 1.1667
        
        # Zero out inner region if specified
        if self.zero_inner_ne > 0.:
            rad = np.sum(xyz**2, axis=0)
            inner = rad < self.zero_inner_ne**2
            if np.any(inner):
                if len(xyz.shape) == 1:
                    ne = 0.
                else:
                    ne[inner] = 0.
        
        return ne
    
    def Ne_Rperp(self, Rperp, step_size=0.1*units.kpc, rmax=1., add_units=True):
        """
        Calculate DM (electron column density) at impact parameter Rperp.
        
        Integrates electron density along a sightline at perpendicular 
        distance Rperp from the halo center.
        
        Parameters
        ----------
        Rperp : Quantity
            Impact parameter in kpc
        step_size : Quantity
            Integration step size. Default: 0.1 kpc
        rmax : float
            Maximum radius in units of r_vir. Default: 1.0
        add_units : bool
            Whether to return with units. Default: True
        
        Returns
        -------
        Ne : float or Quantity
            Dispersion measure contribution in pc/cm³
        """
        dz = step_size.to('kpc').value
        
        # Check if beyond halo
        if Rperp > rmax * self.rvir:
            return 0. * units.pc / units.cm**3 if add_units else 0.
        
        # Generate sightline through halo
        zmax = np.sqrt((rmax * self.rvir)**2 - Rperp**2).to('kpc')
        zval = np.arange(-zmax.value, zmax.value + dz, dz)
        
        # Set xyz coordinates
        xyz = np.zeros((3, zval.size))
        xyz[0, :] = Rperp.to('kpc').value
        xyz[2, :] = zval
        
        # Integrate n_e along sightline
        ne = self.ne(xyz)  # cm⁻³
        Ne = np.sum(ne) * dz * 1000  # pc cm⁻³
        
        if add_units:
            return Ne * units.pc / units.cm**3
        return Ne
    
    def cal_Mb_at_r(self, r):
        """
        Calculate enclosed baryonic mass at radius r.
        
        Parameters
        ----------
        r : Quantity
            Radius in kpc
        
        Returns
        -------
        M_b : Quantity
            Enclosed baryonic mass in solar masses
        """
        y = self.c * r / self.rvir
        return self.M_B(y.value)
    
    def cal_Mdm_at_r(self, r):
        """
        Calculate enclosed dark matter mass at radius r.
        
        Parameters
        ----------
        r : Quantity
            Radius in kpc
        
        Returns
        -------
        M_dm : Quantity
            Enclosed dark matter mass in solar masses
        """
        y = self.c * r / self.rvir
        return self.M_DM(y.value)
    
    # =========================================================================
    # Utility Methods
    # =========================================================================
    
    def __repr__(self):
        return (f"<{self.__class__.__name__}: "
                f"log_Mhalo={self.log_Mhalo:.2f}, "
                f"c={self.c:.2f}, "
                f"k={self.k:.2f}, "
                f"y_out={self.y_out:.2f}, "
                f"g_form='{self.g_form}', "
                f"r_vir={self.rvir:.1f}>")


# =============================================================================
# Convenience subclasses with specific g(y) forms
# =============================================================================

class TanhHaloProfile(GeneralizedHaloProfile):
    """
    Halo profile with g(y) = tanh(y / y_out).
    
    This form ensures:
    - g(0) = 0 (baryon-to-DM mass ratio starts at zero)
    - g(y) → 1 as y → ∞ (cosmic baryon fraction at large r)
    - g(y) < 1 for all finite y (enclosed baryon fraction always below cosmic mean)
    
    The tanh form provides a smooth, monotonic transition from baryon-depleted
    inner regions to cosmic baryon fraction at large radii.
    """
    
    def __init__(self, **kwargs):
        kwargs['g_form'] = 'tanh'
        super().__init__(**kwargs)


class ExpHaloProfile(GeneralizedHaloProfile):
    """
    Halo profile with g(y) = 1 - exp(-y / y_out).
    
    Similar properties to tanh:
    - g(0) = 0 (baryon-to-DM mass ratio starts at zero)
    - g(y) → 1 as y → ∞ (cosmic baryon fraction at large r)
    
    The exponential form approaches unity more slowly than tanh,
    which may better match certain simulation results.
    """
    
    def __init__(self, **kwargs):
        kwargs['g_form'] = 'exp'
        super().__init__(**kwargs)


# =============================================================================
# M31 
# =============================================================================

class M31_GeneralizedHaloProfile(GeneralizedHaloProfile):
    """
    M31-specific generalized halo profile.
    
    This subclass sets default parameters appropriate for M31
    and provides convenience methods for M31 analyses.
    """
    
    def __init__(self,
                 log_Mhalo=12.2,
                 c=7.67,
                 z=0.,
                 cosmo=cosmo,
                 k=2.0,
                 g_form='tanh',
                 zero_inner_ne=0.,
                 tol=1e-6,
                 ):
        super().__init__(log_Mhalo=log_Mhalo,
                         c=c,
                         z=z,
                         cosmo=cosmo,
                         k=k,
                         g_form=g_form,
                         zero_inner_ne=zero_inner_ne,
                         tol=tol)
        
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


    def DM_from_impact_param_b(self, bimpact, **kwargs):
        """
        Calculate DM through M31's halo from the Sun
        given an impact parameter

        Args:
            bimpact: Quantity
                Ratio of the impact parameter to r200
            **kwargs:
               Passed to Ne_Rperp

        Returns:
            DM: Quantity
              Dispersion measure through M31's halo
        """
        a=1
        c=0
        x0, y0 = self.distance.to('kpc').value, 0. # kpc
        # Calculate r200_rad
        r200_rad = (self.r200 / self.distance.to('kpc'))*units.rad

        # Create an Angle object for sep
        sep = Angle(bimpact * r200_rad, unit='radian')

        # More geometry
        atan = np.arctan(sep.radian)
        b = -1 * a / atan

        # Restrct to within 90deg (everything beyond is 0 anyhow)
        if sep > 90.*units.deg:
            return 0 * units.pc / units.cm**3
        # Rperp
        Rperp = np.abs(a * x0 + b * y0 + c) / np.sqrt(a**2 + b**2)  # kpc

        DM = self.Ne_Rperp(Rperp*units.kpc, **kwargs).to('pc/cm**3')
        return DM