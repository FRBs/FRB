""" Module for the modified-modified NFW (mmNFW) baryon density profile.

This implements the generalized halo baryon profile framework described in
the derivation where the baryon mass is defined as:
    M_B(<y) = g(y) * (Omega_B / (Omega_M - Omega_B)) * M_DM(<y)

With g(y) = tanh(y/y_out) to ensure the enclosed mass ratio asymptotes
to the cosmic mean at large radii while remaining below it at all radii.

References:
    - Matthews & Prochaska 2017
    - Prochaska & Zheng 2019
    - Ayromlou+2023 (closure radius concept)
"""
import numpy as np

from astropy import units
from astropy import constants

from frb.halos.models import ModifiedNFW
from frb.defs import frb_cosmo as cosmo


class mmNFW(ModifiedNFW):
    """
    Modified-modified NFW profile implementing the g*M_DM framework.

    This profile defines baryon density such that the enclosed baryon to
    total matter ratio asymptotically approaches the cosmic mean (Omega_B/Omega_M)
    at large radii, while remaining below it at all finite radii.

    The key innovation is introducing g(y) = tanh(y/y_out) which:
    - Ensures M_B/M_total < Omega_B/Omega_M for all y (no spurious crossings)
    - Asymptotes to the cosmic mean as y -> infinity
    - Introduces an outer scale radius y_out beyond which the ratio approaches unity

    Parameters
    ----------
    y_out : float, optional
        Outer scale parameter in units of y = c*r/r_vir. Controls the radius
        beyond which the enclosed mass ratio asymptotes to the cosmic mean.
        Default is 2.0.
    log_Mhalo : float, optional
        log10 of the halo mass in solar masses. Default is 12.2.
    c : float, optional
        NFW concentration parameter. Default is 7.67.
    z : float, optional
        Redshift of the halo. Default is 0.
    cosmo : astropy.cosmology, optional
        Cosmology to use. Default is frb_cosmo.

    Attributes
    ----------
    y_out : float
        Outer scale parameter
    rho0 : Quantity
        Dark matter central density normalization
    M_DM0 : Quantity
        Dark matter mass normalization (4*pi*rho0*(r_vir/c)^3)
    """

    def __init__(self, y_out=2.0, log_Mhalo=12.2, c=7.67, z=0.,
                 cosmo=cosmo, **kwargs):

        self.y_out = y_out

        # Initialize parent class (sets up r200, rho0, etc.)
        # We don't use f_hot here as the baryon mass is set by the g*M_DM framework
        super().__init__(log_Mhalo=log_Mhalo, c=c, z=z,
                        cosmo=cosmo, f_hot=1.0, **kwargs)

        # Compute M_DM0 = 4*pi*rho0*(r_vir/c)^3
        self.M_DM0 = 4 * np.pi * self.rho0 * (self.r200 / self.c)**3

    def rho_b(self, xyz):
        """
        Baryon density profile from the g*M_DM framework.

        Implements the equation:
        rho_B(y) = (Omega_B/(Omega_M - Omega_B)) * [
            (c^3 * M_DM0 / (4*pi*r_vir^3)) * (sech^2(y/y_out) / (y^2 * y_out))
                * (ln(1+y) - y/(1+y))
            + tanh(y/y_out) * rho0 / (y*(1+y)^2)
        ]

        Parameters
        ----------
        xyz : ndarray
            Position array (3, npoints) in kpc, or (3,) for single point

        Returns
        -------
        rho : Quantity
            Baryon mass density in g/cm^3
        """
        from frb.halos.models import rad3d2

        radius = np.sqrt(rad3d2(xyz))
        y = self.c * (radius / self.r200.to('kpc').value)

        # Avoid division by zero at y=0
        y = np.maximum(y, 1e-10)

        # Cosmological ratio
        Omega_B = self.cosmo.Ob0
        Omega_M = self.cosmo.Om0
        cosmo_ratio = Omega_B / (Omega_M - Omega_B)

        # Prefactor for first term: c^3 * M_DM0 / (4*pi*r_vir^3)
        prefactor1 = (self.c**3 * self.M_DM0 /
                     (4 * np.pi * self.r200**3)).to('g/cm**3')

        # sech^2(y/y_out) = 1/cosh^2(y/y_out)
        sech2 = 1.0 / np.cosh(y / self.y_out)**2

        # NFW enclosed mass factor: ln(1+y) - y/(1+y)
        fy_dm = np.log(1 + y) - y / (1 + y)

        # First term: (sech^2(y/y_out) / (y^2 * y_out)) * (ln(1+y) - y/(1+y))
        term1 = prefactor1 * sech2 / (y**2 * self.y_out) * fy_dm

        # Second term: tanh(y/y_out) * rho0 / (y*(1+y)^2)
        tanh_term = np.tanh(y / self.y_out)
        term2 = tanh_term * self.rho0 / (y * (1 + y)**2)

        # Full baryon density
        rho = cosmo_ratio * (term1 + term2)

        return rho.to('g/cm**3')
