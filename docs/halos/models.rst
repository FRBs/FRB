Halo Modeling
=============

.. automodule:: frb.halos.models

Module for DM halo calculations and galaxy halo models.

Constants
---------

.. autodata:: m_p

   Proton mass in CGS units for optimized calculations.

Functions
---------

init_hmf
~~~~~~~~

.. autofunction:: init_hmf

   Initialize the Aemulus Halo Mass Function.

   .. warning::
      Uses the original version which codes Tinker+2008. May be refactored to use the more accurate new version.

   **Returns:**
   
   * ``hmf_emulator`` - Initialized HMF emulator object

frac_in_halos
~~~~~~~~~~~~~

.. autofunction:: frac_in_halos

   Calculate the fraction of matter in collapsed halos over a mass range and at a given redshift.

   .. note::
      The fraction of DM associated with these halos will be scaled down by an additional factor of f_diffuse.

   **Parameters:**
   
   * ``zvals`` : ndarray - Redshift values
   * ``Mlow`` : float - Minimum halo mass in h^-1 units
   * ``Mhigh`` : float - Maximum halo mass in h^-1 units  
   * ``rmax`` : float, optional - Extent of halo in units of rvir (default: 1.0)

   **Returns:**
   
   * ``ratios`` : ndarray - rho_halo / rho_m

   .. warning::
      This calculation assumes a single concentration for all halos.

stellarmass_from_halomass
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: stellarmass_from_halomass

   Stellar mass from Halo Mass using Moster+2013 relation.

   **Parameters:**
   
   * ``log_Mhalo`` : float - log_10 halo mass in solar mass units
   * ``z`` : float, optional - Halo redshift (default: 0)
   * ``params`` : list, optional - Custom model parameters

   **Returns:**
   
   * ``log_mstar`` : float - log_10 galaxy stellar mass in solar mass units

   **Reference:** https://doi.org/10.1093/mnras/sts261

halomass_from_stellarmass
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: halomass_from_stellarmass

   Halo mass from Stellar mass (Moster+2013). Numerically inverts stellarmass_from_halomass.

   **Parameters:**
   
   * ``log_mstar`` : float or ndarray - log_10 stellar mass in solar mass units
   * ``z`` : float, optional - Galaxy redshift (default: 0)
   * ``randomize`` : bool, optional - Add scatter to the relation

   **Returns:**
   
   * ``log_Mhalo`` : float - log_10 halo mass in solar mass units

stellarmass_from_halomass_kravtsov
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: stellarmass_from_halomass_kravtsov

   Stellar mass from Halo Mass using Kravtsov+2018 relation.

   .. caution::
      This relation is valid for low z (z~0). Higher z values may require a scaled relation.

   **Parameters:**
   
   * ``log_mhalo`` : float - log_10 halo mass

   **Returns:**
   
   * ``log_mstar`` : float - log_10 galaxy stellar mass

   **Reference:** https://ui.adsabs.harvard.edu/abs/2018AstL...44....8K/abstract

halomass_from_stellarmass_kravtsov
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: halomass_from_stellarmass_kravtsov

   Inverts stellarmass_from_halomass_kravtsov function.

   **Parameters:**
   
   * ``log_mstar`` : float or ndarray - log_10 stellar mass

   **Returns:**
   
   * ``log_mhalo`` : float - log_10 halo mass

rad3d2
~~~~~~

.. autofunction:: rad3d2

   Calculate radius to x,y,z coordinates. Assumes origin is (0,0,0).

   **Parameters:**
   
   * ``xyz`` : tuple or ndarray - 3D coordinates

   **Returns:**
   
   * ``rad3d`` : float or ndarray - 3D radius squared

Classes
-------

ModifiedNFW
~~~~~~~~~~~

.. autoclass:: ModifiedNFW
   :members:
   :undoc-members:
   :show-inheritance:

   Generate a modified NFW model for hot, virialized gas (e.g. Mathews & Prochaska 2017).

   **Parameters:**

   * ``log_Mhalo`` : float, optional - log10 of halo mass in solar masses (default: 12.2)
   * ``c`` : float, optional - Concentration of the halo (default: 7.67)
   * ``f_hot`` : float, optional - Fraction of baryons in hot phase (default: 0.75)
   * ``alpha`` : float, optional - Parameter to modify NFW profile power-law (default: 0)
   * ``y0`` : float, optional - Parameter to modify NFW profile position (default: 1)
   * ``z`` : float, optional - Redshift of the halo (default: 0)
   * ``cosmo`` : astropy cosmology, optional - Cosmology of the universe

   **Key Attributes:**

   * ``H0`` : Quantity - Hubble constant
   * ``fb`` : float - Cosmic fraction of baryons (default: 0.16)
   * ``r200`` : Quantity - Virial radius
   * ``rho0`` : Quantity - Density normalization
   * ``M_b`` : Quantity - Mass in baryons

   **Methods:**

   .. automethod:: setup_param

      Setup key parameters of the model.

   .. automethod:: Ne_Rperp

      Calculate column density along a perpendicular path through the halo.

ICM
~~~

.. autoclass:: ICM
   :members:
   :undoc-members:
   :show-inheritance:

   Intracluster Medium model, child of ModifiedNFW.

   Implements the Miller & Bregman 2015 ICM model for galaxy clusters.

   **Methods:**

   .. automethod:: nH

      Calculate the number density of Hydrogen.

      **Parameters:**
      
      * ``xyz`` : ndarray - Coordinates in kpc

      **Returns:**
      
      * ``ndarray`` - Number density in units of 1/cm³

MilkyWay
~~~~~~~~

.. autoclass:: MilkyWay
   :members:
   :undoc-members:
   :show-inheritance:

   Fiducial model for the Galaxy. Halo mass follows latest constraints.

   Density profile is similar to Maller & Bullock 2004.

M31
~~~

.. autoclass:: M31
   :members:
   :undoc-members:
   :show-inheritance:

   Preferred model for M31. Mass from van der Marel 2012.

   **Attributes:**

   * ``distance`` : Quantity - Distance from Sun (752 kpc)
   * ``coord`` : SkyCoord - Coordinates of M31

   **Methods:**

   .. automethod:: DM_from_Galactic

      Calculate DM through M31's halo from the Sun given a direction.

      **Parameters:**
      
      * ``scoord`` : SkyCoord - Coordinates of the sightline

      **Returns:**
      
      * ``DM`` : Quantity - Dispersion measure through M31's halo

Examples
--------

Basic halo model usage::

    from frb.halos.models import ModifiedNFW, halomass_from_stellarmass
    from astropy import units as u
    
    # Create a halo model
    halo = ModifiedNFW(log_Mhalo=12.5, z=0.3, f_hot=0.6)
    
    # Calculate DM at 100 kpc offset
    offset = 100 * u.kpc
    dm_contribution = halo.Ne_Rperp(offset)
    print(f"DM contribution: {dm_contribution}")

Stellar-halo mass relations::

    from frb.halos.models import stellarmass_from_halomass, halomass_from_stellarmass
    
    # Convert halo mass to stellar mass
    log_mhalo = 12.0  # log solar masses
    log_mstar = stellarmass_from_halomass(log_mhalo, z=0.5)
    
    # Invert the relation
    log_mhalo_recovered = halomass_from_stellarmass(log_mstar, z=0.5)
    
    print(f"Halo mass: {log_mhalo}")
    print(f"Stellar mass: {log_mstar:.2f}")
    print(f"Recovered halo mass: {log_mhalo_recovered:.2f}")

Galaxy-specific models::

    from frb.halos.models import MilkyWay, M31
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    # Milky Way model
    mw = MilkyWay(log_Mhalo=12.1, f_hot=0.7)
    
    # M31 model with DM calculation
    m31 = M31(log_Mhalo=12.2)
    target_coord = SkyCoord('01h33m51s', '+30d39m37s')
    dm_m31 = m31.DM_from_Galactic(target_coord)
    print(f"DM through M31: {dm_m31}")

Halo mass function calculations::

    from frb.halos.models import frac_in_halos
    import numpy as np
    
    # Calculate matter fraction in halos
    redshifts = np.array([0.1, 0.5, 1.0])
    mass_min = 1e11  # Solar masses
    mass_max = 1e15
    
    fractions = frac_in_halos(redshifts, mass_min, mass_max)
    
    for z, frac in zip(redshifts, fractions):
        print(f"z={z:.1f}: {frac:.3f} of matter in halos")

Advanced usage with scatter::

    # Include scatter in stellar-halo mass relation
    log_mstar = 10.5
    z = 0.3
    
    # Multiple realizations with scatter
    log_mhalos = []
    for i in range(100):
        log_mhalo = halomass_from_stellarmass(log_mstar, z=z, randomize=True)
        log_mhalos.append(log_mhalo)
    
    mean_mhalo = np.mean(log_mhalos)
    std_mhalo = np.std(log_mhalos)
    print(f"Mean log(Mhalo): {mean_mhalo:.2f} ± {std_mhalo:.2f}")