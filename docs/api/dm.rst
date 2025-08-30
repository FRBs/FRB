frb.dm - Dispersion Measure Module
===================================

The ``frb.dm`` module provides comprehensive tools for calculating and modeling dispersion measure contributions from various sources along the line of sight to Fast Radio Bursts.

.. currentmodule:: frb.dm

Module Overview
---------------

The dispersion measure (DM) is a key observable for FRBs, representing the integrated electron column density along the line of sight:

.. math::

   \text{DM} = \int_0^{d} n_e(l) \, dl

This module decomposes the total DM into contributions from:

- **Milky Way**: Galactic disk and halo electrons
- **Intergalactic Medium (IGM)**: Diffuse cosmic electrons  
- **Host Galaxy**: Electrons in the FRB host galaxy
- **Intervening Systems**: Galaxy halos and other structures

Submodules
----------

.. toctree::
   :maxdepth: 1

   dm_igm
   dm_cosmic
   dm_host
   dm_mcmc
   dm_prob_dmz

frb.dm.igm - Intergalactic Medium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: frb.dm.igm
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
^^^^^^^^^^^^^

.. autofunction:: frb.dm.igm.DM_cosmic

.. autofunction:: frb.dm.igm.ne_cosmic

.. autofunction:: frb.dm.igm.z_from_DM

.. autofunction:: frb.dm.igm.DM_halos

frb.dm.cosmic - Cosmic DM Modeling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: frb.dm.cosmic
   :members:
   :undoc-members:
   :show-inheritance:

frb.dm.host - Host Galaxy Contributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: frb.dm.host
   :members:
   :undoc-members:
   :show-inheritance:

frb.dm.mcmc - MCMC Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: frb.dm.mcmc
   :members:
   :undoc-members:
   :show-inheritance:

frb.dm.prob_dmz - DM-Redshift Probabilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: frb.dm.prob_dmz
   :members:
   :undoc-members:
   :show-inheritance:

Detailed Function Documentation
-------------------------------

IGM Functions
~~~~~~~~~~~~~

.. py:function:: ne_cosmic(z, cosmo=defs.frb_cosmo, mu=4./3)

   Calculate the average cosmic electron number density as a function of redshift.
   
   This function computes the physical number density of electrons in the universe,
   accounting for the cosmic baryon density and ionization fraction.

   :param float z: Redshift
   :param astropy.cosmology.Cosmology cosmo: Cosmology in which the calculations are performed
   :param float mu: Reduced mass (default: 4/3 for fully ionized H and He)
   :returns: Average physical number density of electrons in the universe
   :rtype: astropy.units.Quantity
   
   **Example:**
   
   .. code-block:: python
   
      from frb.dm.igm import ne_cosmic
      
      # Calculate electron density at z=1
      z = 1.0
      n_e = ne_cosmic(z)
      print(f"Electron density at z={z}: {n_e}")

.. py:function:: DM_cosmic(z, cosmo=defs.frb_cosmo, cumul=False)

   Calculate the cosmic dispersion measure contribution.
   
   Integrates the cosmic electron density from z=0 to the specified redshift,
   including corrections for cosmological expansion.

   :param float z: Redshift of the FRB
   :param astropy.cosmology.Cosmology cosmo: Cosmology for calculations
   :param bool cumul: Return cumulative DM evolution with redshift
   :returns: Cosmic dispersion measure or array if cumul=True
   :rtype: astropy.units.Quantity
   
   **Example:**
   
   .. code-block:: python
   
      from frb.dm.igm import DM_cosmic
      
      # Calculate cosmic DM to z=0.5
      z_frb = 0.5
      dm_cosmic = DM_cosmic(z_frb)
      print(f"Cosmic DM to z={z_frb}: {dm_cosmic}")

.. py:function:: z_from_DM(DM, cosmo=defs.frb_cosmo, coord=None, corr_nuisance=True)

   Estimate redshift from observed dispersion measure.
   
   Uses the cosmic DM-redshift relation to estimate the redshift of an FRB
   given its observed dispersion measure, after correcting for Galactic
   and other contributions.

   :param astropy.units.Quantity DM: Observed dispersion measure
   :param astropy.cosmology.Cosmology cosmo: Cosmology
   :param astropy.coordinates.SkyCoord coord: Sky coordinates (optional)
   :param bool corr_nuisance: Apply corrections for nuisance parameters
   :returns: Estimated redshift
   :rtype: float
   
   **Example:**
   
   .. code-block:: python
   
      from frb.dm.igm import z_from_DM
      from astropy import units as u
      
      # Estimate redshift from DM
      DM_obs = 500 * u.pc / u.cm**3
      z_est = z_from_DM(DM_obs)
      print(f"Estimated redshift: {z_est:.2f}")

.. py:function:: DM_halos(z, cosmo, f_hot=0.75, rmax=2., logMmin=10.3, logMmax=16., neval=500, cumul=False)

   Calculate dispersion measure contribution from galaxy halos.
   
   Models the DM contribution from hot gas in galaxy halos along the
   line of sight, using a halo mass function approach.

   :param float z: Redshift of the FRB  
   :param astropy.cosmology.Cosmology cosmo: Cosmology for calculations
   :param float f_hot: Fraction of halo baryons in diffuse phase
   :param float rmax: Size of halo in units of r200
   :param float logMmin: Log of minimum halo mass (cannot be much below 10.3)
   :param float logMmax: Log of maximum halo mass (default: 10^16 Msun)
   :param int neval: Number of redshift evaluation points
   :param bool cumul: Return cumulative evaluation
   :returns: DM contribution from halos
   :rtype: astropy.units.Quantity

Utility Functions
~~~~~~~~~~~~~~~~~

.. py:function:: f_diffuse(z, cosmo=defs.frb_cosmo, return_rho=False)

   Calculate the fraction of baryons in the diffuse IGM.

   :param float z: Redshift
   :param astropy.cosmology.Cosmology cosmo: Cosmology  
   :param bool return_rho: Also return the diffuse gas density
   :returns: Diffuse fraction (and density if requested)
   :rtype: float or tuple

.. py:function:: sigma_DM_cosmic(DM_cosmic, rel_err_Mstar=0.1)

   Calculate uncertainty in cosmic DM due to stellar mass uncertainties.

   :param astropy.units.Quantity DM_cosmic: Cosmic dispersion measure
   :param float rel_err_Mstar: Relative error in stellar mass
   :returns: DM uncertainty
   :rtype: astropy.units.Quantity

Physical Models
---------------

Ionization Models
~~~~~~~~~~~~~~~~~

The module includes models for cosmic reionization history:

**Hydrogen Reionization:**
- Reionization redshift: z_HI ~ 6-8
- Affects IGM electron fraction at high redshift

**Helium Reionization:**  
- HeII reionization: z_HeII ~ 3-4
- Increases electron density at intermediate redshift

**Implementation:**

.. code-block:: python

   from frb.dm.igm import f_ionization_He
   
   # Calculate helium ionization fraction
   z = 3.0
   f_He = f_ionization_He(z)
   print(f"Helium ionization fraction at z={z}: {f_He}")

Halo Models
~~~~~~~~~~~

Galaxy halo contributions use the halo mass function and hot gas profiles:

**Key Parameters:**
- Halo mass range: 10^10.3 to 10^16 M_sun
- Hot gas fraction: ~75% of cosmic baryon fraction
- Radial extent: typically 2 × r200

**Scaling Relations:**
- Gas density profiles (e.g., beta models)
- Mass-temperature relations
- Feedback effects on gas distribution

Usage Examples
--------------

Basic DM Analysis
~~~~~~~~~~~~~~~~~

.. code-block:: python

   from frb.dm import igm
   from astropy import units as u
   
   # Calculate total cosmic DM budget
   z_max = 2.0
   dm_cosmic_total = igm.DM_cosmic(z_max)
   
   # Calculate contributions at different redshifts
   redshifts = [0.1, 0.5, 1.0, 1.5, 2.0]
   for z in redshifts:
       dm_z = igm.DM_cosmic(z)
       print(f"DM(z={z:.1f}) = {dm_z:.0f}")

Redshift Estimation
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from frb.dm.igm import z_from_DM
   from astropy import units as u
   from astropy.coordinates import SkyCoord
   
   # Observed FRB parameters
   DM_obs = 750 * u.pc / u.cm**3
   coord = SkyCoord(ra=123.45*u.deg, dec=12.34*u.deg)
   
   # Estimate redshift
   z_est = z_from_DM(DM_obs, coord=coord)
   print(f"Estimated redshift: {z_est:.2f}")
   
   # Calculate expected cosmic DM at this redshift
   dm_cosmic_expected = igm.DM_cosmic(z_est)
   print(f"Expected cosmic DM: {dm_cosmic_expected:.0f}")
   
   # Calculate excess DM (host + local contributions)
   dm_excess = DM_obs - dm_cosmic_expected
   print(f"Excess DM: {dm_excess:.0f}")

Halo Contribution Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from frb.dm.igm import DM_halos
   from astropy.cosmology import Planck18
   
   # Calculate halo DM contribution
   z_frb = 1.0
   dm_halos = DM_halos(z_frb, Planck18)
   
   print(f"Halo DM contribution to z={z_frb}: {dm_halos:.1f}")
   
   # Compare different halo parameters
   f_hot_values = [0.5, 0.75, 1.0]
   for f_hot in f_hot_values:
       dm_h = DM_halos(z_frb, Planck18, f_hot=f_hot)
       print(f"Halo DM (f_hot={f_hot}): {dm_h:.1f}")

Advanced Usage
--------------

Custom Cosmologies
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from astropy.cosmology import FlatLambdaCDM
   from frb.dm.igm import DM_cosmic
   
   # Define custom cosmology
   custom_cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
   
   # Compare DM calculations
   z = 1.0
   dm_planck = DM_cosmic(z)  # Default Planck18
   dm_custom = DM_cosmic(z, cosmo=custom_cosmo)
   
   print(f"DM (Planck18): {dm_planck:.0f}")
   print(f"DM (Custom):   {dm_custom:.0f}")

Monte Carlo Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   from frb.dm.igm import z_from_DM
   from astropy import units as u
   
   # Simulate DM uncertainties
   DM_mean = 500 * u.pc / u.cm**3
   DM_err = 50 * u.pc / u.cm**3
   n_trials = 1000
   
   # Generate DM samples
   DM_samples = np.random.normal(
       DM_mean.value, DM_err.value, n_trials
   ) * DM_mean.unit
   
   # Calculate redshift distribution
   z_samples = [z_from_DM(dm) for dm in DM_samples]
   
   z_mean = np.mean(z_samples)
   z_std = np.std(z_samples)
   print(f"Redshift: {z_mean:.2f} ± {z_std:.2f}")

Error Handling
--------------

The DM module includes robust error handling:

.. code-block:: python

   from frb.dm.igm import z_from_DM
   from astropy import units as u
   
   try:
       # This might fail for very high DM values
       z = z_from_DM(10000 * u.pc / u.cm**3)
   except ValueError as e:
       print(f"DM too high: {e}")
   
   try:
       # This might fail for negative DM
       z = z_from_DM(-100 * u.pc / u.cm**3) 
   except ValueError as e:
       print(f"Invalid DM: {e}")

Performance Notes
-----------------

**Computational Efficiency:**
- DM calculations use vectorized numpy operations
- Cosmological integrals are cached for repeated calls
- Halo mass function integration is optimized for speed

**Memory Usage:**
- Large redshift arrays may require significant memory
- Use ``cumul=False`` for single-point calculations
- Consider chunking for very large datasets

**Accuracy:**
- Integration tolerances balance speed and precision  
- Default parameters provide ~1% accuracy for most applications
- High-precision calculations may require custom tolerances

See Also
--------

- :doc:`../tutorials` - Detailed usage examples
- :doc:`turb_scattering` - Scattering analysis
- :doc:`frbcat` - Catalogue operations
- :doc:`../examples` - Real-world applications

.. note::
   The DM module assumes a standard cosmological model by default.
   For non-standard cosmologies, explicitly pass the cosmology parameter
   to relevant functions.