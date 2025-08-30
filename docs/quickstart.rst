Quick Start Guide
=================

This guide will get you up and running with the FRB package in just a few minutes.

Basic Usage
-----------

Loading FRBs
~~~~~~~~~~~~

The most common starting point is loading a known FRB by name:

.. code-block:: python

   import frb as ffrb
   
   # Load a specific FRB by name
   frb121102 = ffrb.FRB.by_name('FRB121102')
   
   # Access basic properties
   print(f"Coordinates: {frb121102.coord}")
   print(f"Dispersion Measure: {frb121102.DM}")
   print(f"Error ellipse: {frb121102.eellipse}")

Working with Catalogues
~~~~~~~~~~~~~~~~~~~~~~~

Load and explore FRB catalogues:

.. code-block:: python

   from frb.frbcat import FRBCat
   
   # Load the FRB catalogue
   cat = FRBCat()
   
   # Basic catalogue information
   print(f"Number of FRBs: {len(cat.frbcat)}")
   print(f"DM range: {cat.frbcat['DM'].min():.1f} - {cat.frbcat['DM'].max():.1f} pc/cm³")
   
   # Filter high-DM FRBs
   high_dm_frbs = cat.frbcat[cat.frbcat['DM'] > 1000]
   print(f"High-DM FRBs (>1000): {len(high_dm_frbs)}")

Dispersion Measure Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate cosmic and IGM contributions to dispersion measure:

.. code-block:: python

   from frb.dm import igm
   from astropy import units as u
   
   # Calculate cosmic DM contribution at a given redshift
   z_frb = 0.5
   DM_cosmic = igm.DM_cosmic(z_frb)
   print(f"Cosmic DM at z={z_frb}: {DM_cosmic}")
   
   # Estimate redshift from observed DM
   DM_observed = 500 * u.pc / u.cm**3
   z_estimated = igm.z_from_DM(DM_observed)
   print(f"Estimated redshift for DM={DM_observed}: {z_estimated:.2f}")

Scattering Analysis
~~~~~~~~~~~~~~~~~~~

Analyze pulse scattering and broadening:

.. code-block:: python

   from frb import turb_scattering as ts
   from astropy import units as u
   
   # Set up scattering parameters
   n_e = 1e-3 * u.cm**(-3)     # Electron density
   nu_obs = 1.4 * u.GHz        # Observation frequency
   L = 50 * u.kpc              # Structure size
   R = 1 * u.pc                # Cloud size
   
   # Calculate scattering angle
   theta = ts.theta_mist(n_e, nu_obs, L=L, R=R)
   print(f"Scattering angle: {theta.to('microarcsec'):.2f}")
   
   # Calculate temporal broadening
   z_FRB = 1.0    # FRB redshift  
   z_lens = 0.5   # Lens redshift
   tau = ts.tau_mist(n_e, nu_obs, z_FRB, z_lens, L=L, R=R)
   print(f"Temporal broadening: {tau.to('ms'):.2f}")

Host Galaxy Analysis
~~~~~~~~~~~~~~~~~~~~

Analyze FRB host galaxies when available:

.. code-block:: python

   # Load FRB with known host
   frb180924 = ffrb.FRB.by_name('FRB180924')
   
   # Access host galaxy
   host = frb180924.grab_host()
   print(f"Host galaxy properties: {host.derived}")
   
   # Load spectral data if available
   try:
       meta, spec = host.get_metaspec()
       print(f"Spectrum loaded: {len(meta)} spectra available")
   except:
       print("No spectral data available")

Common Workflows
----------------

Workflow 1: Basic FRB Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Complete analysis of a single FRB:

.. code-block:: python

   import frb as ffrb
   from frb.dm import igm
   from astropy import units as u
   
   # Load FRB
   frb_name = 'FRB121102'
   frb_obj = ffrb.FRB.by_name(frb_name)
   
   print(f"=== Analysis of {frb_name} ===")
   print(f"Position: {frb_obj.coord}")
   print(f"Observed DM: {frb_obj.DM}")
   
   # Calculate expected cosmic DM at estimated redshift
   if hasattr(frb_obj, 'z') and frb_obj.z is not None:
       z_frb = frb_obj.z
       DM_cosmic_expected = igm.DM_cosmic(z_frb)
       print(f"Expected cosmic DM at z={z_frb:.2f}: {DM_cosmic_expected}")
       
       # Calculate excess DM
       DM_excess = frb_obj.DM - DM_cosmic_expected
       print(f"Excess DM (host+local): {DM_excess}")

Workflow 2: Population Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyze properties of the FRB population:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from frb.frbcat import FRBCat
   
   # Load catalogue
   cat = FRBCat()
   
   # Get DM values (remove invalid entries)
   dms = cat.frbcat['DM']
   valid_dms = dms[dms > 0]
   
   # Basic statistics
   print(f"DM Statistics:")
   print(f"  Mean: {np.mean(valid_dms):.1f} pc/cm³")
   print(f"  Median: {np.median(valid_dms):.1f} pc/cm³") 
   print(f"  Range: {np.min(valid_dms):.1f} - {np.max(valid_dms):.1f} pc/cm³")
   
   # Plot DM distribution
   plt.figure(figsize=(10, 6))
   plt.hist(valid_dms, bins=20, alpha=0.7, edgecolor='black')
   plt.xlabel('Dispersion Measure (pc/cm³)')
   plt.ylabel('Number of FRBs')
   plt.title('FRB Dispersion Measure Distribution')
   plt.grid(True, alpha=0.3)
   plt.show()

Workflow 3: Scattering Model Comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare different scattering models:

.. code-block:: python

   from frb import turb_scattering as ts
   from astropy import units as u
   import numpy as np
   
   # Set up parameter ranges
   frequencies = np.logspace(0, 2, 50) * u.MHz  # 1 MHz to 100 MHz
   n_e = 1e-3 * u.cm**(-3)
   
   # Calculate scattering for different structure sizes
   L_values = [10, 50, 100] * u.kpc
   
   scattering_angles = {}
   for L in L_values:
       angles = []
       for freq in frequencies:
           theta = ts.theta_mist(n_e, freq, L=L)
           angles.append(theta.to('microarcsec').value)
       scattering_angles[f'{L.value} kpc'] = angles
   
   # Plot results
   plt.figure(figsize=(10, 6))
   for label, angles in scattering_angles.items():
       plt.loglog(frequencies, angles, label=f'L = {label}')
   
   plt.xlabel('Frequency (MHz)')
   plt.ylabel('Scattering Angle (μas)')
   plt.title('Scattering Angle vs Frequency')
   plt.legend()
   plt.grid(True)
   plt.show()

Command Line Tools
------------------

Galaxy Search Tool
~~~~~~~~~~~~~~~~~~

The package includes command-line tools for common tasks:

.. code-block:: bash

   # Search for galaxies near an FRB position
   frb_galaxies J081240.7+320809 --rho 300
   
   # Search by FRB name with plotting
   frb_galaxies FRB180924 --plot
   
   # Specify angular offset instead of physical distance
   frb_galaxies "07:45:00.47,34:17:31.1" --ang_offset 30

Common Parameters
~~~~~~~~~~~~~~~~~

- ``--rho RHO``: Maximum impact parameter in kpc (default: 300)
- ``--ang_offset ANG_OFFSET``: Maximum offset in arcsec (overrides --rho)
- ``--cat``: Only show data from the catalog (not meta)
- ``--specdb SPECDB``: Specify specDB file path
- ``-p, --plot``: Launch a plotting GUI

Working with Your Own Data
--------------------------

Adding Custom FRBs
~~~~~~~~~~~~~~~~~~

You can work with your own FRB data by creating FRB objects:

.. code-block:: python

   from astropy.coordinates import SkyCoord
   from astropy import units as u
   
   # Create a custom FRB object (if supported by the API)
   coord = SkyCoord(ra=123.45*u.deg, dec=-23.67*u.deg)
   DM_value = 750 * u.pc / u.cm**3
   
   # Use the analysis tools on your data
   z_est = igm.z_from_DM(DM_value, coord=coord)
   print(f"Estimated redshift: {z_est:.2f}")

Importing External Catalogues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from frb.frbcat import FRBCat
   
   # Load custom catalogue file
   custom_cat = FRBCat(frbcat_file='path/to/your/catalogue.csv')
   
   # Process the data
   print(f"Loaded {len(custom_cat.frbcat)} FRBs from custom catalogue")

Next Steps
----------

Now that you're familiar with the basics:

1. **Explore Advanced Features**: Check out the :doc:`tutorials` for in-depth examples
2. **Read the API Documentation**: Browse :doc:`api/index` for complete function references  
3. **Study Real Examples**: Look at :doc:`examples` for practical applications
4. **Contribute**: See :doc:`contributing` if you want to help improve the package

Common Gotchas
--------------

**Units**: Always pay attention to astropy units. The package uses physical units throughout:

.. code-block:: python

   # Correct
   DM = 500 * u.pc / u.cm**3
   
   # Incorrect - will cause errors
   DM = 500  # No units

**Data Availability**: Not all FRBs have complete data. Always check for ``None`` values:

.. code-block:: python

   frb_obj = ffrb.FRB.by_name('SomeFRB')
   if frb_obj.z is not None:
       print(f"Redshift: {frb_obj.z}")
   else:
       print("No redshift available")

**Coordinate Systems**: Be aware of coordinate system conventions:

.. code-block:: python

   # The package handles coordinate conversions automatically
   print(f"Galactic coordinates: {frb_obj.coord.galactic}")
   print(f"Equatorial coordinates: {frb_obj.coord.icrs}")

Getting Help
------------

- **Documentation**: This documentation covers most use cases
- **GitHub Issues**: https://github.com/FRBs/FRB/issues
- **Community**: Connect with other users through the FRBs organization
- **Examples**: Check the ``docs/nb/`` directory for Jupyter notebook examples