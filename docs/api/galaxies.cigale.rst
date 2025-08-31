CIGALE
======

.. automodule:: frb.galaxies.cigale
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module provides automation for CIGALE (Code Investigating GALaxy Emission) 
spectral energy distribution fitting. It generates configuration files and runs 
the standard pcigale script for single galaxy analysis.

.. note::
   This module requires pcigale to be installed on the system.

Constants
---------

.. autodata:: frb.galaxies.cigale._DEFAULT_SED_MODULES
   
   Default list of SED modules for CIGALE analysis:
   ('sfhdelayed', 'bc03', 'nebular', 'dustatt_calzleit', 'dale2014', 
    'restframe_parameters', 'redshifting')

Functions
---------

Main Functions
~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.cigale.run

   Input parameters and run CIGALE analysis for a photometry table.
   
   This is the main entry point for running CIGALE on a table of photometric
   measurements. It handles both single galaxy and multi-galaxy analysis.

.. autofunction:: frb.galaxies.cigale.gen_cigale_in

   Generate input data file for CIGALE from photometric measurements.

.. autofunction:: frb.galaxies.cigale._initialise

   Initialize CIGALE configuration with specified parameters.

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.cigale._sed_default_params

   Set the default parameters for CIGALE SED modules.
   
   Provides default parameter grids for different SED modules including:
   
   * **sfhdelayed**: Delayed star formation history
   * **bc03**: Bruzual & Charlot 2003 stellar population models  
   * **nebular**: Nebular emission modeling
   * **dustatt_calzleit**: Calzetti attenuation law
   * **dale2014**: Dust emission templates

Parameters
----------

The module supports extensive customization of CIGALE analysis through various parameters:

**SED Modules**
   - Star formation history models (sfhdelayed, exponential, etc.)
   - Stellar population synthesis (bc03, m05, etc.) 
   - Nebular emission (nebular)
   - Dust attenuation (dustatt_calzleit, dustatt_modified_starburst, etc.)
   - Dust emission (dale2014, casey2012, etc.)

**Analysis Parameters**
   - Redshift estimation (photometric vs spectroscopic)
   - Output variables (stellar mass, SFR, metallicity, etc.)
   - Core usage and computational settings

Examples
--------

Basic CIGALE analysis:

.. code-block:: python

    from frb.galaxies import cigale
    from astropy.table import Table
    
    # Load photometry table
    photom_table = Table.read('galaxy_photom.fits')
    
    # Run CIGALE with default settings
    cigale.run(photom_table, zcol='z_spec', 
               data_file='input.fits',
               config_file='config.ini',
               plot=True)

Custom SED modules:

.. code-block:: python

    # Define custom SED modules
    custom_modules = ['sfhdelayed', 'bc03', 'nebular', 'dustatt_calzleit']
    
    # Run with custom configuration
    cigale.run(photom_table, zcol='z_phot',
               sed_modules=custom_modules,
               save_sed=True,
               cores=4)

Integration with FRBGalaxy:

.. code-block:: python

    from frb.galaxies.frbgalaxy import FRBGalaxy
    
    # Assuming galaxy object exists with photometry
    galaxy.run_cigale(data_file='frb_galaxy.fits',
                      wait_for_input=False,
                      plot=True)