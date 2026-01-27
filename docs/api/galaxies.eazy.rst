EAZY
====

.. automodule:: frb.galaxies.eazy
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module facilitates scripting of EAZY (Easy and Accurate Z from Yale) 
photometric redshift analysis. It provides tools to set up EAZY runs, 
generate input files, and process photometric redshift estimates.

.. note::
   This module requires EAZY to be installed and the EAZYDIR environment 
   variable to be properly set.

Constants and Configuration
---------------------------

Filter Mapping
~~~~~~~~~~~~~~

.. autodata:: frb.galaxies.eazy.frb_to_eazy_filters

   Dictionary mapping FRB filter names to EAZY filter indices. Includes filters from:
   
   * DECaLS/Legacy Survey (g, r, z)
   * DES (u, g, r, i, z, Y) 
   * SDSS (u, g, r, i, z)
   * WISE (W1, W2, W3, W4)
   * Pan-STARRS (g, r, i, z, y)
   * VISTA (Y, J, H, Ks)
   * Various ground-based instruments

Template Sets
~~~~~~~~~~~~~

.. autodata:: frb.galaxies.eazy._template_list

   Available EAZY template sets:
   ('br07_default', 'br07_goods', 'cww+kin', 'eazy_v1.0', 'eazy_v1.1_lines', 
    'eazy_v1.2_dusty', 'eazy_v1.3', 'pegase', 'pegase13')

Prior Options
~~~~~~~~~~~~~

.. autodata:: frb.galaxies.eazy._acceptable_priors

   Available redshift priors: ('prior_R_zmax7', 'prior_K_zmax7', 
   'prior_R_extend', 'prior_K_extend')

Functions
---------

Setup Functions
~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.eazy.eazy_setup

   Set up EAZY input directory with required templates and filter files.

.. autofunction:: frb.galaxies.eazy.eazy_filenames

   Generate standardized filenames for EAZY input files.

Input File Generation
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.eazy.eazy_input_files

   Write complete set of input files needed to run EAZY analysis.
   
   This function creates:
   
   * Catalog file with photometric measurements
   * Translation file mapping columns to EAZY format  
   * Parameter file with analysis configuration

Analysis Functions
~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.eazy.eazy_photoz

   Run complete EAZY photometric redshift analysis.

.. autofunction:: frb.galaxies.eazy.eazy_cat_from_frb_photom

   Convert FRB galaxy photometry to EAZY catalog format.

File I/O Functions
~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.eazy.eazy_getpz

   Read EAZY photometric redshift results from output files.

Configuration Parameters
------------------------

Key parameters for EAZY analysis include:

**Redshift Grid**
   - `zmin`: Minimum redshift (default: 0.050)
   - `zmax`: Maximum redshift (default: 7.000) 
   - `zstep`: Redshift step size (default: 0.001)

**Template Settings**
   - `templates`: Template set to use (default: 'eazy_v1.3')
   - `combo`: Template combination mode (1, 2, 99, -1, 'a')

**Prior Configuration**
   - `prior`: Prior file name (default: 'prior_R_zmax7')
   - `prior_filter`: Filter to use for magnitude prior
   - `prior_ABZP`: AB magnitude zero-point for prior (default: 23.9)

**Quality Control**
   - `n_min_col`: Minimum number of filter detections required
   - `magnitudes`: Use magnitudes instead of fluxes as input

Examples
--------

Basic EAZY setup and analysis:

.. code-block:: python

    from frb.galaxies import eazy
    from frb.galaxies.frbgalaxy import FRBGalaxy
    
    # Set up EAZY working directory
    eazy.eazy_setup('eazy_work/')
    
    # Generate input files from galaxy photometry  
    eazy.eazy_input_files(
        galaxy.photom,
        input_dir='eazy_work/',
        name='FRB180924_host',
        out_dir='output/',
        templates='eazy_v1.3',
        zmax=4.0
    )

Running photometric redshift analysis:

.. code-block:: python

    # Run complete EAZY analysis
    results = eazy.eazy_photoz(
        'galaxy_catalog.fits',
        input_dir='eazy_work/',
        name='survey_field',
        zmax=6.0,
        templates='eazy_v1.2_dusty',
        prior='prior_K_zmax7'
    )

Custom configuration:

.. code-block:: python

    # Advanced configuration with custom parameters
    eazy.eazy_input_files(
        photom_dict,
        input_dir='analysis/',
        name='high_z_candidate', 
        out_dir='results/',
        templates='eazy_v1.1_lines',  # Line emission templates
        combo='a',  # All template combinations
        zmin=0.1,
        zmax=8.0,
        zstep=0.002,
        prior_filter='i',  # Use i-band for prior
        n_min_col=4  # Require 4+ band detections
    )

Integration with FRBGalaxy:

.. code-block:: python

    from frb.frb import FRB
    from frb.galaxies.frbgalaxy import FRBGalaxy
    
    # Create galaxy object with photometry
    frb = FRB.by_name('FRB121102')  
    galaxy = FRBGalaxy(ra=82.998, dec=33.148, frb=frb)
    
    # Populate with photometric measurements
    galaxy.parse_photom(photom_table)
    
    # Run EAZY analysis
    eazy_results = eazy.eazy_photoz(
        galaxy.photom,
        input_dir='frb121102_eazy/',
        name='host_galaxy'
    )

Output Processing:

.. code-block:: python

    # Read EAZY results
    zout = eazy.eazy_getpz('OUTPUT/photz.zout.fits')
    
    # Extract best redshift estimates
    zbest = zout['z_peak']
    z_err_lo = zout['z_err_lo'] 
    z_err_hi = zout['z_err_hi']
    
    print(f"Photo-z: {zbest[0]:.3f} +{z_err_hi[0]:.3f} -{z_err_lo[0]:.3f}")