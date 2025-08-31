Photometric Redshifts
=====================

.. automodule:: frb.halos.photoz

Module for photometric redshift-based halo DM calculations and galaxy analysis.

This module combines photometric redshift estimates with halo models to calculate dispersion measure contributions from galaxy halos along FRB sightlines.

Constants
---------

.. autodata:: DEFAULT_DATA_FOLDER

   Default data folder path: "data"

Functions
---------

get_des_data
~~~~~~~~~~~~

.. autofunction:: get_des_data

   Download photometry for galaxies within an FRB field.

   **Parameters:**
   
   * ``coords`` : SkyCoord - Center coordinates for cone search
   * ``radius`` : Quantity, optional - Search radius (default: 15 arcmin)
   * ``starbright`` : float, optional - Lower r-band magnitude limit (default: 17)
   * ``starflagval`` : float, optional - Star classification upper limit (default: 0.9)
   * ``gaiacat`` : str, optional - Gaia catalog file for star removal
   * ``write`` : bool, optional - Write output table to file (default: False)
   * ``outfile`` : str, optional - Output filename

   **Returns:**
   
   * ``des_data`` : Table - DES galaxies within search radius

dm_grid
~~~~~~~

.. autofunction:: dm_grid

   Produce DM estimates for a 3D grid of redshift, offsets and halo masses.

   **Parameters:**
   
   * ``frb_z`` : float - FRB redshift
   * ``n_z`` : int, optional - Size of redshift grid (default: 100)
   * ``n_o`` : int, optional - Size of offset grid (default: 100) 
   * ``n_m`` : int, optional - Size of halo mass grid (default: 100)
   * ``max_log_mhalo`` : float, optional - Maximum log halo mass (default: 12.8)
   * ``outdir`` : str, optional - Output directory (default: DEFAULT_DATA_FOLDER)
   * ``outfile`` : str, optional - Output .npz filename

   Creates a 3D interpolation grid for DM calculations with dimensions:
   
   * Redshift: np.linspace(0, frb_z, n_z)
   * Offsets: np.linspace(0, 600, n_o) kpc
   * Halo masses: np.linspace(8, 16, n_m) log solar masses

mhalo_lookup_tables
~~~~~~~~~~~~~~~~~~~

.. autofunction:: mhalo_lookup_tables

   For each redshift in z_grid, produces files containing halo mass values corresponding to stellar masses.

   **Parameters:**
   
   * ``z_grid`` : list or ndarray - Redshift values to sample
   * ``datafolder`` : str, optional - Storage directory (default: DEFAULT_DATA_FOLDER)
   * ``n_cores`` : int, optional - CPU threads for parallel processing (default: 8)

   Values are produced by sampling the Moster+13 stellar-halo mass relation (SHMR).

_mhalo_lookup_table
~~~~~~~~~~~~~~~~~~~

.. autofunction:: _mhalo_lookup_table

   Internal function to create halo mass lookup tables for a single redshift.

   **Parameters:**
   
   * ``z`` : float - Redshift
   * ``npz_out`` : str, optional - Output .npz file path
   * ``n_cores`` : int, optional - CPU threads (default: 8)

   .. note::
      This is an internal function. Use mhalo_lookup_tables() directly if you know what you're doing.

_instantiate_intepolators
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: _instantiate_intepolators

   Produce interpolator functions for key quantities required for the analysis.

   **Parameters:**
   
   * ``datafolder`` : str, optional - Folder with interpolation data files (default: DEFAULT_DATA_FOLDER)
   * ``dmfilename`` : str, optional - DM interpolation data filename
   * ``frb_name`` : str, optional - FRB name (default: "FRB180924")

   **Returns:**
   
   * ``dm_interpolator`` : RegularGridInterpolator - DM(z, offset_kpc, log_mhalo)
   * ``mean_interp`` : interp2d - <log_mhalo(log_mstar, z)> based on SHMR
   * ``stddev_interp`` : interp2d - std.dev. log_mhalo(log_mstar, z) based on SHMR  
   * ``ang_dia_interp`` : interp1d - angular_diameter_distance(z) for default cosmology

_mhalo_realizations
~~~~~~~~~~~~~~~~~~~

.. autofunction:: _mhalo_realizations

   Generate halo mass realizations using lookup tables, accounting for both stellar mass uncertainty and SHMR scatter.

   **Parameters:**
   
   * ``log_mstar`` : float - log stellar mass in M_sun
   * ``log_mstar_err`` : float - log error in log_mstar
   * ``z`` : float - Redshift
   * ``mean_interp`` : interp2d - <log_mhalo(log_mstar, z)> interpolator
   * ``stddev_interp`` : interp2d - std.dev. log_mhalo(log_mstar, z) interpolator
   * ``n_mstar`` : int, optional - Number of stellar mass samples (default: 100)
   * ``n_norm`` : int, optional - Number of normal distribution samples (default: 10)
   * ``max_log_mhalo`` : float, optional - Maximum log halo mass (default: 12.8)

   **Returns:**
   
   * ``ndarray`` - Halo mass realizations

_dm_pdf
~~~~~~~

.. autofunction:: _dm_pdf

   Calculate DM realizations for a galaxy using photometric redshift and stellar mass estimates.

   **Parameters:**
   
   * ``cigale_tab`` : Table - CIGALE results for the galaxy
   * ``eazy_outdir`` : str - EAZY output directory
   * ``mean_interp`` : interp2d - Mean SHMR interpolator
   * ``stddev_interp`` : interp2d - SHMR scatter interpolator  
   * ``ang_dia_interp`` : interp1d - Angular diameter distance interpolator
   * ``dm_interpolator`` : RegularGridInterpolator - DM interpolator
   * ``n_cores`` : int, optional - CPU threads

   **Returns:**
   
   * ``dm_values`` : ndarray - DM realizations for the galaxy
   * ``z_draws`` : ndarray - Redshift draws used for DM calculations

full_analysis
~~~~~~~~~~~~~

.. autofunction:: full_analysis

   Perform complete photometric redshift-based halo DM analysis for an FRB field.

   **Parameters:**
   
   * ``frb`` : FRB - FRB object of interest
   * ``input_catfile`` : str - Input photometry catalog path (assumed DES format)
   * ``datafolder`` : str - Results storage directory
   * ``n_cores`` : int, optional - CPU threads (default: varies by function)
   * ``n_gals`` : int, optional - Limit analysis to n_gals galaxies for testing

   **Process:**
   
   1. Runs EAZY photometric redshift estimation
   2. Runs CIGALE SED fitting for stellar masses
   3. Creates interpolation grids
   4. Calculates DM realizations for all galaxies
   5. Saves results to compressed files

   **Outputs:**
   
   * DM_halos_final.npz - Sparse matrix of DM realizations
   * DM_halos_zdraws.npz - Redshift draws for each galaxy

Dependencies
------------

Required packages for full functionality:

* ``pathos`` - For multiprocessing: ``pip install pathos``
* ``progressbar2`` - For progress tracking: ``pip install progressbar2``  
* ``threedhst`` - For EAZY output: ``pip install threedhst``

Examples
--------

Basic DES data retrieval::

    from frb.halos.photoz import get_des_data
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    # Define search parameters
    center = SkyCoord('12h34m56s', '+12d34m56s')
    radius = 10 * u.arcmin
    
    # Get DES photometry
    des_cat = get_des_data(center, radius=radius, write=True)
    print(f"Found {len(des_cat)} galaxies")

Create DM interpolation grid::

    from frb.halos.photoz import dm_grid
    
    # Create 3D DM grid for FRB at z=0.5
    frb_redshift = 0.5
    dm_grid(frb_redshift, n_z=50, n_o=50, n_m=50,
           outdir='./halo_data/', outfile='dm_grid.npz')

Full analysis workflow::

    from frb.frb import FRB
    from frb.halos.photoz import full_analysis
    
    # Load FRB
    frb_obj = FRB.by_name('FRB20180924B')
    
    # Run complete halo DM analysis
    full_analysis(frb_obj, 
                  input_catfile='field_photometry.fits',
                  datafolder='./analysis_results/',
                  n_cores=8,
                  n_gals=100)  # Limit for testing

Setup interpolators::

    from frb.halos.photoz import _instantiate_intepolators
    
    # Create interpolation functions
    dm_interp, mean_interp, std_interp, ang_interp = _instantiate_intepolators(
        datafolder='./halo_data/',
        dmfilename='dm_grid.npz'
    )
    
    # Use interpolators for DM calculations
    import numpy as np
    z_test = 0.3
    offset_test = 100.0  # kpc
    mhalo_test = 12.0   # log solar masses
    
    dm_value = dm_interp((z_test, offset_test, mhalo_test))
    print(f"DM contribution: {dm_value:.2f} pc/cm³")

Generate halo mass lookup tables::

    from frb.halos.photoz import mhalo_lookup_tables
    import numpy as np
    
    # Create lookup tables for multiple redshifts
    z_array = np.linspace(0.1, 1.0, 10)
    mhalo_lookup_tables(z_array, 
                       datafolder='./lookup_tables/',
                       n_cores=4)

Advanced usage with custom parameters::

    from frb.halos.photoz import get_des_data, dm_grid
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    # Custom DES query with star removal
    coords = SkyCoord(ra=123.45, dec=-23.45, unit='deg')
    
    des_data = get_des_data(
        coords, 
        radius=20*u.arcmin,
        starbright=16.0,        # Remove bright stars
        starflagval=0.8,        # More aggressive star removal
        gaiacat='gaia_stars.csv', # Additional star catalog
        write=True,
        outfile='frb_field_photometry.fits'
    )
    
    # High-resolution DM grid
    dm_grid(frb_z=1.2, 
           n_z=200, n_o=150, n_m=120,  # Higher resolution
           max_log_mhalo=13.5,          # Include more massive halos
           outdir='./high_res_grid/',
           outfile='dm_grid_hires.npz')

Integration with FRB analysis::

    from frb.frb import FRB
    from frb.halos.photoz import get_des_data, full_analysis
    from astropy import units as u
    
    # Load FRB and get field data
    frb = FRB.by_name('FRB20180924B')
    
    # Get photometry around FRB location  
    field_data = get_des_data(frb.coord, radius=15*u.arcmin)
    
    # Save for analysis
    field_data.write('frb_field.fits', overwrite=True)
    
    # Run complete halo DM analysis
    full_analysis(frb, 'frb_field.fits', './results/', n_cores=8)
    
    print("Halo DM analysis complete!")
    print("Results saved in ./results/DM_halos_final.npz")