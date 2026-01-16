Galaxy Utilities
================

.. automodule:: frb.galaxies.utils
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module provides utility functions for working with FRB host galaxy data,
including database operations, table building, and various helper functions
for galaxy analysis.

Functions
---------

Database and Loading Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.utils.load_specdb

   Load spectroscopic database for galaxy analysis.

.. autofunction:: frb.galaxies.utils.list_of_hosts

   Generate a list of FRB host galaxies from the database.

.. autofunction:: frb.galaxies.utils.build_table_of_hosts

   Generate a Pandas table of FRB host galaxy data.
   
   This function extracts data from host objects and compiles it into a 
   comprehensive table including photometry, derived quantities, nebular
   line measurements, and morphological parameters.

Analysis Utilities  
~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.utils.load_f_mL

   Generate interpolator from magnitude to luminosity as function of redshift.
   
   Provides approximate magnitude-luminosity relationship up to z=4 for
   galaxy luminosity function analysis.

.. autofunction:: frb.galaxies.utils.load_PATH  

   Load up the PATH (Probabilistic Association of Transients with Hosts) table.

.. autofunction:: frb.galaxies.utils.deredden_spec

   Apply dereddening correction to galaxy spectra using dust extinction models.

Data Processing Functions
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.utils.build_photom_table

   Build standardized photometric table from various survey catalogs.

.. autofunction:: frb.galaxies.utils.parse_galfit_output

   Parse GALFIT morphological analysis output files.

.. autofunction:: frb.galaxies.utils.update_frbgalaxy_coords

   Update galaxy coordinates with improved astrometry.

Output and Table Management
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The module provides several functions for managing host galaxy data tables:

**Table Structure**
   - Coordinates stored as RA_host, DEC_host (degrees)
   - FRB names and objects included for cross-referencing
   - Units tracked in separate dictionary
   - Missing values handled with NaN

**Data Categories**
   - Photometric measurements across multiple surveys
   - Derived physical properties (mass, SFR, metallicity)
   - Nebular emission line fluxes and ratios
   - Morphological parameters from imaging analysis
   - Redshift measurements (spectroscopic and photometric)
   - Offset measurements from FRB positions

Examples
--------

Building host galaxy table:

.. code-block:: python

    from frb.galaxies import utils
    
    # Build comprehensive host table
    host_table, units_dict = utils.build_table_of_hosts(
        attrs=['derived', 'photom', 'neb_lines', 'morphology']
    )
    
    # Display table info
    print(f"Number of hosts: {len(host_table)}")
    print(f"Available columns: {list(host_table.columns)}")
    print(f"Units: {units_dict}")

Working with individual hosts:

.. code-block:: python

    # Get list of host objects
    frbs, hosts = utils.list_of_hosts(verbose=True)
    
    # Access specific host properties
    for host in hosts[:5]:  # First 5 hosts
        print(f"Host: {host.name}")
        print(f"  RA, Dec: {host.coord.ra.deg:.3f}, {host.coord.dec.deg:.3f}")
        print(f"  Redshift: {host.z}")
        if len(host.derived) > 0:
            if 'Mstar' in host.derived:
                print(f"  Stellar mass: {host.derived['Mstar']:.2e} Msun")

Loading and using spectroscopic database:

.. code-block:: python

    # Load spectroscopic database
    specdb = utils.load_specdb(specdb_file='custom_specdb.hdf5')
    
    if specdb is not None:
        # Query for spectra
        meta = specdb.meta_from_coords(coords, radius=5*u.arcsec)
        if len(meta) > 0:
            spectra = specdb.spectra_from_meta(meta)

Dereddening corrections:

.. code-block:: python

    from linetools.spectra.xspectrum1d import XSpectrum1D
    
    # Load spectrum
    spec = XSpectrum1D.from_file('galaxy_spectrum.fits')
    
    # Apply dereddening with AV = 0.5 mag
    corrected_spec = utils.deredden_spec(spec, AV=0.5)

PATH analysis integration:

.. code-block:: python

    # Load PATH results
    path_table = utils.load_PATH('adopted.csv')
    
    # Build host table with PATH probabilities
    host_table, units = utils.build_table_of_hosts()
    
    # PATH probabilities now included as P_Ox, P_O columns
    high_prob_hosts = host_table[host_table['P_Ox'] > 0.8]
    print(f"High probability associations: {len(high_prob_hosts)}")

Working with magnitude-luminosity relations:

.. code-block:: python

    # Load m-L interpolator  
    f_mL = utils.load_f_mL()
    
    # Get characteristic magnitude at different redshifts
    z_array = [0.1, 0.3, 0.5, 1.0, 2.0]
    m_star = f_mL(z_array)
    
    print("Characteristic magnitudes (r-band):")
    for z, m in zip(z_array, m_star):
        print(f"  z = {z:.1f}: m* = {m:.2f}")