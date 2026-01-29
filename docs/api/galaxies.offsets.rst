Galaxy Offsets
==============

.. automodule:: frb.galaxies.offsets
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module provides functions for calculating positional offsets between 
FRBs and their potential host galaxies. It handles both angular and physical
offset measurements, accounting for localization uncertainties and coordinate
transformations.

Functions
---------

Primary Offset Functions
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.offsets.angular_offset

   Calculate angular offset between FRB position and galaxy coordinates.
   
   This is the primary function for computing separations, handling:
   
   * FRB localization error ellipses
   * Galaxy position uncertainties  
   * Statistical error propagation
   * Multiple offset definitions (best estimate vs. averaged)

.. autofunction:: frb.galaxies.offsets.physical_offset

   Convert angular offsets to physical separations using cosmological distances.
   
   Takes angular separations and converts to proper physical distances 
   in kpc, accounting for the galaxy redshift and assumed cosmology.

Error Analysis Functions  
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.offsets.offset_uncertainty

   Calculate uncertainties in offset measurements.
   
   Propagates errors from:
   
   * FRB localization uncertainty ellipse
   * Galaxy astrometric errors
   * Systematic coordinate uncertainties
   * Statistical measurement errors

.. autofunction:: frb.galaxies.offsets.deproject_offset

   Deproject observed offsets to account for galaxy inclination.
   
   For edge-on or highly inclined galaxies, converts sky-plane offsets
   to deprojected separations within the galaxy disk.

Coordinate System Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.offsets.galactic_coords

   Transform offsets to galaxy-centric coordinate system.
   
   Rotates offset vectors to align with galaxy major axis, useful
   for studying offset distributions relative to galaxy structure.

.. autofunction:: frb.galaxies.offsets.position_angle

   Calculate position angle of FRB relative to galaxy center.
   
   Returns the position angle (East of North) of the FRB location
   relative to the galaxy centroid.

Statistical Functions
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.offsets.offset_probability

   Calculate probability of chance alignment given offset distribution.
   
   Uses offset measurements and galaxy number density to assess the
   likelihood that an apparent FRB-galaxy association is coincidental.

.. autofunction:: frb.galaxies.offsets.compare_offset_distributions

   Compare offset distributions between different galaxy populations.
   
   Statistical comparison of FRB offset distributions for different
   host galaxy types, redshift ranges, or other sample cuts.

Offset Types and Definitions
----------------------------

The module implements several offset definitions:

**Angular Offsets**
   * `ang_best`: Offset from FRB localization centroid to galaxy center
   * `ang_avg`: Offset averaged over FRB localization probability distribution  
   * Angular offsets reported in arcseconds

**Physical Offsets**  
   * Proper physical distance in kpc at galaxy redshift
   * Corrected for cosmological expansion
   * Uses `ang_best` by default unless specified

**Deprojected Offsets**
   * Corrected for galaxy inclination angle
   * Represents true separation within galaxy disk
   * Requires morphological information

Error Propagation
-----------------

Comprehensive error handling includes:

**FRB Localization Errors**
   - Error ellipse semi-major and semi-minor axes
   - Position angle of error ellipse
   - Confidence level specification

**Galaxy Position Errors**
   - Astrometric tie uncertainties  
   - Source extraction errors
   - Proper motion corrections (for nearby galaxies)

**Systematic Uncertainties**
   - Absolute astrometric calibration
   - Reference frame differences
   - Coordinate epoch corrections

Examples
--------

Basic offset calculation:

.. code-block:: python

    from frb.galaxies import offsets
    from frb.frb import FRB
    from frb.galaxies.frbgalaxy import FRBGalaxy
    
    # Create FRB and galaxy objects
    frb = FRB.by_name('FRB180924')
    galaxy = FRBGalaxy(ra=349.24, dec=-40.9, frb=frb)
    
    # Calculate angular offsets (done automatically in FRBGalaxy.__init__)
    ang_avg, ang_avg_err, ang_best, ang_best_err = offsets.angular_offset(frb, galaxy)
    
    print(f"Angular offset (best): {ang_best:.2f} ± {ang_best_err:.2f} arcsec")
    print(f"Angular offset (averaged): {ang_avg:.2f} ± {ang_avg_err:.2f} arcsec")

Physical offset calculation:

.. code-block:: python

    # Set galaxy redshift first
    galaxy.set_z(0.3214, 'spec', err=0.0001)
    
    # Calculate physical offset
    phys_offset = offsets.physical_offset(
        ang_best,  # angular offset in arcsec
        galaxy.z,  # redshift
        cosmo=galaxy.cosmo
    )
    
    print(f"Physical offset: {phys_offset:.1f} kpc")

Including position errors:

.. code-block:: python

    # Set galaxy position uncertainties
    galaxy.positional_error['ra_astrometric'] = 0.1  # arcsec
    galaxy.positional_error['dec_astrometric'] = 0.1  # arcsec
    galaxy.positional_error['ra_source'] = 0.05  # arcsec  
    galaxy.positional_error['dec_source'] = 0.05  # arcsec
    
    # Recalculate with position errors included
    ang_avg, ang_avg_err, ang_best, ang_best_err = offsets.angular_offset(frb, galaxy)
    
    print(f"Offset with position errors: {ang_best:.2f} ± {ang_best_err:.2f} arcsec")

Position angle calculation:

.. code-block:: python

    # Calculate position angle of FRB relative to galaxy
    pa = offsets.position_angle(frb.coord, galaxy.coord)
    print(f"Position angle: {pa:.1f} degrees East of North")

Probability assessment:

.. code-block:: python

    # Assess chance alignment probability
    # (requires galaxy surface density information)
    prob_chance = offsets.offset_probability(
        ang_best,  # observed offset
        galaxy_density=1000,  # galaxies per sq. arcmin  
        magnitude_limit=25.0   # survey depth
    )
    
    print(f"Chance alignment probability: {prob_chance:.3f}")

Working with morphology:

.. code-block:: python

    # For galaxies with inclination information
    if 'inclination' in galaxy.morphology:
        # Calculate deprojected offset
        deprojected = offsets.deproject_offset(
            ang_best,
            galaxy.morphology['inclination'],
            galaxy.morphology['position_angle'],
            frb_pa=pa
        )
        
        print(f"Deprojected offset: {deprojected:.2f} arcsec")

Bulk analysis:

.. code-block:: python

    from frb.galaxies.utils import list_of_hosts
    
    # Get all hosts and calculate offset distribution
    frbs, hosts = list_of_hosts()
    
    offsets_list = []
    for host in hosts:
        if host.z is not None:
            phys_off = offsets.physical_offset(
                host.offsets['ang_best'],
                host.z,
                cosmo=host.cosmo
            )
            offsets_list.append(phys_off)
    
    import numpy as np
    median_offset = np.median(offsets_list)
    print(f"Median host offset: {median_offset:.1f} kpc")