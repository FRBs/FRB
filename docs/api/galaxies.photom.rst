Photometry 
==========

.. automodule:: frb.galaxies.photom
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module provides functions for photometric analysis of FRB host galaxies,
including magnitude-to-flux conversions, aperture photometry corrections,
and integration with various survey photometric systems.

Functions
---------

Photometric Conversions
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.photom.mag_to_flux

   Convert magnitudes to flux densities with proper error propagation.
   
   Handles conversion between AB magnitude system and flux densities
   in various units (mJy, μJy, erg/s/cm²/Hz, etc.).

.. autofunction:: frb.galaxies.photom.flux_to_mag

   Convert flux densities to AB magnitudes with error propagation.

.. autofunction:: frb.galaxies.photom.extinction_correct

   Apply extinction corrections to photometric measurements.
   
   Uses various extinction laws:
   
   * Cardelli, Clayton & Mathis (1989) - Milky Way
   * Calzetti et al. (2000) - Starburst galaxies
   * Fitzpatrick (1999) - Updated Milky Way
   * Gordon et al. (2003) - Small Magellanic Cloud

Aperture Corrections
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.photom.aperture_correction

   Apply aperture corrections to photometric measurements.
   
   Corrects photometry measured in fixed apertures to total
   magnitudes using growth curve analysis or model fitting.

.. autofunction:: frb.galaxies.photom.psf_correction

   Correct point source contamination in galaxy photometry.
   
   Removes or accounts for foreground stars or AGN point source
   contributions to integrated galaxy photometry.

Survey Integration
~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.photom.match_survey_photometry

   Cross-match and combine photometry from multiple surveys.
   
   Handles systematic offsets between surveys and provides
   combined photometric datasets with proper error handling.

.. autofunction:: frb.galaxies.photom.synthetic_photometry

   Calculate synthetic photometry from spectra or SED models.
   
   Convolves input spectra with survey filter response functions
   to predict magnitudes in any photometric system.

Color and SED Analysis
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.photom.calculate_colors

   Calculate photometric colors with error propagation.
   
   Computes standard color indices (u-g, g-r, r-i, etc.) and
   handles cases with non-detections or upper limits.

.. autofunction:: frb.galaxies.photom.color_corrections

   Apply K-corrections and evolutionary corrections to photometry.
   
   Corrects observed photometry to rest-frame values accounting
   for redshift effects and cosmological evolution.

.. autofunction:: frb.galaxies.photom.sed_chi_squared

   Calculate chi-squared goodness of fit for SED models.
   
   Compares observed photometry with model predictions,
   properly handling upper limits and systematic uncertainties.

Quality Assessment
~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.photom.photom_quality_flags

   Generate quality flags for photometric measurements.
   
   Identifies potential issues:
   
   * Saturation effects
   * Contamination by nearby sources
   * Poor photometric conditions
   * Systematic calibration problems

.. autofunction:: frb.galaxies.photom.detect_outliers

   Detect outlier photometric measurements using statistical tests.
   
   Flags measurements that are inconsistent with SED expectations
   or show systematic deviations from neighboring bands.

Examples
--------

Basic magnitude-flux conversions:

.. code-block:: python

    from frb.galaxies import photom
    import numpy as np
    
    # Convert AB magnitude to flux density in mJy
    mag = 22.5
    mag_err = 0.1
    
    flux, flux_err = photom.mag_to_flux(mag, mag_err, units='mJy')
    print(f"Flux: {flux:.2f} ± {flux_err:.2f} mJy")
    
    # Convert back to magnitude
    mag_check, mag_err_check = photom.flux_to_mag(flux, flux_err)
    print(f"Magnitude: {mag_check:.2f} ± {mag_err_check:.3f}")

Extinction corrections:

.. code-block:: python

    # Apply Galactic extinction correction
    ebv_gal = 0.05  # E(B-V) from dust maps
    
    # Correct r-band magnitude
    r_obs = 22.8
    r_corr = photom.extinction_correct(
        r_obs, 
        ebv_gal, 
        filter_name='r',
        extinction_law='ccm89',
        rv=3.1
    )
    
    print(f"Observed r: {r_obs:.2f}")
    print(f"Corrected r: {r_corr:.2f}")
    print(f"Correction: {r_corr - r_obs:.3f} mag")

Color calculations:

.. code-block:: python

    # Calculate colors from photometry dictionary
    photom_dict = {
        'DES_g': 23.1, 'DES_g_err': 0.05,
        'DES_r': 22.3, 'DES_r_err': 0.03,  
        'DES_i': 21.9, 'DES_i_err': 0.04
    }
    
    # Calculate g-r color
    gr_color, gr_err = photom.calculate_colors(
        photom_dict, 'DES_g', 'DES_r'
    )
    
    # Calculate r-i color  
    ri_color, ri_err = photom.calculate_colors(
        photom_dict, 'DES_r', 'DES_i'
    )
    
    print(f"g-r = {gr_color:.2f} ± {gr_err:.3f}")
    print(f"r-i = {ri_color:.2f} ± {ri_err:.3f}")

Working with galaxy objects:

.. code-block:: python

    from frb.galaxies.frbgalaxy import FRBGalaxy
    
    # Assuming galaxy object with photometry loaded
    galaxy = FRBGalaxy(ra=180.0, dec=45.0, frb=frb_object)
    
    # Calculate synthetic V-band magnitude from available photometry
    v_synth = photom.synthetic_photometry(
        galaxy.photom, 
        target_filter='V',
        method='interpolation'
    )
    
    print(f"Synthetic V magnitude: {v_synth:.2f}")

SED fitting preparation:

.. code-block:: python

    # Prepare photometry for SED fitting
    clean_photom = photom.detect_outliers(galaxy.photom)
    
    # Apply quality flags
    quality_flags = photom.photom_quality_flags(galaxy.photom)
    
    # Remove flagged measurements
    sed_photom = {}
    for filt, mag in clean_photom.items():
        if quality_flags.get(filt, 0) == 0:  # Good quality
            sed_photom[filt] = mag
            
    print(f"Clean photometry: {len(sed_photom)} measurements")
    print(f"Rejected: {len(galaxy.photom) - len(sed_photom)} measurements")

Multi-survey combination:

.. code-block:: python

    # Combine photometry from multiple surveys
    survey_data = {
        'DES': {'g': 22.1, 'r': 21.5, 'i': 21.2},
        'SDSS': {'g': 22.0, 'r': 21.4, 'i': 21.1},  
        'Pan-STARRS': {'g': 22.05, 'r': 21.45, 'i': 21.15}
    }
    
    combined_photom = photom.match_survey_photometry(
        survey_data,
        weight_by_error=True,
        apply_systematic_corrections=True
    )
    
    print("Combined photometry:")
    for filt, data in combined_photom.items():
        print(f"  {filt}: {data['mag']:.2f} ± {data['err']:.3f}")