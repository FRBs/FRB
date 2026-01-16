Nebular Emission Line Analysis
==============================

.. automodule:: frb.galaxies.nebular
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module provides functions for analyzing nebular emission lines in galaxy spectra,
including line flux measurements, extinction corrections, and star formation rate
calculations from emission line diagnostics.

Functions
---------

Line Analysis Functions
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.nebular.calc_lum

   Calculate emission line luminosity with optional dust extinction correction.
   
   This function converts observed line fluxes to intrinsic luminosities,
   applying distance corrections and optionally correcting for dust extinction
   using the Balmer decrement or other extinction indicators.

.. autofunction:: frb.galaxies.nebular.measure_lines

   Measure emission line fluxes from galaxy spectra.
   
   Automated line fitting routine that identifies and measures common
   nebular emission lines including:
   
   * Hydrogen Balmer series (Hα, Hβ, Hγ, Hδ)
   * Oxygen lines ([OII] λ3727, [OIII] λλ4959,5007)  
   * Nitrogen lines ([NII] λλ6548,6583)
   * Sulfur lines ([SII] λλ6717,6731)

Star Formation Rate Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.nebular.sfr_ha

   Calculate star formation rate from Hα luminosity.
   
   Uses the calibration from Kennicutt (1998) with appropriate corrections
   for dust extinction and metallicity effects.

.. autofunction:: frb.galaxies.nebular.sfr_oii

   Calculate star formation rate from [OII] λ3727 luminosity.
   
   Useful for higher redshift galaxies where Hα is redshifted out of
   optical wavelength range.

Extinction and Reddening
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.nebular.calc_extinction

   Calculate dust extinction from Balmer line ratios.
   
   Uses the Balmer decrement (Hα/Hβ ratio) to determine the dust
   extinction affecting nebular emission, assuming case B recombination.

.. autofunction:: frb.galaxies.nebular.get_ebv

   Get Galactic extinction E(B-V) for given coordinates.
   
   Queries dust maps to obtain Milky Way foreground extinction values
   using Schlegel, Finkbeiner & Davis (1998) or other dust maps.

Diagnostic Functions
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.nebular.bpt_classification

   Classify galaxies using BPT (Baldwin, Phillips & Terlevich) diagnostics.
   
   Uses emission line ratios to distinguish between:
   
   * Star-forming regions
   * Active galactic nuclei (AGN) 
   * Low-ionization nuclear emission regions (LINERs)
   * Composite systems

.. autofunction:: frb.galaxies.nebular.metallicity_diagnostics

   Calculate gas-phase metallicity from emission line ratios.
   
   Implements various metallicity calibrations:
   
   * N2 method ([NII]/Hα)
   * O3N2 method ([OIII]/Hβ vs [NII]/Hα)
   * R23 method (([OII]+[OIII])/Hβ)

Physical Properties
~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.nebular.electron_density

   Calculate electron density from [SII] doublet ratio.
   
   Uses the [SII] λ6717/λ6731 ratio to determine the electron density
   in HII regions, sensitive in the range ~10-10^4 cm^-3.

.. autofunction:: frb.galaxies.nebular.ionization_parameter

   Estimate ionization parameter from emission line diagnostics.
   
   Uses various line ratio diagnostics to constrain the ionization
   parameter in HII regions and AGN narrow-line regions.

Constants and Calibrations
--------------------------

The module includes various physical constants and calibration factors:

**Recombination Constants**
   - Case B recombination coefficients
   - Temperature and density dependent line ratios
   - Intrinsic Balmer line ratios

**Extinction Laws**
   - Cardelli, Clayton & Mathis (1989) extinction curve
   - Calzetti et al. (2000) starburst attenuation law
   - Fitzpatrick (1999) Milky Way extinction

**SFR Calibrations**
   - Kennicutt (1998) Hα-SFR relation
   - Modern IMF-corrected calibrations
   - Metallicity-dependent corrections

Examples
--------

Basic line analysis:

.. code-block:: python

    from frb.galaxies import nebular
    from frb.galaxies.frbgalaxy import FRBGalaxy
    
    # Assuming galaxy object with emission line measurements
    galaxy = FRBGalaxy(ra=180.0, dec=45.0, frb=frb_object)
    
    # Set emission line fluxes (in units of erg/s/cm^2)
    galaxy.neb_lines['Ha_flux'] = 5.2e-16
    galaxy.neb_lines['Ha_flux_err'] = 0.3e-16
    galaxy.neb_lines['Hb_flux'] = 1.1e-16
    galaxy.neb_lines['OIII_5007_flux'] = 2.8e-16
    
    # Calculate luminosity
    ha_lum = nebular.calc_lum(galaxy, 'Ha')
    print(f"Hα luminosity: {ha_lum:.2e} erg/s")

Star formation rate calculation:

.. code-block:: python

    # Calculate extinction from Balmer decrement
    extinction = nebular.calc_extinction(
        galaxy.neb_lines['Ha_flux'],
        galaxy.neb_lines['Hb_flux']
    )
    print(f"A_V = {extinction:.2f} mag")
    
    # Calculate extinction-corrected SFR
    sfr = nebular.sfr_ha(ha_lum, extinction=extinction)
    print(f"Star formation rate: {sfr:.2f} Msun/yr")

BPT classification:

.. code-block:: python

    # Calculate line ratios for BPT diagram
    line_ratios = {
        'NII_Ha': galaxy.neb_lines['NII_6583_flux'] / galaxy.neb_lines['Ha_flux'],
        'OIII_Hb': galaxy.neb_lines['OIII_5007_flux'] / galaxy.neb_lines['Hb_flux'],
        'SII_Ha': galaxy.neb_lines['SII_6717_flux'] / galaxy.neb_lines['Ha_flux'],
        'OI_Ha': galaxy.neb_lines['OI_6300_flux'] / galaxy.neb_lines['Ha_flux']
    }
    
    # Classify source type
    classification = nebular.bpt_classification(line_ratios)
    print(f"BPT classification: {classification}")

Metallicity analysis:

.. code-block:: python

    # Calculate metallicity using N2 method
    n2_ratio = galaxy.neb_lines['NII_6583_flux'] / galaxy.neb_lines['Ha_flux']
    metallicity_n2 = nebular.metallicity_diagnostics(n2_ratio, method='N2')
    
    # Convert to 12 + log(O/H) scale
    oh_abundance = 8.69 + metallicity_n2
    print(f"12 + log(O/H) = {oh_abundance:.2f}")

Electron density measurement:

.. code-block:: python

    # Calculate electron density from [SII] doublet
    sii_ratio = (galaxy.neb_lines['SII_6717_flux'] / 
                 galaxy.neb_lines['SII_6731_flux'])
    
    n_e = nebular.electron_density(sii_ratio)
    print(f"Electron density: {n_e:.1f} cm^-3")

Galactic extinction correction:

.. code-block:: python

    from astropy.coordinates import SkyCoord
    
    # Get Galactic extinction
    coord = SkyCoord(ra=180.0*u.deg, dec=45.0*u.deg)
    ebv_gal = nebular.get_ebv(coord)
    
    # Apply foreground correction to observed fluxes
    extinction_corr = 10**(0.4 * 2.5 * ebv_gal)  # Hα extinction
    intrinsic_flux = galaxy.neb_lines['Ha_flux'] * extinction_corr
    
    print(f"Galactic E(B-V): {ebv_gal:.3f}")
    print(f"Corrected Hα flux: {intrinsic_flux:.2e} erg/s/cm^2")