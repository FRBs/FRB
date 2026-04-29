Galaxy Data Definitions
=======================

.. automodule:: frb.galaxies.defs
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module defines constants, valid field names, and data structures used 
throughout the FRB galaxies package. It serves as a central registry for
allowed values in galaxy data dictionaries and provides validation lists
for data integrity.

Constants and Valid Fields
--------------------------

Photometric Data
~~~~~~~~~~~~~~~~

.. autodata:: frb.galaxies.defs.valid_photom

   List of valid photometric filter names including:
   
   * **Survey filters**: DES (u,g,r,i,z,Y), SDSS (u,g,r,i,z), DECaLS (g,r,z)
   * **Near-IR**: VISTA (Y,J,H,Ks), 2MASS (J,H,K), WISE (W1,W2,W3,W4)  
   * **Space-based**: HST various filters, Spitzer IRAC channels
   * **Ground-based**: Pan-STARRS (g,r,i,z,y), LSST projected bands

.. autodata:: frb.galaxies.defs.valid_flux

   Corresponding flux measurements for each photometric band, in mJy units.

.. autodata:: frb.galaxies.defs.valid_filters

   Complete list of recognized filter systems across all surveys.

Derived Physical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: frb.galaxies.defs.valid_derived_photom

   Derived quantities from photometric SED fitting:
   
   * **Stellar properties**: Mstar, Mtotal, age_mass, Z_stellar
   * **Star formation**: SFR_photom, SFR_SED, lg_sSFR  
   * **Dust extinction**: EBV_photom, AV_young, AV_old
   * **AGN contribution**: f_AGN, agn_tau
   * **Rest-frame properties**: u-r, M_r, Lnu_r

.. autodata:: frb.galaxies.defs.valid_derived_nebular

   Properties derived from nebular emission line analysis:
   
   * **Extinction**: AV_nebular from Balmer decrement
   * **Star formation**: SFR_nebular from Hα, [OII] 
   * **Metallicity**: Z_gas from emission line ratios
   * **Excitation**: Ionization parameter, electron density

Nebular Emission Lines  
~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: frb.galaxies.defs.valid_neb_lines

   Recognized emission lines for nebular analysis:
   
   * **Hydrogen Balmer series**: Ha, Hb, Hg, Hd, H8, H9, H10, H11
   * **Oxygen lines**: [OII]_3727, [OIII]_4959, [OIII]_5007, [OI]_6300
   * **Nitrogen lines**: [NII]_6548, [NII]_6583  
   * **Sulfur lines**: [SII]_6717, [SII]_6731, [SIII]_6312
   * **Other ions**: [NeIII]_3869, [ArIII]_7136, HeI_5876, HeII_4686

Redshift Information
~~~~~~~~~~~~~~~~~~~~

.. autodata:: frb.galaxies.defs.valid_z

   Valid redshift measurement types:
   
   * **z**: Preferred redshift value
   * **z_spec**: Spectroscopic redshift  
   * **z_phot**: Photometric redshift
   * **z_SED**: SED fitting redshift
   * Error values: z_err, z_spec_err, z_phot_err

Morphological Properties
~~~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: frb.galaxies.defs.valid_morphology

   Morphological parameters from imaging analysis:
   
   * **Profile fitting**: n (Sersic index), Re (effective radius)
   * **Ellipticity**: ellip, position_angle, inclination  
   * **Surface brightness**: mu_e (effective surface brightness)
   * **Concentration**: C (concentration index), A (asymmetry)
   * **Multi-component**: bulge/disk decomposition parameters

Positional Information
~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: frb.galaxies.defs.valid_offsets

   Offset measurements between FRB and galaxy positions:
   
   * **ang_best**: Angular offset from localization centroid (arcsec)
   * **ang_avg**: Averaged angular offset over error distribution (arcsec)
   * **physical**: Physical offset in kpc (uses ang_best)

.. autodata:: frb.galaxies.defs.valid_positional_error

   Position uncertainty components:
   
   * **Astrometric errors**: ra_astrometric, dec_astrometric (arcsec)
   * **Source errors**: ra_source, dec_source (arcsec)
   * Includes systematic and random error contributions

Reference Information
~~~~~~~~~~~~~~~~~~~~~

.. autodata:: frb.galaxies.defs.valid_ref

   Reference tags for measurements, following ADS bibcode format.

.. autodata:: frb.galaxies.defs.valid_neb_ref

   References specific to nebular line measurements.

.. autodata:: frb.galaxies.defs.valid_derived_ref

   References for derived physical property calculations.

Survey-Specific Constants
-------------------------

.. autodata:: frb.galaxies.defs.DES_bands

   Dark Energy Survey filter bands: ['u', 'g', 'r', 'i', 'z', 'Y']

.. autodata:: frb.galaxies.defs.SDSS_bands  

   SDSS filter bands: ['u', 'g', 'r', 'i', 'z']

.. autodata:: frb.galaxies.defs.PanSTARRS_bands

   Pan-STARRS filter bands: ['g', 'r', 'i', 'z', 'y']

.. autodata:: frb.galaxies.defs.WISE_bands

   WISE infrared bands: ['W1', 'W2', 'W3', 'W4']

.. autodata:: frb.galaxies.defs.VISTA_bands

   VISTA near-infrared bands: ['Y', 'J', 'H', 'Ks']

Data Validation
---------------

The definitions in this module are used throughout the package for:

**Input Validation**
   - Checking that dictionary keys correspond to recognized fields
   - Ensuring consistent naming conventions across modules
   - Preventing typos in field names

**Data Integrity**  
   - Validating measurements against expected ranges
   - Cross-checking reference formats
   - Maintaining compatibility with external databases

**Documentation**
   - Providing complete lists of supported measurements
   - Defining units and conventions for each quantity
   - Enabling automatic documentation generation

Usage Examples
--------------

Validation in FRBGalaxy objects:

.. code-block:: python

    from frb.galaxies import defs
    from frb.galaxies.frbgalaxy import FRBGalaxy
    
    # Create galaxy object
    galaxy = FRBGalaxy(ra=180.0, dec=45.0, frb=frb_object)
    
    # Add photometry - keys must be in valid_photom
    galaxy.photom['DES_g'] = 22.5
    galaxy.photom['DES_r'] = 21.8  
    galaxy.photom['DES_i'] = 21.3
    
    # Validate photometry dictionary
    valid_keys = set(galaxy.photom.keys()).issubset(set(defs.valid_photom))
    print(f"All photometry keys valid: {valid_keys}")

Checking available measurements:

.. code-block:: python

    # See what emission lines are supported
    print("Supported emission lines:")
    for line in defs.valid_neb_lines[:10]:  # First 10
        print(f"  {line}")
    
    # Check derived property options  
    print("\\nPhysical properties from SED fitting:")
    for prop in defs.valid_derived_photom[:8]:
        print(f"  {prop}")

Working with survey bands:

.. code-block:: python

    # Build photometry for specific survey
    des_filters = ['DES_' + band for band in defs.DES_bands]
    print(f"DES filters: {des_filters}")
    
    # Check if galaxy has complete DES photometry
    has_des = all(filt in galaxy.photom for filt in des_filters)
    print(f"Complete DES photometry: {has_des}")

Validation during data ingestion:

.. code-block:: python

    def validate_galaxy_data(data_dict, data_type):
        """Validate galaxy data against definitions"""
        
        if data_type == 'photom':
            valid_list = defs.valid_photom + defs.valid_flux + defs.valid_ref
        elif data_type == 'derived':
            valid_list = defs.valid_derived_photom + defs.valid_derived_nebular
        elif data_type == 'neb_lines':
            valid_list = defs.valid_neb_lines + defs.valid_neb_ref
        elif data_type == 'morphology':
            valid_list = defs.valid_morphology
        else:
            return False
            
        invalid_keys = set(data_dict.keys()) - set(valid_list)
        if invalid_keys:
            print(f"Invalid keys found: {invalid_keys}")
            return False
        return True
    
    # Example usage
    photom_data = {'DES_g': 22.1, 'DES_r': 21.5, 'invalid_filter': 20.0}
    is_valid = validate_galaxy_data(photom_data, 'photom')