FRB Hosts
=========

.. automodule:: frb.galaxies.hosts
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module provides specialized functionality for FRB host galaxies, extending
the basic FRBGalaxy class with additional methods specific to confirmed host
associations and enhanced analysis capabilities.

.. note::
   This module builds upon `frb.galaxies.frbgalaxy` and provides host-specific
   analysis tools and database interfaces.

Classes
-------

FRBHost
~~~~~~~

.. autoclass:: frb.galaxies.hosts.FRBHost
   :members:
   :special-members: __init__
   :show-inheritance:

   Specialized class for confirmed FRB host galaxies, extending FRBGalaxy
   with additional functionality for:
   
   * Enhanced database integration
   * Host-specific analysis methods  
   * Association probability calculations
   * Literature compilation features

Host Analysis Functions
-----------------------

Association Analysis
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.hosts.calc_association_prob

   Calculate the probability that a galaxy is the true host of an FRB.
   
   Uses multiple factors including:
   
   * Angular offset from FRB position
   * Galaxy surface density in field
   * Magnitude-dependent number counts
   * Redshift compatibility with FRB dispersion measure

.. autofunction:: frb.galaxies.hosts.host_candidate_ranking

   Rank potential host galaxies in an FRB field by association probability.
   
   Produces ranked list considering:
   
   * Positional offsets and uncertainties
   * Galaxy properties (magnitude, morphology)
   * Field galaxy density
   * Prior expectations from FRB population studies

Database Integration
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.hosts.load_host_database

   Load the FRB host galaxy database with all confirmed associations.
   
   Provides access to:
   
   * Photometric measurements across multiple surveys
   * Spectroscopic redshifts and derived properties
   * Morphological parameters from imaging
   * Literature references and discovery papers

.. autofunction:: frb.galaxies.hosts.update_host_database

   Update host database with new measurements or revised values.

.. autofunction:: frb.galaxies.hosts.query_hosts_by_property

   Query host database by specific galaxy properties or FRB characteristics.

Population Analysis
~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.hosts.host_mass_function

   Calculate the stellar mass function of FRB host galaxies.
   
   Compares host mass distribution to field galaxy populations,
   accounting for survey selection effects and completeness.

.. autofunction:: frb.galaxies.hosts.host_sfr_distribution

   Analyze star formation rate distribution of host galaxies.

.. autofunction:: frb.galaxies.hosts.offset_distribution_analysis

   Statistical analysis of FRB-host offset distributions.
   
   Includes:
   
   * Comparison with galaxy light profiles
   * Offset vs. host properties correlations  
   * Population synthesis modeling

Literature Compilation
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: frb.galaxies.hosts.compile_literature_data

   Compile measurements from literature for specific host galaxies.
   
   Searches and consolidates:
   
   * Published photometry across papers
   * Spectroscopic measurements and redshifts
   * Morphological analyses
   * Derived physical properties

.. autofunction:: frb.galaxies.hosts.cross_match_catalogs

   Cross-match host positions with major survey catalogs.

Host-Specific Properties
------------------------

The FRBHost class includes additional attributes:

**Association Metadata**
   * Discovery paper reference
   * Association method (statistical, spectroscopic confirmation)
   * Confidence level or probability
   * Alternative host candidates

**Enhanced Measurements**
   * Compiled literature photometry
   * Multiple redshift estimates with references
   * Morphological measurements from different studies
   * Environmental context (group/cluster membership)

**Analysis Results**  
   * SED fitting results from multiple codes
   * Spectral line analysis summaries
   * Host-normalized offset measurements
   * Population comparison statistics

Examples
--------

Creating FRBHost objects:

.. code-block:: python

    from frb.galaxies.hosts import FRBHost
    from frb.frb import FRB
    
    # Create FRB object
    frb = FRB.by_name('FRB180924')
    
    # Create host object with enhanced functionality
    host = FRBHost(ra=349.24, dec=-40.9, frb=frb)
    
    # Set confirmed host status
    host.association_prob = 0.99
    host.discovery_ref = '2019Sci...365..565B'

Loading from host database:

.. code-block:: python

    from frb.galaxies.hosts import load_host_database
    
    # Load complete host database
    host_db = load_host_database()
    
    # Access specific host
    frb180924_host = host_db['FRB180924']
    
    print(f"Host redshift: {frb180924_host.z:.4f}")
    print(f"Stellar mass: {frb180924_host.derived['Mstar']:.2e} Msun")

Association probability calculation:

.. code-block:: python

    from frb.galaxies.hosts import calc_association_prob
    
    # Calculate association probability for candidate
    prob = calc_association_prob(
        offset_arcsec=1.2,
        galaxy_mag=23.1,
        field_density=1500,  # galaxies per sq arcmin to this depth
        survey_depth=25.0
    )
    
    print(f"Association probability: {prob:.3f}")

Host population analysis:

.. code-block:: python

    from frb.galaxies.hosts import host_mass_function
    import matplotlib.pyplot as plt
    
    # Calculate host stellar mass function
    masses, phi, phi_err = host_mass_function(
        completeness_limit=1e9,  # Msun
        volume_correction=True
    )
    
    # Plot comparison with field galaxies
    plt.errorbar(masses, phi, yerr=phi_err, label='FRB hosts')
    plt.xlabel('Stellar Mass [Msun]')
    plt.ylabel('Φ [Mpc^-3 dex^-1]')
    plt.yscale('log')
    plt.legend()

Literature compilation:

.. code-block:: python

    from frb.galaxies.hosts import compile_literature_data
    
    # Compile all literature data for specific host
    lit_data = compile_literature_data('FRB121102')
    
    print("Literature photometry:")
    for paper, data in lit_data['photometry'].items():
        print(f"  {paper}: {len(data)} measurements")
    
    print("\\nRedshift measurements:")  
    for z_entry in lit_data['redshifts']:
        print(f"  z = {z_entry['z']:.4f} ± {z_entry['z_err']:.4f} ({z_entry['ref']})")

Candidate ranking:

.. code-block:: python

    from frb.galaxies.hosts import host_candidate_ranking
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    # Define FRB position and error
    frb_coord = SkyCoord(ra=82.998, dec=33.148, unit='deg')
    frb_error = 0.1 * u.arcsec  # localization uncertainty
    
    # List of galaxy candidates with positions and magnitudes  
    candidates = [
        {'coord': SkyCoord(ra=82.999, dec=33.149, unit='deg'), 'mag': 22.1},
        {'coord': SkyCoord(ra=83.001, dec=33.146, unit='deg'), 'mag': 23.8},
        {'coord': SkyCoord(ra=82.995, dec=33.151, unit='deg'), 'mag': 24.2}
    ]
    
    # Rank by association probability
    ranked_candidates = host_candidate_ranking(
        frb_coord, frb_error, candidates,
        field_density=2000,
        magnitude_limit=25.0
    )
    
    print("Ranked host candidates:")
    for i, candidate in enumerate(ranked_candidates):
        print(f"  {i+1}. P = {candidate['prob']:.3f}, "
              f"offset = {candidate['offset']:.2f}\", "
              f"mag = {candidate['mag']:.1f}")

Cross-matching with surveys:

.. code-block:: python

    from frb.galaxies.hosts import cross_match_catalogs
    
    # Cross-match host position with major surveys
    matches = cross_match_catalogs(
        host.coord,
        radius=2.0 * u.arcsec,
        surveys=['DES', 'Pan-STARRS', 'WISE', 'GALEX']
    )
    
    print("Survey matches:")
    for survey, data in matches.items():
        if len(data) > 0:
            print(f"  {survey}: {len(data)} sources")
            print(f"    Closest: {data[0]['separation']:.2f}\" ")

Statistical analysis:

.. code-block:: python

    from frb.galaxies.hosts import offset_distribution_analysis
    import numpy as np
    
    # Analyze offset distribution for all confirmed hosts
    analysis_results = offset_distribution_analysis(
        normalize_by_size=True,  # Normalize by galaxy effective radius
        compare_to_light=True,   # Compare with surface brightness profiles
        bootstrap_errors=True
    )
    
    print(f"Median normalized offset: {analysis_results['median_norm']:.2f} R_e")
    print(f"Fraction within 1 R_e: {analysis_results['frac_1Re']:.2f}")
    print(f"KS test vs exponential profile: p = {analysis_results['ks_pvalue']:.3f}")