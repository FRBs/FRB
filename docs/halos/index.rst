frb.halos Package Documentation
===============================

Modules in the ``frb.halos`` folder provide tools 
for modeling galaxy halos and calculating 
their contributions to Fast Radio Burst (FRB) dispersion measures. 
This package includes halo models, stellar-halo mass relations, 
photometric redshift analysis, and statistical tools for halo populations.

Overview
--------

The halos package enables:

* **Halo Modeling**: Modified NFW profiles for galaxy halos, including specialized models for the Milky Way, M31, and galaxy clusters
* **Mass Relations**: Stellar-halo mass relations from literature (Moster+2013, Kravtsov+2018)
* **DM Calculations**: Dispersion measure contributions from individual halos and halo populations
* **Photo-z Analysis**: Integration with photometric redshift and SED fitting tools
* **Statistical Analysis**: Halo mass functions and encounter rates along FRB sightlines

Modules
-------

Core Models
~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   models

Main module containing halo density profiles, stellar-halo mass relations, and galaxy-specific models.

**Key Features:**
* ModifiedNFW class for customizable halo profiles
* Stellar-halo mass relations with scatter
* Specialized models for Milky Way, M31, and galaxy clusters
* DM calculation methods

Photometric Analysis  
~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   photoz

Tools for photometric redshift-based halo analysis of FRB fields.

**Key Features:**
* DES photometry retrieval and processing
* EAZY and CIGALE integration for photo-z and stellar masses
* 3D interpolation grids for efficient DM calculations
* Full pipeline analysis for FRB fields

Statistical Tools
~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   hmf

Halo mass function calculations and statistical analysis.

.. warning::
   This module is deprecated. Functionality moved to frb.halos.models.

**Key Features:**
* Matter fraction in collapsed halos
* Halo encounter statistics along sightlines
* Redshift evolution of halo populations

Quick Start
-----------

Basic halo DM calculation::

    from frb.halos.models import ModifiedNFW
    from astropy import units as u
    
    # Create halo model
    halo = ModifiedNFW(log_Mhalo=12.0, z=0.3)
    
    # Calculate DM at 50 kpc impact parameter
    dm = halo.Ne_Rperp(50 * u.kpc)
    print(f"DM contribution: {dm}")

Stellar-halo mass conversion::

    from frb.halos.models import halomass_from_stellarmass
    
    log_mstar = 10.5  # log solar masses
    z = 0.5
    log_mhalo = halomass_from_stellarmass(log_mstar, z=z)
    print(f"Halo mass: {log_mhalo:.2f}")

Complete field analysis::

    from frb.frb import FRB
    from frb.halos.photoz import full_analysis
    
    # Load FRB and analyze field
    frb = FRB.by_name('FRB20180924B')
    full_analysis(frb, 'field_photometry.fits', './results/')

Common Workflows
----------------

**Individual Galaxy Analysis:**

1. Load galaxy photometry and redshift
2. Estimate stellar mass (from SED fitting)
3. Convert to halo mass using SHMR
4. Create ModifiedNFW model 
5. Calculate DM contribution

**FRB Field Analysis:**

1. Retrieve field photometry (get_des_data)
2. Run photo-z analysis (EAZY)  
3. Fit SEDs for stellar masses (CIGALE)
4. Generate halo mass realizations
5. Calculate statistical DM contributions

**Population Studies:**

1. Define mass and redshift ranges
2. Calculate halo fractions (frac_in_halos)
3. Estimate encounter rates (halo_incidence)
4. Model cumulative effects

Dependencies
------------

**Required:**
* astropy
* numpy  
* scipy

**Optional for full functionality:**
* hmf_emulator (for halo mass functions)
* pathos (for multiprocessing)
* progressbar2 (for progress tracking)
* threedhst (for EAZY integration)

**External codes:**
* EAZY (photometric redshifts)
* CIGALE (SED fitting)

Related Modules
---------------

The halos package integrates with:

* :mod:`frb.frb` - FRB object definitions
* :mod:`frb.galaxies` - Galaxy analysis tools
* :mod:`frb.surveys` - Survey data access
* :mod:`frb.dm` - DM calculations and components

References
----------

Key papers implemented in this package:

* **Moster+2013**: Stellar-halo mass relations
* **Kravtsov+2018**: Alternative SHMR at low redshift  
* **Mathews & Prochaska 2017**: Modified NFW profiles
* **Miller & Bregman 2015**: ICM models for clusters
* **Tinker+2008**: Halo mass function (via Aemulus)

Examples
--------

Advanced halo modeling::

    from frb.halos.models import ModifiedNFW, halomass_from_stellarmass
    from astropy import units as u
    import numpy as np
    
    # Galaxy parameters
    log_mstar = 10.8
    z_gal = 0.4
    
    # Convert to halo mass with scatter
    log_mhalos = []
    for i in range(100):
        log_mh = halomass_from_stellarmass(log_mstar, z=z_gal, randomize=True)
        log_mhalos.append(log_mh)
    
    # Create halo model with mean mass
    mean_mhalo = np.mean(log_mhalos)
    halo = ModifiedNFW(log_Mhalo=mean_mhalo, z=z_gal, 
                       f_hot=0.6, alpha=2, y0=2)
    
    # Calculate DM vs impact parameter
    offsets = np.logspace(0, 3, 50) * u.kpc  # 1-1000 kpc
    dms = []
    for offset in offsets:
        dm = halo.Ne_Rperp(offset)
        dms.append(dm.to('pc/cm**3').value)
    
    # Plot results
    import matplotlib.pyplot as plt
    plt.loglog(offsets.value, dms)
    plt.xlabel('Impact parameter [kpc]')
    plt.ylabel('DM contribution [pc/cm³]')
    plt.title(f'Halo DM profile (log Mhalo = {mean_mhalo:.1f})')
    plt.show()

Field-wide statistical analysis::

    from frb.halos.photoz import get_des_data, dm_grid
    from frb.halos.models import frac_in_halos
    from frb.frb import FRB
    from astropy import units as u
    import numpy as np
    
    # Load FRB
    frb = FRB.by_name('FRB20180924B')
    
    # Get field photometry
    field_cat = get_des_data(frb.coord, radius=10*u.arcmin)
    print(f"Field contains {len(field_cat)} galaxies")
    
    # Create DM interpolation grid
    dm_grid(frb.z, n_z=100, n_o=100, n_m=100)
    
    # Calculate halo statistics
    z_range = np.linspace(0.1, frb.z, 50)
    mass_ranges = [(1e11, 1e12), (1e12, 1e13), (1e13, 1e14)]
    
    for m_low, m_high in mass_ranges:
        fractions = frac_in_halos(z_range, m_low, m_high)
        print(f"Mass range {m_low:.0e}-{m_high:.0e}:")
        print(f"  Peak fraction: {np.max(fractions):.3f} at z={z_range[np.argmax(fractions)]:.2f}")

See Also
--------

* :doc:`../api/galaxies` - Galaxy analysis tools
* :doc:`../api/dm` - DM calculation modules
* :doc:`../quickstart` - Quick start guide