Halo Mass Function
==================

.. automodule:: frb.halos.hmf

Module for halo mass function calculations and related statistical analysis.

.. warning::
   This module is deprecated. The functionality has been moved to :mod:`frb.halos.models`.
   Please use that module instead for new code.

Functions
---------

frac_in_halos
~~~~~~~~~~~~~

.. autofunction:: frac_in_halos

   Calculate the fraction of matter in collapsed halos over a mass range and at a given redshift.

   **Parameters:**
   
   * ``zvals`` : ndarray - Redshift values
   * ``Mlow`` : float - Minimum halo mass in h^-1 units  
   * ``Mhigh`` : float - Maximum halo mass in h^-1 units
   * ``rmax`` : float, optional - Extent of halo in units of rvir (default: 1.0)

   **Returns:**
   
   * ``ratios`` : ndarray - rho_halo / rho_m

   .. note::
      The fraction of DM associated with these halos will be scaled down by an additional factor of f_diffuse.

   .. warning::
      This calculation assumes a single concentration for all halos.

   **Requirements:**
   
   * Requires Aemulus HMF to be installed

halo_incidence
~~~~~~~~~~~~~~

.. autofunction:: halo_incidence

   Calculate the (approximate) average number of intersections to halos of a given minimum mass to a given FRB redshift.

   **Parameters:**
   
   * ``Mlow`` : float - Mass of minimum halo in Solar masses (minimum: 2e10)
   * ``zFRB`` : float - Redshift of the FRB
   * ``radius`` : Quantity, optional - Physical separation from sightline
   * ``hmfe`` : object, optional - HMF emulator instance
   * ``Mhigh`` : float, optional - Maximum halo mass (default: 1e16)
   * ``nsample`` : int, optional - Number of samples (default: 20)
   * ``cumul`` : bool, optional - Return cumulative values (default: False)

   **Returns:**
   
   * Average number of halo intersections

   .. note::
      The code handles h^-1 factors automatically. If radius is not specified, 
      it uses rvir derived from Mlow.

   **Requirements:**
   
   * Requires Aemulus HMF to be installed

Examples
--------

Calculate matter fraction in halos::

    from frb.halos.hmf import frac_in_halos
    import numpy as np
    
    # Define redshift range
    redshifts = np.linspace(0.1, 1.0, 10)
    
    # Mass range for galaxy-scale halos
    mass_min = 1e11  # Solar masses  
    mass_max = 1e13
    
    # Calculate fractions
    fractions = frac_in_halos(redshifts, mass_min, mass_max)
    
    for z, frac in zip(redshifts, fractions):
        print(f"z={z:.1f}: {frac:.3f} of matter in galaxy halos")

Calculate halo incidence along sightline::

    from frb.halos.hmf import halo_incidence
    from astropy import units as u
    
    # Parameters for FRB sightline
    min_mass = 1e12  # Solar masses
    frb_redshift = 0.8
    impact_radius = 100 * u.kpc
    
    # Calculate average number of halo intersections
    n_halos = halo_incidence(min_mass, frb_redshift, 
                           radius=impact_radius)
    
    print(f"Expected {n_halos:.2f} halo intersections")

Compare different mass ranges::

    from frb.halos.hmf import frac_in_halos
    import numpy as np
    import matplotlib.pyplot as plt
    
    z_array = np.linspace(0, 2, 50)
    
    # Different mass ranges
    dwarf_range = frac_in_halos(z_array, 1e9, 1e11)   # Dwarf galaxies
    galaxy_range = frac_in_halos(z_array, 1e11, 1e13) # Normal galaxies  
    cluster_range = frac_in_halos(z_array, 1e13, 1e15) # Clusters
    
    plt.figure(figsize=(10, 6))
    plt.plot(z_array, dwarf_range, label='Dwarf halos (1e9-1e11)')
    plt.plot(z_array, galaxy_range, label='Galaxy halos (1e11-1e13)')  
    plt.plot(z_array, cluster_range, label='Cluster halos (1e13-1e15)')
    
    plt.xlabel('Redshift')
    plt.ylabel('Matter fraction in halos')
    plt.legend()
    plt.title('Halo matter fraction vs redshift')
    plt.show()

Extended halo calculations::

    from frb.halos.hmf import frac_in_halos
    import numpy as np
    
    # Calculate with extended halo profiles
    z_test = 0.5
    mass_min = 1e12
    mass_max = 1e14
    
    # Standard virial radius
    frac_1rvir = frac_in_halos([z_test], mass_min, mass_max, rmax=1.0)[0]
    
    # Extended to 2 virial radii  
    frac_2rvir = frac_in_halos([z_test], mass_min, mass_max, rmax=2.0)[0]
    
    # Extended to 3 virial radii
    frac_3rvir = frac_in_halos([z_test], mass_min, mass_max, rmax=3.0)[0]
    
    print(f"Matter fraction at z={z_test}:")
    print(f"  1 rvir: {frac_1rvir:.4f}")
    print(f"  2 rvir: {frac_2rvir:.4f}")  
    print(f"  3 rvir: {frac_3rvir:.4f}")
    
    boost_factor = frac_2rvir / frac_1rvir
    print(f"Boost factor (2 rvir): {boost_factor:.2f}")

Statistical analysis of halo encounters::

    from frb.halos.hmf import halo_incidence
    from astropy import units as u
    import numpy as np
    
    # Range of FRB redshifts
    frb_redshifts = np.linspace(0.1, 2.0, 20)
    
    # Different halo mass thresholds
    mass_thresholds = [1e11, 1e12, 1e13]  # Solar masses
    radius = 50 * u.kpc
    
    results = {}
    for mass_min in mass_thresholds:
        encounters = []
        for z_frb in frb_redshifts:
            n_enc = halo_incidence(mass_min, z_frb, radius=radius)
            encounters.append(n_enc)
        results[mass_min] = encounters
    
    # Plot results
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    
    for mass_min, encounters in results.items():
        plt.plot(frb_redshifts, encounters, 
                label=f'M > {mass_min:.0e} Msun')
    
    plt.xlabel('FRB Redshift')
    plt.ylabel('Average Number of Halo Encounters')
    plt.title(f'Halo Encounters vs FRB Redshift (R < {radius})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

Cumulative analysis::

    from frb.halos.hmf import halo_incidence
    from astropy import units as u
    
    # Parameters
    mass_min = 1e12
    z_frb = 1.0
    radius = 100 * u.kpc
    
    # Get cumulative encounters
    cumul_encounters = halo_incidence(mass_min, z_frb, radius=radius,
                                    cumul=True, nsample=100)
    
    print(f"Cumulative halo encounters: {cumul_encounters}")

Redshift evolution study::

    from frb.halos.hmf import frac_in_halos
    import numpy as np
    
    # Study evolution from z=0 to z=3
    z_range = np.linspace(0, 3, 100)
    
    # Different mass ranges
    mass_ranges = [
        (1e10, 1e11, 'Low mass'),
        (1e11, 1e12, 'Intermediate mass'), 
        (1e12, 1e13, 'High mass'),
        (1e13, 1e15, 'Cluster mass')
    ]
    
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.ravel()
    
    for i, (m_low, m_high, label) in enumerate(mass_ranges):
        fractions = frac_in_halos(z_range, m_low, m_high)
        
        axes[i].plot(z_range, fractions, 'b-', linewidth=2)
        axes[i].set_title(f'{label} halos')
        axes[i].set_xlabel('Redshift')  
        axes[i].set_ylabel('Matter fraction')
        axes[i].grid(True, alpha=0.3)
        axes[i].text(0.1, 0.9, f'{m_low:.0e} - {m_high:.0e} Msun', 
                    transform=axes[i].transAxes, 
                    bbox=dict(boxstyle='round', facecolor='wheat'))
    
    plt.tight_layout()
    plt.suptitle('Halo Matter Fraction Evolution', y=1.02)
    plt.show()