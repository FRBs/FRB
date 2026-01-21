**************
FRBHost Class
**************

The ``FRBHost`` class represents the host galaxy of an FRB. It inherits from
``FRBGalaxy`` and stores photometry, redshifts, morphology, nebular line
measurements, and derived physical properties.

Overview
========

The ``FRBHost`` class provides methods to:

* Load host galaxy data from the repository
* Access photometric measurements and derived quantities
* Calculate nebular properties (SFR, extinction)
* Interface with SED fitting tools (CIGALE, EAZY)
* Retrieve spectra from the specDB archive

Loading Host Galaxies
=====================

by_frb
------

The recommended way to load a host galaxy is through its associated FRB::

    from frb.frb import FRB

    # Load the FRB first
    frb = FRB.by_name('FRB20180924B')

    # Get the host galaxy
    host = frb.grab_host()

    # Print basic info
    print(host)
    # <FRBHost: 21:44:25.25 -40:54:00.80, FRB=FRB20180924B z=0.3212>

Alternatively, use the ``FRBHost.by_frb()`` class method directly::

    from frb.frb import FRB
    from frb.galaxies.frbgalaxy import FRBHost

    frb = FRB.by_name('FRB20180924B')
    host = FRBHost.by_frb(frb)

from_json
---------

Load from a specific JSON file::

    from frb.frb import FRB
    from frb.galaxies.frbgalaxy import FRBHost

    frb = FRB.by_name('FRB20180924B')
    host = FRBHost.from_json(frb, '/path/to/FRB20180924B_host.json')

Key Attributes
==============

Host galaxies have several attribute dictionaries containing measurements
and derived quantities.

Coordinates
-----------

::

    # Host galaxy coordinates
    host.coord
    # <SkyCoord (ICRS): (ra, dec) in deg (321.60521, -40.90022)>

    # Associated FRB
    host.frb.frb_name
    # 'FRB20180924B'

    # Host name
    host.name
    # 'HG20180924B'

Redshift
--------

::

    # Best redshift
    host.z
    # 0.3212

    # Redshift error
    host.z_err

    # Full redshift dict
    host.redshift
    # {'z': 0.3212, 'z_spec': 0.3212, 'z_FRB': 0.3212}

Photometry
----------

Photometric measurements are stored in the ``photom`` dict::

    # All photometry
    host.photom
    # {'SDSS_u': 21.45, 'SDSS_u_err': 0.15, 'SDSS_g': 20.23, ...}

    # Individual bands
    host.photom['SDSS_r']
    host.photom['SDSS_r_err']

    # Fluxes (in mJy) are also stored
    host.photom['SDSS_r_flux']
    host.photom['SDSS_r_flux_err']

Derived Quantities
------------------

Physical properties derived from SED fitting or spectroscopy::

    # All derived quantities
    host.derived
    # {'Mstar': 1.2e10, 'Mstar_err': 2e9, 'SFR_photom': 1.5, ...}

    # Stellar mass (solar masses)
    host.derived['Mstar']

    # Star formation rate (Msun/yr)
    host.derived['SFR_photom']    # From SED fitting
    host.derived['SFR_nebular']   # From emission lines

    # Extinction
    host.derived['AV_nebular']
    host.derived['EBV_photom']

Nebular Emission Lines
----------------------

Emission line fluxes (erg/s/cm^2)::

    # All line measurements
    host.neb_lines
    # {'Halpha': 1.2e-16, 'Halpha_err': 1e-17, '[OIII] 5007': 5e-17, ...}

    # Individual lines
    host.neb_lines['Halpha']
    host.neb_lines['Hbeta']
    host.neb_lines['[NII] 6584']

Morphology
----------

Galaxy structural parameters (typically from GALFIT)::

    host.morphology
    # {'reff_ang': 0.85, 'reff_kpc': 3.2, 'n': 1.5, 'b/a': 0.7, ...}

    # Effective radius
    host.morphology['reff_ang']  # arcsec
    host.morphology['reff_kpc']  # kpc

    # Sersic index
    host.morphology['n']

Offsets
-------

Angular and physical offsets between FRB and host::

    host.offsets
    # {'ang_avg': 0.5, 'ang_best': 0.45, 'physical': 2.1, ...}

    # Physical offset in kpc
    host.offsets['physical']
    host.offsets['physical_err']

Calculating Derived Quantities
==============================

Nebular SFR
-----------

Calculate star formation rate from emission lines::

    # First calculate extinction (optional but recommended)
    host.calc_nebular_AV(method='Ha/Hb')
    print(f"AV = {host.derived['AV_nebular']}")

    # Calculate SFR from H-alpha
    host.calc_nebular_SFR(method='Ha')
    print(f"SFR = {host.derived['SFR_nebular']} Msun/yr")

Line Luminosities
-----------------

Calculate emission line luminosities::

    Lum, Lum_err = host.calc_nebular_lum('Halpha')
    print(f"L(Ha) = {Lum}")

Halo DM
-------

Calculate the halo contribution to DM::

    DM_halo = host.calc_dm_halo()
    print(f"DM_halo = {DM_halo}")

Retrieving Spectra
==================

If spectra are available in the specDB archive::

    # Get spectrum and metadata
    meta, spec = host.get_metaspec()

    # meta is an astropy Table with spectrum info
    print(meta)

    # spec is an XSpectrum1D object
    spec.wavelength  # Wavelength array
    spec.flux        # Flux array

    # Get all spectra (if multiple exist)
    meta, spec = host.get_metaspec(return_all=True)

    # Specify instrument
    meta, spec = host.get_metaspec(instr='MUSE')

Running SED Fitting
===================

CIGALE
------

Run CIGALE SED fitting directly::

    host.run_cigale(
        data_file='cigale_input.fits',
        config_file='pcigale.ini',
        wait_for_input=False,
        save_sed=True,
        plot=True,
        outdir='cigale_output/'
    )

Parse CIGALE results::

    host.parse_cigale('cigale_output/results.txt')
    print(host.derived['Mstar'])
    print(host.derived['SFR_photom'])

Parsing External Results
========================

GALFIT
------

Parse GALFIT output for morphology::

    host.parse_galfit('galfit_output.fits')
    print(host.morphology['reff_kpc'])
    print(host.morphology['n'])

pPXF
----

Parse pPXF spectral fitting results::

    host.parse_ppxf('ppxf_results.ecsv')
    print(host.neb_lines['Halpha'])

Building a Table of Hosts
=========================

Create a pandas DataFrame with all host galaxy properties::

    from frb.galaxies import utils

    # Build the table
    host_tbl, tbl_units = utils.build_table_of_hosts()

    # View columns
    print(host_tbl.columns.tolist())
    # ['Host', 'FRBname', 'RA_host', 'DEC_host', 'FRBobj', 'Mstar', 'SFR_photom', ...]

    # Filter by stellar mass
    massive_hosts = host_tbl[host_tbl['Mstar'] > 1e10]

    # Get hosts with spectroscopic redshifts
    spec_z = host_tbl[host_tbl['z_spec'].notna()]

    # Merge with FRB table
    from frb.frb import build_table_of_frbs
    frb_tbl, _ = build_table_of_frbs()

    import pandas as pd
    merged = pd.merge(frb_tbl, host_tbl, left_on='FRB', right_on='FRBname')

Default columns include data from:

* **derived**: Mstar, SFR_photom, SFR_nebular, AV_nebular, M_r, etc.
* **photom**: All photometric bands and fluxes
* **neb_lines**: Emission line fluxes
* **offsets**: Angular and physical offsets
* **morphology**: Structural parameters
* **redshift**: z, z_spec, z_phot, z_FRB

Listing All Hosts
=================

::

    from frb.galaxies.utils import list_of_hosts

    # Get all hosts
    frbs, hosts = list_of_hosts()

    print(f"Number of hosts: {len(hosts)}")

    # Iterate
    for frb, host in zip(frbs, hosts):
        print(f"{frb.frb_name}: z={host.z}, Mstar={host.derived.get('Mstar', 'N/A')}")

Loading PATH Results
====================

Load probabilistic association results::

    from frb.galaxies import utils

    path_tbl = utils.load_PATH()
    print(path_tbl.columns)
    # ['FRB', 'RA', 'Dec', 'P_Ox', 'P_O', 'ang_size', ...]

Writing to JSON
===============

Save a host galaxy object::

    host.write_to_json(path='./output/')
    # Writes to FRB20180924B_host.json

Setting Redshifts
=================

::

    # Set spectroscopic redshift
    host.set_z(0.3212, 'spec', err=0.0001)

    # Set photometric redshift
    host.set_z(0.35, 'phot', err=0.05)

API Reference
=============

FRBGalaxy (Parent Class)
------------------------

.. autoclass:: frb.galaxies.frbgalaxy.FRBGalaxy
   :members:
   :undoc-members:
   :show-inheritance:

FRBHost Class
-------------

.. autoclass:: frb.galaxies.frbgalaxy.FRBHost
   :members:
   :undoc-members:
   :show-inheritance:

Utility Functions
-----------------

.. autofunction:: frb.galaxies.utils.build_table_of_hosts

.. autofunction:: frb.galaxies.utils.list_of_hosts

.. autofunction:: frb.galaxies.utils.load_PATH

See Also
========

* :doc:`frb_class` - FRB class documentation
* :doc:`galaxies` - Galaxy analysis overview
* :doc:`database` - Database utilities
