*********
FRB Class
*********

The ``FRB`` class is the core object representing an observed Fast Radio Burst event.
It stores coordinates, dispersion measure, rotation measure, redshift, and other
key observables for well-localized FRBs in the repository.

Overview
========

The ``FRB`` class inherits from ``GenericFRB`` and provides methods to:

* Load FRB data from the repository by name
* Access measured properties (DM, RM, coordinates, error ellipse, etc.)
* Retrieve associated host galaxy information
* Export data to JSON format

Loading FRBs
============

by_name
-------

The most common way to load an FRB is by its name using the ``by_name()`` class method::

    from frb.frb import FRB

    # Load FRB by name
    frb121102 = FRB.by_name('FRB20121102A')

    # Or without the 'FRB' prefix for older naming
    frb180924 = FRB.by_name('FRB20180924B')

    # Print basic info
    print(frb180924)
    # <FRB: FRB20180924B J214425.26-405400.1 DM=362.16 pc / cm3 z=0.3212>

from_json
---------

You can also load from a specific JSON file::

    from frb.frb import FRB

    frb = FRB.from_json('/path/to/FRB20180924B.json')

Key Attributes
==============

Once loaded, an FRB object provides access to:

Coordinates
-----------

::

    # SkyCoord object
    frb.coord
    # <SkyCoord (ICRS): (ra, dec) in deg (321.60525, -40.90003)>

    # Individual components
    frb.coord.ra   # Right Ascension
    frb.coord.dec  # Declination

    # Different formats
    frb.coord.to_string('hmsdms')
    # '21h26m25.26s -40d54m00.11s'

Dispersion Measure
------------------

::

    # Observed DM (with units)
    frb.DM
    # <Quantity 362.16 pc / cm3>

    # DM error (if available)
    frb.DM_err

    # ISM contribution from NE2001
    frb.DMISM

Rotation Measure
----------------

::

    # RM (if measured)
    frb.RM
    # <Quantity 14.5 rad / m2>

    frb.RM_err

Redshift
--------

::

    # Redshift (if known)
    frb.z
    # 0.3212

Error Ellipse
-------------

The localization uncertainty is stored in an error ellipse::

    # Error ellipse dict
    frb.eellipse
    # {'a': 0.07, 'b': 0.06, 'theta': 0.0, 'cl': 68.0, ...}

    # Combined semi-major axis (statistical + systematic)
    frb.sig_a  # arcsec

    # Combined semi-minor axis
    frb.sig_b  # arcsec

Pulse Properties
----------------

::

    # Pulse properties dict
    frb.pulse
    # Contains: Wi (intrinsic width), tscatt (scattering time), etc.

Other Properties
----------------

::

    # FRB name
    frb.frb_name
    # 'FRB20180924B'

    # Is it a repeater?
    frb.repeater
    # False

    # References
    frb.refs
    # ['Bannister2019']

Grabbing the Host Galaxy
========================

For FRBs with identified host galaxies, use the ``grab_host()`` method::

    from frb.frb import FRB

    # Load FRB
    frb = FRB.by_name('FRB20180924B')

    # Get the host galaxy object
    host = frb.grab_host()

    # Access host properties
    print(host.z)           # Redshift
    print(host.derived)     # Derived quantities (Mstar, SFR, etc.)
    print(host.photom)      # Photometry

See :doc:`frbhost_class` for details on the FRBHost class.

Listing All FRBs
================

To get a list of all FRBs in the repository::

    from frb.frb import list_of_frbs

    # Get all FRBs
    all_frbs = list_of_frbs()
    print(f"Total FRBs: {len(all_frbs)}")

    # Only those with redshifts
    frbs_with_z = list_of_frbs(require_z=True)
    print(f"FRBs with redshift: {len(frbs_with_z)}")

Building a Table of FRBs
========================

The ``build_table_of_frbs()`` function creates a pandas DataFrame
with all FRB properties::

    from frb.frb import build_table_of_frbs

    # Build the table
    frb_tbl, tbl_units = build_table_of_frbs()

    # View basic info
    print(frb_tbl.columns.tolist())
    # ['FRB', 'RA', 'DEC', 'ee_a', 'ee_b', ..., 'DM', 'z', 'RM', ...]

    # Check units
    print(tbl_units['DM'])
    # 'pc / cm3'

    # Filter by DM
    high_dm = frb_tbl[frb_tbl['DM'] > 500]
    print(f"High DM FRBs: {len(high_dm)}")

    # Get repeaters
    repeaters = frb_tbl[frb_tbl['repeater'] == True]

Default columns include:

* **FRB**: FRB name
* **RA, DEC**: Coordinates (deg)
* **ee_a, ee_b, ee_theta**: Error ellipse parameters
* **DM, DM_err**: Dispersion measure
* **z**: Redshift
* **RM, RM_err**: Rotation measure
* **DMISM**: Galactic ISM DM contribution
* **fluence**: Burst fluence
* **repeater**: Boolean flag
* **pulse_Wi, pulse_tscatt**: Pulse properties
* **refs**: References

Creating Custom FRBs
====================

You can also create FRB objects for your own data::

    from frb.frb import FRB
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    # Create coordinate
    coord = SkyCoord(ra=123.456, dec=-45.678, unit='deg')

    # Create FRB object
    my_frb = FRB('FRB20230101A', coord, DM=500*u.pc/u.cm**3)

    # Set additional properties
    my_frb.z = 0.5
    my_frb.set_ee(a=0.1, b=0.08, theta=45, cl=68)

Writing to JSON
===============

Save an FRB object to a JSON file::

    frb.write_to_json(outfile='my_frb.json', path='./')

API Reference
=============

.. autoclass:: frb.frb.FRB
   :members:
   :undoc-members:
   :show-inheritance:

.. autofunction:: frb.frb.list_of_frbs

.. autofunction:: frb.frb.build_table_of_frbs

See Also
========

* :doc:`frbhost_class` - Host galaxy class documentation
* :doc:`database` - Database access utilities
* :doc:`dm` - DM calculations
