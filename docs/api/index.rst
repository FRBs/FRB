API Reference
=============

This section contains detailed documentation for all modules, classes, and functions in the FRB package.

Core Modules
------------

The FRB package is organized into several key modules:

.. toctree::
   :maxdepth: 2

   frb
   dm
   galaxies
   halos
   surveys
   associate
   frb_surveys
   turb_scattering
   rm
   mw
   em
   dlas
   frbcat
   experiment

Quick Reference
---------------

Most Common Functions
~~~~~~~~~~~~~~~~~~~~~

**Loading FRBs:**

.. code-block:: python

   from frb.frb import FRB
   frb_obj = FRB.by_name('FRB20121102A')

**DM Calculations:**

.. code-block:: python

   from frb.dm import igm
   DM_cosmic = igm.DM_cosmic(z=0.5)
   z_est = igm.z_from_DM(DM_obs)

**Building FRB Tables:**

.. code-block:: python

   from frb.frb import build_table_of_frbs
   frb_tbl, tbl_units = build_table_of_frbs()

**Scattering Analysis:**

.. code-block:: python

   from frb import turb_scattering as ts
   theta = ts.theta_mist(n_e, nu_obs)

Function Index
--------------

Core Functions by Category
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Dispersion Measure:**

.. autosummary::
   :nosignatures:

   frb.dm.igm.DM_cosmic
   frb.dm.igm.ne_cosmic  
   frb.dm.igm.z_from_DM
   frb.dm.igm.DM_halos
   frb.dm.cosmic.DM_cosmic
   frb.dm.host.DM_host

**Scattering:**

.. autosummary::
   :nosignatures:

   frb.turb_scattering.theta_mist
   frb.turb_scattering.tau_mist
   frb.turb_scattering.ne_from_tau_mist

**FRB Tables:**

.. autosummary::
   :nosignatures:

   frb.frb.build_table_of_frbs
   frb.frb.list_of_frbs

**FRB Objects:**

.. autosummary::
   :nosignatures:

   frb.FRB.by_name
   frb.FRB.grab_host

Class Hierarchy
---------------

.. code-block:: none

   FRB
   ├── FRB.by_name()
   ├── FRB.grab_host()
   └── FRB.set_ee()

   FRBHost
   ├── FRBHost.by_frb()
   ├── FRBHost.derived
   ├── FRBHost.photom
   └── FRBHost.get_metaspec()

Constants and Defaults
----------------------

The package uses several default values and constants:

**Cosmological Parameters:**

.. code-block:: python

   # Default cosmology (defined in frb.defs)
   frb_cosmo = astropy.cosmology.Planck18

**Default Scattering Parameters:**

.. code-block:: python

   # Default structure sizes
   L_default = 50 * u.kpc      # Structure size
   R_default = 1 * u.pc        # Cloud size  
   fV_default = 1.0           # Filling factor

**Physical Constants:**

The package uses astropy constants throughout for physical calculations.

Error Handling
--------------

The FRB package uses standard Python exception handling:

**Common Exceptions:**

- ``FileNotFoundError``: When FRB data files are not found
- ``ValueError``: When invalid parameters are passed to functions
- ``AttributeError``: When accessing properties not available for a particular FRB
- ``ImportError``: When optional dependencies are not available

**Error Handling Example:**

.. code-block:: python

   try:
       frb_obj = ffrb.FRB.by_name('NonExistentFRB')
   except FileNotFoundError:
       print("FRB data not found")
   except ValueError as e:
       print(f"Invalid FRB name: {e}")

Data Types
----------

The package primarily works with astropy quantities and coordinates:

**Common Data Types:**

- ``astropy.coordinates.SkyCoord`` - Sky positions
- ``astropy.units.Quantity`` - Physical quantities with units
- ``astropy.table.Table`` - Tabular data (catalogues)
- ``numpy.ndarray`` - Numerical arrays
- ``dict`` - Configuration and metadata

**Unit Conventions:**

- Dispersion Measure: ``pc / cm³``
- Distances: ``kpc``, ``Mpc``, ``Gpc``  
- Frequencies: ``MHz``, ``GHz``
- Time: ``s``, ``ms``, ``μs``
- Angles: ``deg``, ``arcmin``, ``arcsec``, ``mas``, ``μas``

Contributing to the API
-----------------------

If you're contributing new functions or classes:

1. Follow the existing documentation style
2. Include comprehensive docstrings
3. Add type hints where appropriate
4. Include examples in docstrings
5. Update this API documentation

.. seealso::

   - :doc:`../quickstart` - Get started with basic usage
   - :doc:`../dm` - Dispersion measure calculations
   - :doc:`../halos/index` - Halo modeling tools