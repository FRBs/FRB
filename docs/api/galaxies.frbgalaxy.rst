FRBGalaxy Class
===============

.. automodule:: frb.galaxies.frbgalaxy
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module provides the core functionality for handling galaxies related to Fast Radio Bursts (FRBs).
It contains the main `FRBGalaxy` class which serves as a parent class for galaxies in FRB fields,
providing a simple object to hold key observable and derived quantities.

Classes
-------

FRBGalaxy
~~~~~~~~~

.. autoclass:: frb.galaxies.frbgalaxy.FRBGalaxy
   :members:
   :special-members: __init__
   :show-inheritance:

   The `FRBGalaxy` class is designed to hold key observable and derived quantities
   for galaxies associated with FRB events. 

   **Key Attributes:**

   * `redshift` (dict): Redshift measurements and estimates
   * `photom` (dict): Photometric data across multiple bands  
   * `morphology` (dict): Morphological properties
   * `neb_lines` (dict): Nebular emission line measurements
   * `kinematics` (dict): Kinematic measurements
   * `derived` (dict): Derived physical quantities
   * `offsets` (dict): Positional offsets from FRB coordinates
   * `positional_error` (dict): Astrometric and source position errors

   .. warning::
      Generating hundreds of these objects will likely be slow, especially 
      due to SkyCoord generation. A new class will be warranted for that use case.

Key Methods
-----------

Class Methods
~~~~~~~~~~~~~

.. automethod:: frb.galaxies.frbgalaxy.FRBGalaxy.from_dict

   Instantiate an FRBGalaxy object from a dictionary containing galaxy parameters.

.. automethod:: frb.galaxies.frbgalaxy.FRBGalaxy.from_json

   Load an FRBGalaxy object from a JSON file.

Instance Methods
~~~~~~~~~~~~~~~~

.. automethod:: frb.galaxies.frbgalaxy.FRBGalaxy.set_z

   Set the redshift value(s) with specified origin (spectroscopic or photometric).

.. automethod:: frb.galaxies.frbgalaxy.FRBGalaxy.calc_nebular_lum

   Calculate line luminosity with optional dust extinction correction.

.. automethod:: frb.galaxies.frbgalaxy.FRBGalaxy.run_cigale

   Run CIGALE SED fitting analysis on the galaxy's photometry.

.. automethod:: frb.galaxies.frbgalaxy.FRBGalaxy.parse_photom

   Parse photometry from an input table and populate the photom dictionary.

.. automethod:: frb.galaxies.frbgalaxy.FRBGalaxy.vet_one

   Validate one of the main attribute dictionaries against allowed values.

Properties
----------

.. autoproperty:: frb.galaxies.frbgalaxy.FRBGalaxy.z

   Return the preferred redshift of the galaxy.

.. autoproperty:: frb.galaxies.frbgalaxy.FRBGalaxy.z_err

   Return the error in the preferred redshift.

Examples
--------

Creating an FRBGalaxy instance:

.. code-block:: python

    from frb.frb import FRB
    from frb.galaxies.frbgalaxy import FRBGalaxy
    
    # Create FRB object
    frb = FRB.by_name('FRB180924')
    
    # Create galaxy object
    galaxy = FRBGalaxy(ra=349.24, dec=-40.9, frb=frb)
    
    # Set redshift
    galaxy.set_z(0.3214, 'spec', err=0.0001)

Loading from dictionary:

.. code-block:: python

    # Load from dictionary
    galaxy_dict = {
        'ra': 349.24, 
        'dec': -40.9,
        'cosmo': 'Planck18'
    }
    galaxy = FRBGalaxy.from_dict(frb, galaxy_dict)