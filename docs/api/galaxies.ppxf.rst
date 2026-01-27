pPXF
====

.. automodule:: frb.galaxies.ppxf
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

This module provides functionality for running pPXF (Penalized Pixel-Fitting) 
analyses on galaxy spectra. pPXF is used to extract stellar kinematics and 
stellar population information from absorption-line spectra.

Functions
---------

.. autofunction:: frb.galaxies.ppxf.run

   Main wrapper function for running and handling pPXF outputs.
   
   This function processes galaxy spectra through the pPXF pipeline, handling
   input/output and providing a simplified interface to the underlying pPXF code.

Parameters
----------

The `run` function accepts the following key parameters:

* **spec_file** (str or XSpectrum1D): Input spectrum file or object
* **R** (float): Spectral resolution
* **zgal** (float): Galaxy redshift
* **results_file** (str, optional): Output results filename
* **spec_fit** (str, optional): Fitted spectrum output filename  
* **chk** (bool, optional): Enable checking/validation
* **flux_scale** (float, optional): Flux scaling factor
* **atmos** (list, optional): Atmospheric absorption regions to mask
* **gaps** (list, optional): Detector gaps or bad regions to ignore
* **wvmnx** (tuple, optional): Wavelength range limits

Masking Options
~~~~~~~~~~~~~~~

The module provides flexible masking capabilities:

**Atmospheric Lines**
   Regions affected by atmospheric absorption can be masked during analysis:
   
   .. code-block:: python
   
       atmos = [[7150., 7300.], [7594., 7621.]]  # O2 bands

**Detector Gaps**  
   Bad regions or detector gaps can be excluded:
   
   .. code-block:: python
   
       gaps = [[6675., 6725.]]  # CCD gap

Usage Examples
--------------

Basic pPXF analysis:

.. code-block:: python

    from frb.galaxies import ppxf
    
    # Run pPXF on a spectrum file
    ppxf.run('galaxy_spectrum.fits', 
             R=3000.,  # Resolution
             zgal=0.1,  # Redshift
             results_file='ppxf_results.fits')

With masking:

.. code-block:: python

    # Define regions to mask
    atmos_lines = [[7594., 7621.], [6864., 6884.]]  # Telluric features
    detector_gaps = [[6675., 6725.]]  # Bad detector region
    
    ppxf.run('spectrum.fits',
             R=2500.,
             zgal=0.25,
             atmos=atmos_lines,
             gaps=detector_gaps,
             wvmnx=(4000., 9000.))  # Wavelength limits

Output Processing:

.. code-block:: python

    # Run with custom output files
    ppxf.run('galaxy.fits',
             R=3500.,
             zgal=0.15,
             results_file='kinematic_results.fits',
             spec_fit='best_fit_spectrum.fits',
             flux_scale=1e-17,  # Scale factor for flux units
             chk=True)  # Enable validation checks

Dependencies
------------

This module requires:

* ppxf package (Cappellari 2017, 2023)
* linetools for spectrum handling
* astropy for units and constants
* matplotlib for plotting capabilities
* numpy for numerical computations

.. note::
   The module uses MILES stellar library templates by default for 
   stellar population fitting.