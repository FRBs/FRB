.. _surveys:

frb.surveys — Survey Data Access
=================================

The `frb.surveys` package provides a uniform interface for querying dozens of
astronomical surveys and catalogs. It handles the backend complexity of
accessing data from various archives (NOIRLab Data Lab, HEASARC, IRSA, MAST,
SkyView) and returns standardized `astropy` objects.

For a hands-on guide with detailed examples for every survey class, see the
companion tutorial notebook:

.. toctree::
   :maxdepth: 1

   ../nb/Surveys







General Usage
-------------

All survey classes are instantiated with a sky coordinate and a search radius.
The primary methods are `get_catalog()` to retrieve a source catalog as an
`astropy.table.Table` and `get_image()` to download a FITS image as an
`astropy.io.fits.HDUList`.

Here is a typical example using Pan-STARRS:

.. code-block:: python

   from astropy.coordinates import SkyCoord
   from astropy import units as u
   from frb.surveys.panstarrs import Pan_STARRS_Survey

   # Define coordinate and search radius
   coord = SkyCoord('J081240.68+320809.0', unit=(u.hourangle, u.deg))
   radius = 10 * u.arcsec

   # Instantiate the survey object
   ps1 = Pan_STARRS_Survey(coord, radius)

   # Get the source catalog
   catalog = ps1.get_catalog()

   # Get a FITS image cutout
   image_hdu = ps1.get_image()


Utility Modules
---------------

The package includes two key utility modules for survey-agnostic operations
and catalog manipulation.

.. toctree::
   :maxdepth: 1

   surveys.survey_utils
   surveys.catalog_utils


Survey-Specific Submodules
--------------------------

The following modules contain the individual survey classes.

.. toctree::
   :maxdepth: 1

   surveys.sdss
   surveys.des
   surveys.decals
   surveys.delve
   surveys.nsc
   surveys.hsc
   surveys.panstarrs
   surveys.wise
   surveys.vista
   surveys.twomass
   surveys.galex
   surveys.euclid
   surveys.first
   surveys.nvss
   surveys.wenss
   surveys.desi
   surveys.nedlvs
   surveys.psrcat
   surveys.heasarc
   surveys.dlsurvey
   surveys.images
   surveys.cluster_search
   surveys.surveycoord
   surveys.survey_io
   surveys.tns_util
   surveys.utils_crossmatching
   surveys.defs

.. note::
   ``frb.surveys.dlsurvey`` requires the optional dependency
   ``datalab-client`` for NOIRLab Data Lab access.


.. note::
   ``frb.surveys.psrcat`` requires optional ``FRB-pulsars`` support.
