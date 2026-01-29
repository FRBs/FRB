.. highlight:: rest

*******
Surveys
*******

This document describes software used to
parse catalogs and data from public surveys.
In the context of FRB research, of course.

These codes have extra dependencies that
are detailed in the Installation notes.

Usage
=====

Catalog
-------

Here is an example of grabbing a catalog of sources around
an input coordinate from the :ref:`surveys-des` survey::

    from astropy.coordinates import SkyCoord
    from astropy import units
    #
    coord = SkyCoord('J214425.25-403400.81', unit=(units.hourangle, units.deg))
    search_r = 10 * units.arcsec
    #
    des_srvy = survey_utils.load_survey_by_name('DES', coord, search_r)
    des_tbl = des_srvy.get_catalog(print_query=True)


Cutout
------

Here is how to retrieve a cut-out image::

    cutout, cutout_hdr = des_srvy.get_cutout(search_r, band='r')

The header may be used to generate a WCS object and
overlay coordinates, etc.


Available Surveys
=================

.. _surveys-des:

DES
---

The Dark Energy Survey (DES) is an imaging survey
with the Dark Energy Camera on the CTIO telescope
in Chile.  The current public dataset is
Data Release 1 (DR1).  The software will slurp
catalogs and images from DES-DR1.

DECaL
-----

The DECaL survey imaged the majority of the Northern,
extragalactic sky in the g,r, and z bands for the
DESI survey.

HEASARC and SkyView
-------

Any of the catalogs maintained by HEASARC and
images maintained by SkyView are
in principle available.  These are the ones
that have been integrated thus far:

NVSS
++++

The  `NVSS <https://www.cv.nrao.edu/nvss/>`_ survey
imaged the entire Northern sky north of -40 deg with the VLA.
Catalogs and cutouts are available.

FIRST
+++++

The `FIRST <http://sundog.stsci.edu/>`_
survey imaged the Northern sky at 1.4 GHz with the VLA.
Catalogs and cutouts are available.

WENSS
+++++

The `WENSS <https://heasarc.gsfc.nasa.gov/w3browse/all/wenss.html/>`_
survey imaged the sky north of +30 deg at 325 MHz with the WSRT.
Catalogs and cutouts are available.

PSRCAT
------

One can also access the pulsar catalog known as
PSRCAT which is regularly slurped into the FRB
`GitHub organization <https://github.com/FRBs/pulsars>`_.
