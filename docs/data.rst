****
Data
****

Overview
========

One purpose of this repository is to provide data
related to FRBs, their host galaxies, and the galaxies
foreground to them.  This includes measurements
(e.g. photometry), derived quantities (e.g. star formation
rate), and observational data (e.g. spectra).

FRBs
====

For FRBs included in this repository, currently those that
are well localized and published, we have archived a set
of basic measaurements, e.g. coordinates.  These are
saved as a set of JSON files in the data/FRBs/ folder
of the repository.

One can load these data into an FRB object as follows::

    frb121102 = ffrb.FRB.by_name('FRB121102')
    # Coordinate
    frb121102.coord
    # Error ellipse
    frb121101.eelipse
    # DM
    frb121101.DM

RM and other measurements are also included
when available.  See the FRB_Event.ipynb notebook
for a bit more.

Host Galaxies
=============

When a high probability association has been made for an
FRB to its host galaxy, we also include data on the latter.

Here is an example for 180924::

    # Instantiate the FRB
    frb180924 = frb.FRB.by_name('FRB180924')
    # Grab its host
    hg180924 = frb180924.grab_host()
    # Derived properties
    hg180924.derived
    # Grab the spectrum (see data access below)
    meta, spec = hg180924.get_metaspec()

Future code will generate tables of the key quantities
for the galaxies.

Spectra
=======

SpecDB
------

As galaxy spectra related to FRB surveys becomes available,
we intend to archive these within a
`specdb <https://specdb.readthedocs.io/en/latest/>`_
database file.

Here is the
`public specdb <https://drive.google.com/file/d/14Wx4ctpxHRDEI9joVzHGidtiO3spg5fb/view?usp=sharing>`_
which currently includes galaxy spectra related
to all of the published ASKAP/CRAFT FRBs.

You will need to:

#. Install `specdb <https://specdb.readthedocs.io/en/latest/>`_
#. Place the `public specdb <https://drive.google.com/file/d/14Wx4ctpxHRDEI9joVzHGidtiO3spg5fb/view?usp=sharing>`_ file in a folder
#. Point the environmental variable SPECDB to that folder

Galaxy Spectrum
---------------

The easiest way (perhaps) to load up a spectrum is
by first instantiating an FRB galaxy object.  Here
is an example for the host galaxy of FRB 180924::

    # Load the FRB
    frb180924 = frb.FRB.by_name('FRB180924')
    # Load the host galaxy
    hg180924 = frb180924.grab_host()
    # Load a meta data Table and the spectra
    meta, spec = hg180924.get_metaspec()

*meta* is an astropy Table describing all of the archived spectra
for this galaxy (here only 1 spectrum).  *spec* is an
XSpectrum1D object from `linetools <https://github.com/linetools/linetools>`_.

Galaxy script
-------------

The FRB repo also provides a basic script -- frb_galaxies -- for accessing galaxy spectra
in the *specdb* archive.  Here is the usage::

    usage: frb_galaxies [-h] [--rho RHO] [--ang_offset ANG_OFFSET] [--cat]
                    [--specdb SPECDB] [-p]
                    coord

    Script to fuss with FRB galaxies [v1.1]

    positional arguments:
      coord                 Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322
                            or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)

    optional arguments:
      -h, --help            show this help message and exit
      --rho RHO             Maximum impact parameter in kpc [default=300.]
      --ang_offset ANG_OFFSET
                            Maximum offset in arcsec [over-rides --rho if set]
      --cat                 Only show data from the catalog (not meta)
      --specdb SPECDB       specDB file; defaults to $SPECDB/FRB_specdb.hdf5
      -p, --plot            Launch a plotting GUI?

And here is an example call::

    frb_galaxies FRB180924

This prints a brief summary of the spectra available
in the field surrounding FRB180924 (default is a 300kpc
radius).  You can plot spectra by adding the -p option.