************
FRB Galaxies
************

A significant portion of FRB science will flow
from studying galaxies related to these events.

This will include both the galaxy that hosted
the event and galaxies foreground to the event
which may imprint signatues in the signal itself.

We have thus far generated a class FRBGalaxy to
hold, manipulate, and derive key observed and
physical quantities.

Once a large enough sample exists, there will be other
objects to generate a easy-to-access database
(e.g. Table) of galaxy properties.

FRBGalaxy
=========

This is the parent class and is not intended to be instantiated
on its own, at least not often.  Instead, instantiate
one of its children.  Nevertheless, it contains the majority
of methods and attributes used by the children, described here.

Attributes
----------

The FRBGalaxy class contains a set of attributes, listed in
self.main_attr, that are a series of *dict* and are intended
to hold the primary quantities of the object.  To set a
common structure, the allowed keys within each *dict* are defined
in the frb.galaxies.defs module.  For example, the *morphology*
attribute may contain the entries listed in defs.valid_morphology
which is reproduced here::

    valid_morphology = [
        'reff_ang',   # Effective radius in arcsec; Galfit
        'reff_kpc',   # Effective radius in kpc; Galfit
        'n',          # Sersic index; Galfit
        'PA',         # Position angle (deg); Galfit
        'b/a',        # Ellipticity; Galfit
    ]

The uncertainty in each of these values is also permitted, e.g.
reff_ang_err and n_err.

Methods
-------

There are a number of methods used to modify the class and
also to calculate quantities of interest from existing
photometry and line flux measurements.

set_z
+++++

This sets the redshift of the galaxy and/or the related FRB,
e.g.::

    frbgalaxy.set_z(0.192, 'spec')

Here *spec* indicates the measurement is spectroscopic.

parse_galfit
++++++++++++

There are a series of methods that will parse output files
from standard galaxy analysis codes, e.g. CIGALE, GALFIT, pPXF.
Here is the sample call for GALFIT::

    frbgalaxy.parse_galfit('frb121102_galfit.log', 0.15)

The first item is the filename and second is the plate-scale
of the image analyzed (arcsec per pixel).  The values measured
are ingested into the frbgalaxy.morphology *dict*.

calc_nebular_SFR
++++++++++++++++

There are a few methods for deriving quantities of scientific
interest for the galaxy.  The example provided here is for the
star-formation rate (SFR) using a nebular line flux.
Here is an example call::

    frbgalaxy.calc_nebular_SFR(method='Ha')

If successful, the 'SFR_nebular' key of the frbgalaxy.derived *dict*
will be filled with the value (using units of Msun/yr).
By default, an extinction correction will be applied to the measurement
if the 'AV_nebular' was filled previously.
This method also requires that the redshift have been set previously.

I/O
---

One can write the main contents of the FRBGalaxy object
to disk with the write_to_json() method::

    frbgalaxy.write_to_json()

If not outfile name is given (as above), one is auto-generated
based on the class and FRB.

Similarly, one can instantiate the class with a JSON file::

    frbgalaxy.from_json('name_of_file.json')

FRBHost
=======

This is a child of FRBGalaxy and is intended to be used
for the galaxy which hosts a given FRB.

by_name
-------

One particulaly useful method is by_name() which lets
you Instantiate the class by the FRB name::

    host121102 = frbgalaxy.FRBHost.by_name('121102')

Of course, a previously generated JSON file must already
have been archived.

PATH
----

A subset of the FRB Host galaxies have been analyzed using the
Probabilistic Assignment of Transients to their Hosts (PATH) framework
as described in `Aggarawal et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021ApJ...911...95A/abstract>`_.
This relies on the `astropath <>`_ code base.

You can load up the PATH results using::

    from frb.galaxies import utils
    path_table = utils.load_PATH()

The standard file uses the adopted priors with results
as given in the `adopted.csv` file
under frb/data/Galaxies/PATH.  See the README there
for further details.

Developers can use *frb_build* to generate new PATH csv files.