.. highlight:: rest

******************
Dispersion Measure
******************

This document describes
calculations related to the
Dispersion Measure (DM) of FRBs.

Galactic ISM
============

The repository contains a HEALPix map of the
Galactic ISM DM values based on the NE2001 model.
This map is used to calculate the DM contribution
from the Milky Way ISM.  The map is in the
`frb/data/DM/ne2001_dm_healpix_map_128.fits` file.

Here is an example usage:

```python
from frb.dm import dm_ism_healpix_map

b, l = 5., 50.  # Galactic coordinates in degrees

dm_map = dm_ism_healpix_map.get_dm_map()
dm_ism = dm_ism_healpix_map.dm_ism_from_healpix_map(l, b, dm_map)

```

One can also input arrays of b,l coordinates.


Intervening Galaxies
====================

A set of estimates for the DM from intervening
galaxies is provided.  These are based on the
analysis given in `PN17`_.  All of the analysis
is based on statistics of the damped Lya systems
(DLAs).

Approximate
-----------

An approximate average DM to a given
redshift may be calculated with the
frb.dlas.approx_avgDM() method::

    DM = approx_avgDM(1.)

This may be calucated a single or array of redshifts.
The return value is an astropy Quantity with default
unites of pc/cm^3.

.. _PN17: http://coming.soon

MonteCarlo
----------

This algorithm generates realizations of the DM values
by randomly placing intervening galaxies along a
given sightline or sightlines.  Here is an example call::

    zeval = np.array([0.,1.,2.])
    DMs = monte_DM(np.array(zeval), nrand=1000)  # Returned without units

The values returned are in units of pc/cm^3 but
without astropy units attached.
One can then perform stats on these outputs to the
heart's desire.
