.. highlight:: rest

****************
Rotation Measure
****************

This document describes
calculations related to the
Rotation Measure (RM) of FRBs.

Galactic RM
===========

`Oppermann et al. (2014) <https://arxiv.org/abs/1404.3701>`_ have generated an all-sky
map of RM estimates and have provided a public Healpix description of it.
We have archived a copy and have added a simple wrapper to generate
RM estimates from it.  Here is a call to the main method::

    from astropy.coordinates import SkyCoord
    from frb import rm
    repeater_coord = SkyCoord('05h31m58.698s +33d8m52.59s', frame='icrs')
    #
    RM, RM_err = rm.galactic_rm(repeater_coord)

The value and uncertainty are returned as astropy.units.Quantity objects
with units of rad/m^2.
