.. highlight:: rest

******************
Dispersion Measure
******************

This document calculations related to the
Dispersion Measure (DM) of FRBs

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
unites of pc/cm^2.

.. _PN17: http://coming.soon
