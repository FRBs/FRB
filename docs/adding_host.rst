********************
Adding a Host Galaxy
********************

Overview
========

This document describes how to add a host galaxy
to the Repo.

Steps
=====

1. Add the FRB to the Repository (JSON file)
1. Measure the galaxy redshift, RA, DEC
1. Add the Host to the Google sheet [optional]
1. Download the Google sheet as CSV [optional]
1. Move/edit the public_hosts.csv file in frb/data/Galaxies
1. Run frb_build (e.g. *frb_build Hosts --frb 20201123*)

ppxf
----

1. Add spectrum to the specdb [challenging]
   1. Add spectrum to the Galaxy_DB with name JXXXXXXX.XX+XXXXXX.XX_INSTR_A_spec.fits
   1. Add redshift to the z_hand.ascii table 
   1. Edit + Run build_specdb.py 
1. Add "cuts" to the Google Sheet (regions of the spectrum to avoid)
1. Add spectrum folder to Projects and References
1. Add instrument name to the Google sheet (Spectrum column)
1. Download the Google sheet as CSV [optional]
1. Move/edit the public_hosts.csv file in frb/data/Galaxies