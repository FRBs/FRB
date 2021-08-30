********************
Adding a Host Galaxy
********************

Overview
========

This document describes how to add a host galaxy
to the Repo.

Steps
=====

Basics
------

#. Add the FRB to the Repository (JSON file)
#. Measure the galaxy redshift, RA, DEC
#. Add the Host to the Google sheet [optional]
#. Download the Google sheet as CSV [optional]
#. Move/edit the public_hosts.csv file in frb/data/Galaxies
#. Run frb_build (e.g. *frb_build Hosts --frb 20201123*)

Literature
----------

ppxf
----

#. Add spectrum to the specdb [challenging]

   #. Add spectrum to the Galaxy_DB with name JXXXXXXX.XX+XXXXXX.XX_INSTR_A_spec.fits

   #. Add redshift to the z_hand.ascii table 
   
   #. Edit + Run build_specdb.py 

#. Add "cuts" to the Google Sheet (regions of the spectrum to avoid)
#. Add spectrum folder to Projects and References
#. Add instrument name to the Google sheet (Spectrum column)
#. Download the Google sheet as CSV [optional]
#. Move/edit the public_hosts.csv file in frb/data/Galaxies