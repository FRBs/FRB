*******
Scripts
*******

Overview
========

This document describes scripts that come with the Repo.

frb_summary
===========

This script prints a simple summary of a given FRB and its
host galaxy (when that exists) to the screen.

Here is the usage::

    usage: frb_summary [-h] [--verbose] frb_name

    Script to print a summary of an FRB to the screen [v1.0]

    positional arguments:
      frb_name    FRB name, e.g. FRB180924 or simply 180924

    optional arguments:
      -h, --help  show this help message

Here is an example::

    frb_summary 180924

    FRB180924
    J214425.26-405400.1
    ee={
        "a": 0.07,
        "a_sys": 0.09,
        "b": 0.06,
        "b_sys": 0.07,
        "cl": 68.0,
        "cl_sys": 68.0,
        "theta": 0.0,
        "theta_sys": 0.0
    }
    DM=362.16 pc / cm3
    =========================================================

    Host

    J214425.25-405400.8
    z:
     {
        "z": 0.3212,
        "z_FRB": 0.3212,
        "z_spec": 0.3212
    }

Enjoy

frb_mag_limit
=============

This script takes as input the FRB DM and a magnitude limit for 
the data and then estimates the redshift range assuming the Macquart relation
(and simple host + MW contributions).  It then converts the 
magnitude into fraction of L* for the host.

Here is the usage::

    usage: frb_mag_limit [-h] [--filter FILTER] [--dm_hostmw DM_HOSTMW]
                        coord DM_FRB mag_limit

    Script to print a summary of an FRB to the screen [v1.0]

    positional arguments:
    coord                 Coordinates, e.g. J081240.7+320809 or
                            122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB
                            name (FRB180924)
    DM_FRB                FRB DM
    mag_limit             Magnitude limit in filter *without* extinction
                            correction

    optional arguments:
    -h, --help            show this help message and exit
    --filter FILTER       Filter -- only used for extinction correction. Must
                            be a Repo approved choice
    --dm_hostmw DM_HOSTMW
                            Assumed DM contribution from MW and Host


And an example::

    frb_mag_limit J151849.52+122235.8 200. 23. 
    EBV = 0.0366
    NE2001 = 27.019541038348713 pc / cm3
    Loading P(DM,z) from /data/Projects/FRB_Software/FRB/frb/data/DM/PDM_z.npz
    -----------------------------------------------------
    For z_10=0.06, the limiting magnitude corresponds to L=0.00152L*
    For z_90=0.14, the limiting magnitude corresponds to L=0.01805L*

The first time you run this script, it will generate a npz file and place
it in the repo for future use.  That takes a few minutes.

frb_pz_dm
=========

This script takes as input the FRB DM and its coordinates (approximate
are fine) and then estimates the redshift range assuming 
the Macquart relation (and simple host + MW contributions, optionally 
input).  For an (optionally input; tuple) confidence interval, 
it reports back the putative redshift range for the FRB.

Here is the usage::

    usage: frb_pz_dm [-h] [--dm_hostmw DM_HOSTMW] [--cl CL] coord DM_FRB

    Script to print a summary of an FRB to the screen [v1.0]

    positional arguments:
    coord                 Coordinates, e.g. J081240.7+320809 or
                            122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB
                            name (FRB180924)
    DM_FRB                FRB DM (pc/cm^3)

    optional arguments:
    -h, --help            show this help message and exit
    --dm_hostmw DM_HOSTMW
                            Assumed DM contribution from the Milky Way Halo (ISM
                            is calculated from NE2001) and Host. Default = 100
    --cl CL               Confidence limits for the z estimate [default is a 95
                            percent c.l., (2.5,97.5)]


frb_sightline
=============

Simple script to derive a few items along a given sightline
including a listing of the public surveys covering that location.  
Input is the coordinates.  Here is the usage::

    usage: frb_sightline [-h] [-v] coord

    Script to print a summary of an FRB to the screen [v1.0]

    positional arguments:
    coord          Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or
                    07:45:00.47,34:17:31.1 or FRB name (FRB180924)

    optional arguments:
    -h, --help     show this help message and exit
    -v, --verbose  Overwhelm the screen?
