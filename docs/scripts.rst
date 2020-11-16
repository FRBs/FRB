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
