*******
Scripts
*******

Overview
========

This document describes scripts that come with the Repo.
All scripts are installed to the ``bin/`` directory.

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

frb_pzdm_mag
=============

This script takes as input the FRB DM and its coordinates (approximate
are fine) and then estimates the redshift range assuming
the Macquart relation (and host + MW halo contributions, optionally
input).  For an (optionally input; tuple) confidence interval,
it reports back the putative redshift range for the FRB. It also
allows for plotting the host redshift range on the magnitude vs redshift
evolution and setting a title for the figure.  These calculations can be
done assuming a few different telescope models (CHIME, DSA, Parkes, FAST,
CRAFT, CRAFT_ICS_892/1300/1632) or a perfect telescope model (default).
The telescope models are used to determine the DM-z grids that have been
computed with the zdm code/repository.

Here is the usage::

    usage: frb_pzdm_mag [-h] [--mag_limit MAG_LIMIT] [--filter FILTER]
                        [--dm_host DM_HOST] [--dm_mwhalo DM_MWHALO]
                        [--cl CL] [--telescope TELESCOPE]
                        [--magdm_plot] [--fig_title FIG_TITLE]
                        [--fig_name FIG_NAME] [--zmin ZMIN] [--zmax ZMAX]
                        coord DM_FRB

    positional arguments:
    coord                 Coordinates, e.g. J081240.7+320809 or
                            122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB
                            name (FRB180924)
    DM_FRB                FRB DM (pc/cm^3)

    optional arguments:
    -h, --help            show this help message and exit
    --mag_limit MAG_LIMIT
                            Magnitude limit without extinction correction.
                            Default = 20
    --filter FILTER       Filter for extinction correction. Must be a Repo
                            approved choice. Default = DECaL_r
    --dm_host DM_HOST     Assumed DM contribution from the Host. Default = 50
    --dm_mwhalo DM_MWHALO
                            Assumed DM contribution from the MW halo. Default = 50
    --cl CL               Confidence limits for the z estimate [default is a 95
                            percent c.l., (2.5,97.5)]
    --telescope TELESCOPE
                            telescope model for the DM-z grid: CHIME, DSA, Parkes,
                            FAST, CRAFT, CRAFT_ICS_892/1300/1632, perfect.
                            Default = perfect
    --magdm_plot          Plot the host redshift range given DM on the magnitude
                            vs redshift evolution. Default=False.
    --fig_title FIG_TITLE
                            title for the figure; e.g., FRBXXXXX
    --fig_name FIG_NAME   name of the output figure. Default = fig_r_vs_z.png


frb_sightline
=============

Simple script to derive a few items along a given sightline
including a listing of the public surveys covering that location.
Input is the coordinates.  Here is the usage::

    usage: frb_sightline [-h] [-v] coord

    positional arguments:
    coord          Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or
                    07:45:00.47,34:17:31.1 or FRB name (FRB180924)

    optional arguments:
    -h, --help     show this help message and exit
    -v, --verbose  Overwhelm the screen?

frb_build
=========

Build FRB data products including FRB JSON files, Host JSON files,
specDB, foreground galaxy data, and PATH association results::

    usage: frb_build [-h] [--flag FLAG] [--options OPTIONS] [--frb FRB]
                     [--data_file DATA_FILE] [--lit_refs LIT_REFS]
                     [--override]
                     item

    positional arguments:
    item                  Item to build ['FRBs', 'Hosts', 'specDB', 'FG', 'PATH'].
                            Case insensitive

    optional arguments:
    --flag FLAG           Flag passed to the build
    --options OPTIONS     Options for the build, e.g. fg/host building
                            (cigale,ppxf); PATH (write_indiv)
    --frb FRB             Full TNS FRB name, e.g. FRB20191001A
    --data_file DATA_FILE
                            Alternate file for data than the default (public)
    --lit_refs LIT_REFS   Alternate file for literature sources than all_refs.csv
    --override            Over-ride errors (as possible)? Not recommended

frb_galaxies
=============

Script to access FRB galaxy data and spectra from the specDB archive.
Here is the usage::

    usage: frb_galaxies [-h] [--rho RHO] [--ang_offset ANG_OFFSET] [--cat]
                    [--specdb SPECDB] [-p] [--dust]
                    coord

    positional arguments:
      coord                 Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322
                            or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)

    optional arguments:
      --rho RHO             Maximum impact parameter in kpc [default=300.]
      --ang_offset ANG_OFFSET
                            Maximum offset in arcsec [over-rides --rho if set]
      --cat                 Only show data from the catalog (not meta)
      --specdb SPECDB       specDB file; defaults to $SPECDB/FRB_specdb.hdf5
      -p, --plot            Launch a plotting GUI?
      --dust                Dust correct the spectrum?

Here is an example call::

    frb_galaxies FRB180924

frb_image
=========

Script to make a quick image figure from a FITS file with WCS::

    usage: frb_image [-h] [--imsize IMSIZE] [--vmnx VMNX] [--outfile OUTFILE]
                     fits_file frb_coord

    positional arguments:
    fits_file             Image FITS file with WCS
    frb_coord             FRB Coordinates, e.g. J081240.7+320809 or FRB name

    optional arguments:
    --imsize IMSIZE       Image size in arcsec [default=30]
    --vmnx VMNX          Image scale: vmin,vmax
    --outfile OUTFILE     Output filename [default=image.png]

frb_dmism
=========

Script for DM ISM HEALPix map operations including
generating maps, querying DM_ISM values, and plotting.

frb_macquart
============

Script to collate F4 Watchlist CSV tables and generate a Macquart-style
DM versus redshift plot. It computes DM_ISM, DM_EG, and DM_cosmic_est for
each FRB, writes a collated CSV table, and saves the figure.

Here is the usage::

        usage: frb_macquart [-h] [--table-loc TABLE_LOC] [--outfile OUTFILE]
                                                [--force-rebuild] [--fig-outfile FIG_OUTFILE]
                                                [--color COLOR] [--emoji-file EMOJI_FILE]
                                                [--emoji-zoom EMOJI_ZOOM] [--font-size FONT_SIZE]
                                                [--special-frbs [SPECIAL_FRBS ...]]
                                                [--special-color SPECIAL_COLOR]
                                                [--plot-dm-eg | --no-plot-dm-eg]
                                                [--plot-running-mean | --no-plot-running-mean]
                                                [--show-plot | --no-show-plot]

        optional arguments:
            -h, --help            show this help message and exit
            --table-loc TABLE_LOC
                                                        Directory containing watchlist CSVs named like
                                                        F4_watchlist*.csv (default: ./)
            --outfile OUTFILE     Path to collated watchlist table CSV
                                                        (default: watchlist_collated.csv)
            --force-rebuild       Rebuild collated table even if --outfile exists
            --fig-outfile FIG_OUTFILE
                                                        Output filename for figure
                                                        (default: fig_macquart_with_script.png)
            --color COLOR         Scatter color for FRB points (default: pink)
            --emoji-file EMOJI_FILE
                                                        Path to emoji image for custom point markers
            --emoji-zoom EMOJI_ZOOM
                                                        Zoom level for emoji markers (default: 0.02)
            --font-size FONT_SIZE
                                                        Final axis/legend font size (default: 17)
            --special-frbs [SPECIAL_FRBS ...]
                                                        Optional list of FRB names to highlight
            --special-color SPECIAL_COLOR
                                                        Color for highlighted FRBs (default: blue)
            --plot-dm-eg | --no-plot-dm-eg
                                                        Plot DM_EG instead of DM_cosmic_est (default: True)
            --plot-running-mean | --no-plot-running-mean
                                                        Overplot binned running mean with errors
                                                        (default: True)
            --show-plot | --no-show-plot
                                                        Display figure interactively (default: False)

Here is an example call::

        frb_macquart --table-loc ../DESI/FRB_tables/ \
            --outfile watchlist_collated_with_script.csv \
            --fig-outfile fig_macquart_with_script.png \
            --color pink --emoji-file Cherry-Blossom-Emoji.png \
            --emoji-zoom 0.02 --plot-dm-eg --plot-running-mean --show-plot

frb_search_for_halos
====================

Script to search for foreground halos along an FRB sightline.

frb_tns
=======

Script for querying the Transient Name Server (TNS) for FRB entries.

