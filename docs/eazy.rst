****
Eazy
****

Overview
========

This document describes how to run the Eazy code.

Setup
=====

Before beginning, you will need to have installed EAZY
as per the :doc:`installing` instructions.

EAZY runs on a set of photometric data.  These need to be
provides as a `dict`, preferentially in an FRBGalaxy object.

Here is an example of creating the object::


    from astropy.table import Table
    from frb.galaxies.frbgalaxy import FRBHost

    photom = Table()
    photom['Name'] = ['G_TEST']
    photom['ra'] = 123.422
    photom['dec'] = 23.222
    # These are observed
    photom['LRISb_V'] = 25.86
    photom['LRISb_V_err'] = 0.25
    photom['GMOS_S_r'] = 23.61
    photom['GMOS_S_r_err'] = 0.15
    photom['LRISr_I'] = 23.09
    photom['LRISr_I_err'] = 0.1
    photom['NOT_z'] = 23.35
    photom['NOT_z_err'] = 0.3
    photom['NIRI_J'] = 21.75 + 0.91
    photom['NIRI_J_err'] = 0.2

    #
    host = FRBHost(photom['ra'], photom['dec'], photom['Name'])
    host.parse_photom(photom)
    host.name = 'G_TEST'

Now you are ready to generate the input files.  Here is
a standard call::

    from frb.galaxies import eazy as frbeazy
    input_dir = './eazy_running/inputs'  # Sets where the input files land
    output_dir = './eazy_running/outputs'  # Sets where the output files land
    frbeazy.eazy_input_files(host.photom, input_dir,
                             host.name, output_dir
                             templates='br07_default',
                             prior_filter='GMOS_S_r')

See the doc string of `eazy_input_files()` for all the options.

If run successfully, you will have 3 files in the `input_dir`
and a softlink above that folder to the `templates` folder of EAZY.

Running
=======

You are now ready to run EAZY.  Here is a standard call::

    frbeazy.run_eazy(input_dir, host.name,
                     os.path.join(output_dir 'logfile'))

This will run EAZY and place the outputs in the `output_dir`.

Parsing
=======

You can read in the outputs using scripts taken from
threedhst by Brummer.  Here is an example::

    zgrid, pzi, prior = frbeazy.getEazyPz(-1, MAIN_OUTPUT_FILE='photz',
                                          OUTPUT_DIRECTORY=out_dir,
                                          CACHE_FILE='Same', binaries=None, get_prior=True)
    zphot, sig_zphot = frbeazy.eazy_stats(zgrid, pzi)

