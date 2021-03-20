#!/usr/bin/env python
"""
Build FRB bits and pieces
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Build parts of the CASBAH database; Output_dir = $CASBAH_GALAXIES [v1.1]')
    parser.add_argument("item", type=str, help="Item to build ['FRBs', 'Hosts', 'specDB', 'FG']. Case insensitive")
    parser.add_argument("--flag", type=str, default='all', help="Flag passed to the build")
    parser.add_argument("-g", "--galaxy_options", type=str, help="Options for fg/host building (photom,cigale,ppxf)")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    """
    import warnings
    from frb.builds import build_specdb
    from frb.builds import build_frbs
    from frb.builds import build_hosts
    from frb.builds import build_fg

    # Parse
    item = pargs.item.lower()
    if item == 'frbs':
        build_frbs.main(inflg=pargs.flag)
    elif item == 'hosts':
        build_hosts.main(inflg=pargs.flag, options=pargs.galaxy_options)
    elif item == 'specdb':
        build_specdb.main(inflg=pargs.flag)
    elif item == 'fg':
        build_fg.main(inflg=pargs.flag, options=pargs.galaxy_options)
    else:
        raise IOError("Bad build item {:s}".format(item))


