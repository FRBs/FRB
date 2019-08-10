#!/usr/bin/env python
"""
Build FRB bits and pieces
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Build parts of the CASBAH database; Output_dir = $CASBAH_GALAXIES [v1.1]')
    parser.add_argument("item", type=str, help="Item to build ['FRBs', 'Hosts', 'specDB', 'FG']")
    parser.add_argument("--flag", type=str, default='all', help="Flag passed to the build")

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
    item = pargs.item
    if item == 'FRBs':
        build_frbs.main(inflg=pargs.flag)
    elif item == 'Hosts':
        build_hosts.main(inflg=pargs.flag)
    elif item == 'specDB':
        build_specdb.main(inflg=pargs.flag)
    elif item == 'FG':
        build_fg.main(inflg=pargs.flag)
    else:
        raise IOError("Bad build item {:s}".format(item))


