#!/usr/bin/env python
"""
Build FRB bits and pieces
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Build parts of the CASBAH database; Output_dir = $CASBAH_GALAXIES [v1.1]')
    parser.add_argument("item", type=str, help="Item to build ['FRBs', 'Hosts', 'specdb']")
    parser.add_argument("--flag", default='all', action="store_true", help="Flag passed to the build")

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

    # Parse
    item = pargs.item
    if item == 'FRBs':
        build_frbs.main(flg=pargs.flag)
    elif item == 'specdb':
        v = '12'  # Current version
        pargs = mk_specdb.parser([v])
        mk_specdb.main(pargs)
    else:
        raise IOError("Bad build item {:s}".format(item))


