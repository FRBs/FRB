#!/usr/bin/env python
"""
Build FRB bits and pieces
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Build parts of the CASBAH database; Output_dir = $CASBAH_GALAXIES [v1.1]')
    parser.add_argument("item", type=str, help="Item to build ['FRBs', 'Hosts', 'specDB', 'FG', 'PATH']. Case insensitive")
    parser.add_argument("--flag", type=str, default='all', help="Flag passed to the build")
    parser.add_argument("--options", type=str, help="Options for the build, e.g. fg/host building (cigale,ppxf)")
    parser.add_argument("--frb", type=str, help="FRB name, e.g. FRB191001, FRB20191001, 20191001")
    parser.add_argument("--data_file", type=str, help="Alternate file for data than the default (public)")
    parser.add_argument("--lit_refs", type=str, help="Alternate file for literature sources than all_refs.csv")
    parser.add_argument("--override", default=False, action='store_true',
                        help="Over-ride errors (as possible)? Not recommended")

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
    from frb.builds import build_path

    # Parse
    item = pargs.item.lower()
    if item in ['frbs', 'hosts']:
        if pargs.frb is None:
            print("You must specify --frb")
            return
        # 
        frbs = pargs.frb.split(',')
        frbs = [ifrb.strip() for ifrb in frbs]
        if item == 'frbs':
            build_frbs.main(frbs, data_file=pargs.data_file)
        else:
            build_hosts.main(frbs, options=pargs.options, 
                             hosts_file=pargs.data_file,
                             lit_refs=pargs.lit_refs,
                             override=pargs.override) 
    elif item == 'specdb':
        build_specdb.main(inflg=pargs.flag)
    elif item == 'fg':
        build_fg.main(inflg=pargs.flag, options=pargs.options)
    elif item == 'path':
        build_path.main(options=pargs.options)
    else:
        raise IOError("Bad build item {:s}".format(item))


