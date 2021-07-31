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
    parser.add_argument("--frb", type=str, help="FRB name, e.g. FRB191001, FRB20191001, 20191001")
    parser.add_argument("--hosts_file", type=str, help="Alternate file for hosts than the default public_hosts.csv")
    parser.add_argument("--lit_refs", type=str, help="Alternate file for literature sources than all_refs.csv")

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
    from frb.builds import new_build_hosts
    from frb.builds import build_fg

    # Parse
    item = pargs.item.lower()
    if item == 'frbs':
        build_frbs.main(inflg=pargs.flag)
    elif item == 'hosts':
        if pargs.frb is None:
            print("You must specify --frb")
            return
        # 
        frbs = pargs.frb.split(',')
        frbs = [ifrb.strip() for ifrb in frbs]
        new_build_hosts.main(frbs, options=pargs.galaxy_options, 
                             hosts_file=pargs.hosts_file,
                             lit_refs=pargs.lit_refs) 
    elif item == 'specdb':
        build_specdb.main(inflg=pargs.flag)
    elif item == 'fg':
        build_fg.main(inflg=pargs.flag, options=pargs.galaxy_options)
    else:
        raise IOError("Bad build item {:s}".format(item))


