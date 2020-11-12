#!/usr/bin/env python
"""
Script to print a summary of a given FRB to the terminal
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to print a summary of an FRB to the screen [v1.0]')
    parser.add_argument("frb_name", type=str, help="FRB name, e.g. FRB180924 or simply 180924")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    """
    import warnings
    import numpy as np
    import json

    from linetools import utils as ltu

    from frb import frb

    # Load
    frb_name = pargs.frb_name if pargs.frb_name[0] == 'F' else 'FRB'+pargs.frb_name
    FRB = frb.FRB.by_name(frb_name)

    def pjson(obj):
        return json.dumps(obj, sort_keys=True, indent=4,
                         separators=(',', ': '))

    # Coords
    print(frb_name)
    print(ltu.name_from_coord(FRB.coord))
    print('ee={}'.format(pjson(FRB.eellipse)))
    print('DM={}'.format(FRB.DM))

    # Host
    hg = FRB.grab_host()
    if hg is not None:
        print("Host")
        print(ltu.name_from_coord(hg.coord))
        print('z={}'.format(pjson(hg.redshift)))
        # photometry
        print('photom:{}'.format(pjson(hg.photom)))

