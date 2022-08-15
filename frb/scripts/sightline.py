#!/usr/bin/env python
"""
Script to print a summary of a given FRB to the terminal
"""
from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to print a summary of an FRB to the screen [v1.0]')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)")
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="Overwhelm the screen?")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    """
    import json

    from linetools import utils as ltu
    from linetools.scripts.utils import coord_arg_to_coord

    from frb.galaxies import nebular
    from frb import mw
    from frb.surveys import survey_utils
    from frb import frb

    # Deal with coord
    if 'FRB' in pargs.coord:
        FRB = frb.FRB.by_name(pargs.coord)
        icoord = FRB.coord
    else:
        icoord = ltu.radec_to_coord(coord_arg_to_coord(pargs.coord))

    # EBV
    EBV = nebular.get_ebv(icoord)['meanValue']  #
    AV = EBV * 3.1  # RV

    print("AV = {}".format(AV))

    # NE 2001
    DM_ISM = mw.ismDM(icoord)
    print(f"NE2001 = {DM_ISM}")

    # Surveys
    print("Checking the imaging surveys...")
    inside = survey_utils.in_which_survey(icoord)
    print(inside)