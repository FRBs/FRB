#!/usr/bin/env python
"""
Script to print a summary of a given FRB to the terminal
"""
from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to print a summary of an FRB to the screen [v1.0]')
    parser.add_argument("frb_name", type=str, help="FRB name, e.g. FRB180924 or simply 180924")
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

    from frb import frb

    # Load
    frb_name = pargs.frb_name if pargs.frb_name[:3].upper() == 'FRB' else 'FRB'+pargs.frb_name
    try:
        FRB = frb.FRB.by_name(frb_name.upper())
    except FileNotFoundError:
        print("Your desired FRB is not in this Repo.  Sorry..")
        return

    def pjson(obj):
        return json.dumps(obj, sort_keys=True, indent=4,
                         separators=(',', ': '))

    # Coords
    print("=========================================================\n")
    print(frb_name)
    print(ltu.name_from_coord(FRB.coord, precision=(4,4)), '\n     ra,dec = {:0.7f},{:0.7f} deg'.format(FRB.coord.ra.value,
                                                                           FRB.coord.dec.value))
    print('sig_a = {:0.3f} arcsec'.format(FRB.sig_a))
    print('sig_b = {:0.3f} arcsec'.format(FRB.sig_b))
    print('PA = {:0.1f} deg'.format(FRB.eellipse['theta']))
    print('')
    print('DM={}'.format(FRB.DM))

    print('Repeater={}'.format(FRB.repeater))

    if pargs.verbose:
        print('ee={}'.format(pjson(FRB.eellipse)))
        print('DM_ISM(ne2001)={:0.1f}'.format(FRB.DMISM))


    # Host
    hg = FRB.grab_host()
    if hg is not None:
        print("")
        print("=========================================================\n")
        print("Host\n")
        print(ltu.name_from_coord(hg.coord))
        print('z: \n {}'.format(pjson(hg.redshift)))
        if pargs.verbose:
            # photometry
            print('photom: \n {}'.format(pjson(hg.photom)))
            print('derived: \n {}'.format(pjson(hg.derived)))
            print('morphology: \n {}'.format(pjson(hg.morphology)))
            print('neb_lines: \n {}'.format(pjson(hg.neb_lines)))


