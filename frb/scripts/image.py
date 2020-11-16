#!/usr/bin/env python
"""
Script generate an image of an FRB
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to make a quick image figure [v1.0]')
    parser.add_argument("fits_file", type=str, help="Image FITS file with WCS")
    parser.add_argument("frb_coord", type=str, help="FRB Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)")
    parser.add_argument("--imsize", default=30., type=float, help="Image size in arcsec")
    parser.add_argument("--vmnx", type=str, help="Image scale:  vmin,vmax")
    parser.add_argument("--outfile", default='image.png', type=str, help="Output filename")

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
    from matplotlib import pyplot as plt
    from astropy.io import fits

    from astropy import units

    from frb import frb
    from frb.figures import galaxies as ffgal
    from frb.figures import utils as ffutils

    from linetools.scripts.utils import coord_arg_to_coord

    # Load up
    hdu = fits.open(pargs.fits_file)
    icoord = coord_arg_to_coord(pargs.frb_coord)

    # Parse
    if pargs.vmnx is not None:
        tstr = pargs.vmnx.split(',')
        vmnx = (float(tstr[0]), float(tstr[1]))
    else:
        vmnx = (None,None)

    # Dummy FRB object
    FRB = frb.FRB('TMP', icoord, 0.)
    FRB.set_ee(1.0, 1.0, 0., 95.)

    fig = plt.figure(figsize=(7, 7))
    ffutils.set_mplrc()
    ffgal.sub_image(fig, hdu, FRB, vmnx=vmnx, cmap='gist_heat',
                    frb_clr='white', imsize=pargs.imsize) #img_center=HG190608.coord,

    # Layout and save
    plt.tight_layout(pad=0.2, h_pad=0.1, w_pad=0.1)
    plt.savefig(pargs.outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(pargs.outfile))
