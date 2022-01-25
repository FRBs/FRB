#!/usr/bin/env python
"""
Estimate p(z|DM) for an assumed location on the sky and DM_FRB
"""
from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to print a summary of an FRB to the screen [v1.0]')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)")
    parser.add_argument("DM_FRB", type=float, help="FRB DM (pc/cm^3)")
    parser.add_argument("--dm_hostmw", type=float, default=100., help="Assumed DM contribution from the Milky Way Halo (ISM is calculated from NE2001) and Host. Default = 100")
    parser.add_argument("--cl", type=tuple, default=(2.5,97.5), 
                        help="Confidence limits for the z estimate [default is a 95 percent c.l., (2.5,97.5)]")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    """
    import numpy as np

    from linetools import utils as ltu
    from linetools.scripts.utils import coord_arg_to_coord

    from frb import mw
    from frb.dm import prob_dmz

    # Deal with coord
    icoord = ltu.radec_to_coord(coord_arg_to_coord(pargs.coord))

    # NE 2001
    DM_ISM = mw.ismDM(icoord)
    print("")
    print("-----------------------------------------------------")
    print(f"NE2001 = {DM_ISM:.2f}")

    # DM Cosmic
    DM_cosmic = pargs.DM_FRB - DM_ISM.value - pargs.dm_hostmw

    # Redshift estimates

    # Load
    sdict = prob_dmz.grab_repo_grid()
    PDM_z = sdict['PDM_z']
    z = sdict['z']
    DM = sdict['DM']

    # Do it
    iDM = np.argmin(np.abs(DM - DM_cosmic))
    PzDM = PDM_z[iDM, :] / np.sum(PDM_z[iDM, :])

    cum_sum = np.cumsum(PzDM)
    limits = pargs.cl

    z_min = z[np.argmin(np.abs(cum_sum-limits[0]/100.))]
    z_max = z[np.argmin(np.abs(cum_sum-limits[1]/100.))]

    # Finish
    print("")
    print(f"The redshift range for your confidence interval {pargs.cl} is:")
    print(f"z = [{z_min:.3f}, {z_max:.3f}]")

    return z_min, z_max
