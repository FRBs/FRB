#!/usr/bin/env python
"""
Estimate the limiting luminosity given the magnitude limit and DM and coord
"""
from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to print a summary of an FRB to the screen [v1.0]')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)")
    parser.add_argument("DM_FRB", type=float, help="FRB DM")
    parser.add_argument("mag_limit", type=float, help="Magnitude limit in filter *without* extinction correction")
    parser.add_argument("--filter", type=str, default='DECaL_r', help="Filter -- only used for extinction correction.  Must be a Repo approved choice")
    parser.add_argument("--dm_hostmw", type=float, default=100., help="Assumed DM contribution from MW and Host")
    #parser.add_argument("-v", "--verbose", default=False, action="store_true", help="Overwhelm the screen?")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    """
    import json
    import os
    import numpy as np
    from pkg_resources import resource_filename

    from linetools import utils as ltu
    from linetools.scripts.utils import coord_arg_to_coord

    from frb.galaxies import nebular
    from frb.galaxies import photom
    from frb.galaxies import utils as frb_gal_u
    from frb import mw
    from frb.dm import prob_dmz

    # Deal with coord
    icoord = ltu.radec_to_coord(coord_arg_to_coord(pargs.coord))

    # EBV
    EBV = nebular.get_ebv(icoord)['meanValue']  #
    print(f"EBV = {EBV}")

    # NE 2001
    DM_ISM = mw.ismDM(icoord)
    print(f"NE2001 = {DM_ISM}")

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
    limits = [10, 90]

    z_min = z[np.argmin(np.abs(cum_sum-limits[0]/100.))]
    z_max = z[np.argmin(np.abs(cum_sum-limits[1]/100.))]

    # Setup Luminosity

    # Extinction correct
    dust_correct = photom.extinction_correction(pargs.filter, EBV)
    mag_dust = 2.5 * np.log10(1. / dust_correct)
    mag_corr = pargs.mag_limit + mag_dust

    # ##########################3
    # Convert to L

    # Load f_mL
    f_mL = frb_gal_u.load_f_mL()
    # m_r(L*)
    m_r_Lstar_min = float(f_mL(z_min))
    m_r_Lstar_max = float(f_mL(z_max))

    frac_Lstar_min = 10**(-0.4*(mag_corr-m_r_Lstar_min))
    frac_Lstar_max = 10**(-0.4*(mag_corr-m_r_Lstar_max))

    # Finish
    print("-----------------------------------------------------")
    print(f"For z_{limits[0]}={z_min:.2f}, the limiting magnitude corresponds to L={frac_Lstar_min:.5f}L*")
    print(f"For z_{limits[1]}={z_max:.2f}, the limiting magnitude corresponds to L={frac_Lstar_max:.5f}L*")

    return frac_Lstar_min, frac_Lstar_max

# frb_mag_limit J151849.52+122235.8 200. 23.        
    