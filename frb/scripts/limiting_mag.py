#!/usr/bin/env python
"""
Estimate the limiting luminosity given the magnitude limit and DMs and coord and telescope
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
    parser.add_argument("--dm_host", type=float, default=50., help="Assumed DM contribution from the Host. Default = 50")
    parser.add_argument("--dm_mwhalo", type=float, default=50., help="Assumed DM contribution from the MW halo. Default = 50")
    #parser.add_argument("-v", "--verbose", default=False, action="store_true", help="Overwhelm the screen?")
    parser.add_argument("--telescope", type=str, default='perfect', help="telescope model for the DM-z grid: CHIME, DSA, Parkes, FAST, CRAFT, \
                        CRAFT_ICS_892/1300/1632, perfect. Default = perfect")
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

    # DM cosmic and EG
    DM_extragalactic = pargs.DM_FRB - DM_ISM.value - pargs.dm_mwhalo
    DM_cosmic = DM_extragalactic - pargs.dm_host

    # Redshift estimates

    # Load the telescope specific grid
    telescope_dict = {
        'CHIME': 'CHIME_pzdm.npz',
        'DSA': 'DSA_pzdm.npy',
        'Parkes': 'parkes_mb_class_I_and_II_pzdm.npy',
        'CRAFT': 'CRAFT_class_I_and_II_pzdm.npy',
        'CRAFT_ICS_1300': 'CRAFT_ICS_1300_pzdm.npy',
        'CRAFT_ICS_892': 'CRAFT_ICS_892_pzdm.npy',
        'CRAFT_ICS_1632': 'CRAFT_ICS_1632_pzdm.npy',
        'FAST': 'FAST_pzdm.npy',
        'perfect': 'PDM_z.npz'
    }

    # Get the perfect telescope grid (default)
    sdict = prob_dmz.grab_repo_grid(telescope_dict['perfect'])
    PDM_z = sdict['PDM_z']
    z = sdict['z']
    DM = sdict['DM']

    # Grab the right entry
    iDM = np.argmin(np.abs(DM - DM_cosmic))
    PzDM = PDM_z[iDM, :] / np.sum(PDM_z[iDM, :])


    # Get the telescope specific PZDM grid
    if pargs.telescope and pargs.telescope != 'CHIME' and pargs.telescope != 'perfect':
        if pargs.telescope not in telescope_dict:
            raise ValueError(f"Unknown telescope: {pargs.telescope}")
        zdict = prob_dmz.grab_repo_grid(telescope_dict['CHIME'])
        z = zdict['z']
        DM = zdict['DM']
        PDM_z = prob_dmz.grab_repo_grid(telescope_dict[pargs.telescope])
        iDM = np.argmin(np.abs(DM - DM_extragalactic))
        PzDM = PDM_z[:,iDM] / np.sum(PDM_z[:,iDM])


    if pargs.telescope and pargs.telescope == 'CHIME':
        sdict = prob_dmz.grab_repo_grid(telescope_dict['CHIME'])
        PDM_z = sdict['pzdm']
        z = sdict['z']
        DM = sdict['DM']
        iDM = np.argmin(np.abs(DM - DM_extragalactic))
        PzDM = PDM_z[:,iDM] / np.sum(PDM_z[:,iDM])

    cum_sum = np.cumsum(PzDM)
    limits = pargs.cl

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
    print("")
    print(f"Allowing for the MW halo, DM_MW_halo = {int(pargs.dm_mwhalo)} pc/cm^3")
    print(f"Allowing for the Host, DM_host = {int(pargs.dm_host)} pc/cm^3")
    print("")
    if not pargs.telescope or pargs.telescope == 'perfect':
        print("WARNING: This all assumes a perfect telescope and a model of the scatter in DM_cosmic (Macquart+2020)")
    else:
        print("This assumes the "+(str(pargs.telescope))+" telescope and a model of the scatter in DM_cosmic (Macquart+2020)")
    print("")
    
    print(f"For z_({limits[0]} %)={z_min:.2f}, the limiting magnitude corresponds to L={frac_Lstar_min:.5f}L*")
    print(f"For z_({limits[1]} %)={z_max:.2f}, the limiting magnitude corresponds to L={frac_Lstar_max:.5f}L*")

    return frac_Lstar_min, frac_Lstar_max

# frb_mag_limit J151849.52+122235.8 200. 23.        
    