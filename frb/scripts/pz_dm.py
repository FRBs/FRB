
"""
Estimate p(z|DM) for an assumed location on the sky and DM_FRB
Defaults to using a perfect telescope model for the DM-z grid
"""
from IPython import embed


def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to print a summary of an FRB to the screen [v1.0]')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)")
    parser.add_argument("DM_FRB", type=float, help="FRB DM (pc/cm^3)")
    parser.add_argument("--dm_host", type=float, default=50., help="Assumed DM contribution from the Host. Default = 50")
    parser.add_argument("--dm_mwhalo", type=float, default=50., help="Assumed DM contribution from the MW halo. Default = 50")
    parser.add_argument("--cl", type=tuple, default=(2.5,97.5), 
                        help="Confidence limits for the z estimate [default is a 95 percent c.l., (2.5,97.5)]")
    parser.add_argument("--magdm_plot", default=False, action='store_true', 
                        help="Plot the host redshift range given DM on the magnitude vs redshift evolution")
    parser.add_argument("--telescope", type=str, default='perfect', help="telescope model for the DM-z grid: CHIME, DSA, Parkes, FAST, CRAFT, \
                        CRAFT_ICS_892/1300/1632, perfect. Default = perfect")
    parser.add_argument("--fig_title", type=str,  help="title for the figure; e.g., FRBXXXXX")

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
    from frb.galaxies import mag_dm


    # Deal with coord
    icoord = ltu.radec_to_coord(coord_arg_to_coord(pargs.coord))
   
    # NE 2001
    DM_ISM = mw.ismDM(icoord)
    print("")
    print("-----------------------------------------------------")
    print(f"NE2001 = {DM_ISM:.2f}")

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
    
    z_50 = z[np.argmin(np.abs(cum_sum-50./100.))]
    z_mode = z[np.argmax(PzDM)]

    # Finish
    print("")
    print(f"Allowing for the MW halo, DM_MW_halo = {int(pargs.dm_mwhalo)} pc/cm^3")
    print(f"Allowing for the Host, DM_host = {int(pargs.dm_host)} pc/cm^3")
    print("")
    print("")
    print(f"The mean redshift value is: {z_50:.3f}")
    print(f"The mode redshift value is: {z_mode:.3f}")
    print("")
    print(f"The redshift range for your confidence interval {pargs.cl} is:")
    print(f"z = [{z_min:.3f}, {z_max:.3f}]")
    print("")
    if not pargs.telescope or pargs.telescope == 'perfect':
        print("WARNING: This all assumes a perfect telescope and a model of the scatter in DM_cosmic (Macqurt+2020)")
    else:
        print("This assumes the "+(str(pargs.telescope))+" telescope and a model of the scatter in DM_cosmic (Macqurt+2020)")
    print("-----------------------------------------------------")



    # make the magnitude vs redshift plot with z-range if requested
    if pargs.magdm_plot:
        mag_dm.r_vs_dm_figure(z_min, z_max, z, PzDM, outfile='fig_r_vs_z.png',
               flipy=True, known_hosts=False, title=pargs.fig_title, logz_scale=False)

    return z_min, z_max, z_50, z_mode

