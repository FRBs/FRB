
"""
Estimate p(z|DM) for an assumed location on the sky and DM_FRB 
as well as the limiting magnitude for the host galaxy
Defaults to using a perfect telescope model for the DM-z grid
"""
from IPython import embed


def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to print a summary of an FRB to the screen [v1.0]')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)")
    parser.add_argument("DM_FRB", type=float, help="FRB DM (pc/cm^3)")
    parser.add_argument("--mag_limit", type=float, default=20., help="Magnitude limit in filter *without* extinction correction. Default = 20")
    parser.add_argument("--filter", type=str, default='DECaL_r', help="Filter -- only used for extinction correction.  Must be a Repo approved choice")
    parser.add_argument("--dm_host", type=float, default=50., help="Assumed DM contribution from the Host. Default = 50")
    parser.add_argument("--dm_mwhalo", type=float, default=50., help="Assumed DM contribution from the MW halo. Default = 50")
    parser.add_argument("--cl", type=str, default="2.5,97.5", 
                        help="Confidence limits for the z estimate [default is a 95 percent c.l., (2.5,97.5)]")
    parser.add_argument("--telescope", type=str, default='perfect', help="telescope model for the DM-z grid: CHIME, DSA, Parkes, FAST, CRAFT, \
                        CRAFT_ICS_892/1300/1632, perfect. Default = perfect")
    parser.add_argument("--magdm_plot", default=False, action='store_true', 
                        help="Plot the host redshift range given DM on the magnitude vs redshift evolution")
    parser.add_argument("--fig_title", type=str,  help="title for the figure; e.g., FRBXXXXX")
    parser.add_argument("--fig_name", type=str, default='fig_r_vs_z.png', help="name of the output figure")
    parser.add_argument("--zmin", type=float, required=False,  help="Minimum redshift for the plot")
    parser.add_argument("--zmax", type=float, required=False, help="Maximum redshift for the plot")

    

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
    from frb.galaxies import nebular
    from frb.galaxies import photom
    from frb.galaxies import utils as frb_gal_u


    # Deal with coord
    icoord = ltu.radec_to_coord(coord_arg_to_coord(pargs.coord))

    # EBV
    EBV = nebular.get_ebv(icoord)['meanValue']  
    print("EBV = ", EBV)
   
    # NE 2001
    DM_ISM = mw.ismDM(icoord)
    print("")
    print("-----------------------------------------------------")
    print(f"NE2001 = {DM_ISM:.2f}")

    # DM cosmic and EG
    DM_extragalactic = pargs.DM_FRB - DM_ISM.value - pargs.dm_mwhalo
    DM_cosmic = DM_extragalactic - pargs.dm_host
     

    # Load the telescope specific grid
    telescope_dict = prob_dmz.telescope_dict

    # Get the perfect telescope grid (default)
    if not pargs.telescope or pargs.telescope == 'perfect':
        sdict = prob_dmz.grab_repo_grid(telescope_dict['perfect'])
        DM = sdict['DM']
        PDM_z = sdict['PDM_z']
        PDM_z = PDM_z.T # Perfect is opposite of the rest 
        iDM = np.argmin(np.abs(DM - DM_cosmic))
    else: # Grab a non-perfect telescope grid, which use DM_extragalactic
        sdict = prob_dmz.grab_repo_grid(telescope_dict[pargs.telescope])
        PDM_z = sdict['pzdm']
        DM = sdict['DM']
        iDM = np.argmin(np.abs(DM - DM_extragalactic))

    # Grab the right entry
    z = sdict['z']
    PzDM = PDM_z[:,iDM] / np.sum(PDM_z[:,iDM])
    
    cum_sum = np.cumsum(PzDM)
    limits = [float(item) for item in pargs.cl.split(',')]

    z_min = z[np.argmin(np.abs(cum_sum-limits[0]/100.))]
    z_max = z[np.argmin(np.abs(cum_sum-limits[1]/100.))]
    
    z_50 = z[np.argmin(np.abs(cum_sum-50./100.))]
    z_mode = z[np.argmax(PzDM)]


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
    z_min_capped = max([z_min,0.02])  # Capped at 0.02
    z_max_capped = max([z_max,0.02])  # Capped at 0.02
    m_r_Lstar_min = float(f_mL(z_min_capped))
    m_r_Lstar_max = float(f_mL(z_max_capped))

    frac_Lstar_min = 10**(-0.4*(mag_corr-m_r_Lstar_min))
    frac_Lstar_max = 10**(-0.4*(mag_corr-m_r_Lstar_max))


    # Finish
    print("")
    print(f"Allowing for the MW halo, DM_MW_halo = {int(pargs.dm_mwhalo)} pc/cm^3")
    print(f"Allowing for the Host, DM_host = {int(pargs.dm_host)} pc/cm^3")
    print("")
    print("")
    print(f"The mean redshift value is: {z_50:.3f}")
    print(f"The mode redshift value is: {z_mode:.3f}")
    print("")
    print(f"The redshift range for your confidence interval [{pargs.cl}] is:")
    print(f"z = [{z_min:.3f}, {z_max:.3f}]")
    print("")
    if not pargs.telescope or pargs.telescope == 'perfect':
        print("WARNING: This all assumes a perfect telescope and a model of the scatter in DM_cosmic (Macquart+2020)")
    else:
        print("This assumes the "+(str(pargs.telescope))+" telescope and a model of the scatter in DM_cosmic (Macquart+2020)")
    print("-----------------------------------------------------")

    print(f"For z_({limits[0]} %)={z_min:.2f}, the limiting magnitude corresponds to L={frac_Lstar_min:.5f}L*")
    print(f"For z_({limits[1]} %)={z_max:.2f}, the limiting magnitude corresponds to L={frac_Lstar_max:.5f}L*")

    # make the magnitude vs redshift plot with z-range if requested
    # Default values
    z_min_plot = z_min
    z_max_plot = z_max
    # If user provides values
    if pargs.zmin:
        z_min_plot = pargs.zmin
    if pargs.zmax:
        z_max_plot = pargs.zmax
    if pargs.magdm_plot:
        fig_name = pargs.fig_name
        mag_dm.r_vs_dm_figure(z_min_plot, z_max_plot, z, PzDM, outfile=fig_name,
               flipy=True, known_hosts=False, title=pargs.fig_title, logz_scale=False)


    return z_min, z_max, z_50, z_mode, frac_Lstar_min, frac_Lstar_max

