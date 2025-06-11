""" Calculate p(z|DM) for a given DM and survey"""

# It should be possible to remove all the matplotlib calls from this
# but in the current implementation it is not removed.
import argparse
import imp
import numpy as np
import os

from zdm import iteration as it
from zdm import io
from zdm.craco import loading

from IPython import embed

def main(pargs):
    
    import numpy as np
    from matplotlib import pyplot as plt

    from linetools import utils as ltu
    from linetools.scripts.utils import coord_arg_to_coord

    from frb import mw
    from frb.figures.utils import set_fontsize 

    from zdm import survey
    from zdm import parameters
    from zdm import cosmology as cos
    from zdm import misc_functions

    limits = (2.5, 97.5)

    # Deal with coord
    icoord = ltu.radec_to_coord(coord_arg_to_coord(pargs.coord))

    # NE 2001
    DM_ISM = mw.ismDM(icoord)
    print("")
    print("-----------------------------------------------------")
    print(f"NE2001 = {DM_ISM:.2f}")

    # DM EG
    DM_EG = pargs.DM_FRB - DM_ISM.value - pargs.dm_mwhalo

    # Cosmology
    state = parameters.State()
    cos.set_cosmology(state)
    cos.init_dist_measures()

    # Load Survey

    # get the grid of p(DM|z)
    zDMgrid, zvals,dmvals = misc_functions.get_zdm_grid(
        state, new=True, plot=False, method='analytic')

    # Suvey
    isurvey = survey.load_survey(pargs.survey, state, dmvals)

    # Grid
    igrid = misc_functions.initialise_grids(
        [isurvey], zDMgrid, zvals, dmvals, state, wdist=True)[0]
    PDM_z = igrid.rates # z, DM

    # Fuss
    iDM = np.argmin(np.abs(dmvals - DM_EG))
    PzDM = PDM_z[:, iDM] / np.sum(PDM_z[:, iDM])

    # Set zmax
    izmax = np.max(np.where(PzDM > 1e-10)[0])
    zmax = zvals[izmax]

    # Limits
    cum_sum = np.cumsum(PzDM)
    z_min = zvals[np.argmin(np.abs(cum_sum-limits[0]/100.))]
    z_max = zvals[np.argmin(np.abs(cum_sum-limits[1]/100.))]

    # Plot
    plt.clf()
    ax = plt.gca()
    ax.plot(zvals, PzDM)

    # Limits
    for z in [z_min, z_max]:
        ax.axvline(z, color='red', ls='--')

    ax.set_xlim(0, zmax)

    ax.set_xlabel('z')
    ax.set_ylabel('P(z|DM) [Normalized]')
    set_fontsize(ax, 15.)
    plt.show()

        

def parse_args(options=None):
    # test for command-line arguments here
    parser = argparse.ArgumentParser()
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)")
    parser.add_argument("DM_FRB", type=float, help="FRB DM (pc/cm^3)")
    parser.add_argument('-s','--survey',type=str, default='CRAFT/ICS', help="Name of survey [CRAFT/ICS, PKS/Mb]")
    parser.add_argument("--dm_mwhalo", type=float, default=40., help="Assumed DM contribution from the Milky Way Halo (ISM is calculated from NE2001). Default = 40")
    args = parser.parse_args()
    return args

def run():
    pargs = parse_args()
    main(pargs)

'''
# Test
python py/craco_H0_Emax_cube.py -n 1 -m 100 -o tmp.out --clobber

# 
python py/craco_H0_Emax_cube.py -n 1 -m 250 -o Cubes/craco_H0_Emax_cube0.out --clobber '''
