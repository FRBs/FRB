""" Estimate magnitude range for a host FRB given DM """

import numpy as np
import glob, os, sys, json
import pdb

from scipy.interpolate import interp1d

#import healpy as hp

from cycler import cycler
import matplotlib as mpl
import seaborn as sns

import pandas

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.coordinates import Distance
from astropy import visualization as vis
from astropy.coordinates import SkyCoord

from linetools import utils as ltu
from linetools.scripts.utils import coord_arg_to_coord

from frb.frb import FRB, build_table_of_frbs, list_of_frbs
from frb.figures import utils as frb_fig_u
from frb.figures import dm as frb_fig_dm
from frb.galaxies import utils as frb_gal_u
from frb import mw
from frb.dm import prob_dmz, igm

from IPython import embed

# user parameters
FRB = 'FRB20240304'
frb_coords = SkyCoord(ra=182.9972083, dec=11.8130667, unit='deg')
DM_FRB=2641.5 
host_DM = 80.
MW_halo = 50. 


def grab_zlim():
    """
    Computes the minimum and maximum FRB redshift limits within a confidence 
    interval from DM
    
    Returns z_min, z_max
    
    """

    # NE 2001
    DM_ISM = mw.ismDM(frb_coords)

    # DM Cosmic: total - MW_ISM from ne2001 - MW_halo from user - host_DM from user
    DM_cosmic = DM_FRB - DM_ISM.value - MW_halo - host_DM 


    # Load z-DM distribution
    sdict = prob_dmz.grab_repo_grid()
    PDM_z = sdict['PDM_z']
    z = sdict['z']
    DM = sdict['DM']

    # find probability corresponding to calculated DM_cosmic
    iDM = np.argmin(np.abs(DM - DM_cosmic))
    PzDM = PDM_z[iDM, :] / np.sum(PDM_z[iDM, :])

    # create cumulative dist and set confidence limits
    cum_sum = np.cumsum(PzDM)
    limits = (2.5, 97.5)

    # find z corresponding to confidence limits in the CDF
    z_min = z[np.argmin(np.abs(cum_sum-limits[0]/100.))]
    z_max = z[np.argmin(np.abs(cum_sum-limits[1]/100.))]

    return z_min, z_max


def r_vs_dm(outfile='fig_r_vs_z.png',
               flipy=True, known_hosts = False):
    """
    Plots the intersection of galaxy apparent magnitude evolution with redshift 
    and redshift range of the FRB and saves it. 

    Args:
        outfile (str) : name of the output file. Default "fig_r_vs_z.png"
        flipy (bool) : flip order of y-axis (mag) values. Default True
        known_hosts (bool) : wish to overplot m_r - z for known FRB hosts. Default False
    
    """
    # get redshift limits given DM
    z_min, z_max = grab_zlim()

    # Function to convert redshift to m_r for Luminosity
    f_mL = frb_gal_u.load_f_mL()


    # set up the figure
    plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0])
    ax.grid(alpha=0.2)



    # if you wish to overplot known FRB hosts 
    if known_hosts == True:

        # Load up hosts
        host_tbl, _ = frb_gal_u.build_table_of_hosts()

        # Cut on PATH
        cut_path = (host_tbl.P_Ox > 0.8) & np.isfinite(host_tbl.P_Ox)
        # Cut on M_r
        cut_Mr = np.isfinite(host_tbl.M_r)

        good_hosts = host_tbl[cut_Mr & cut_path].copy()

        # Calculate apparaent magnitude
        m_r = []
        for kk in range(len(good_hosts)):
            ihost = good_hosts.iloc[kk]
            # Distance
            d = Distance(z=ihost.z)
            m_r.append(ihost.M_r + d.distmod.value)
        good_hosts['m_r'] = m_r 

        # plot the known hosts
        ax.plot(good_hosts.z, good_hosts.m_r, 'ok',
                label='Secure hosts')

    # plot
    # L curves
    zvals = np.linspace(0.021, 4.0, 200)
    m_Lstar = f_mL(zvals)
    ax.plot(zvals, m_Lstar, '-r', label='L*')

    m_01Lstar = m_Lstar + 2.5
    ax.plot(zvals, m_01Lstar, '--r', label='0.1 L*')

    m_001Lstar = m_Lstar + 5
    ax.plot(zvals, m_001Lstar, ':r', label='0.01 L*')


    # Add P(z|DM)
    xmnx = (0., 4.)
    ymnx = (14, 33.)
    if flipy:
        ymnx = (ymnx[1], ymnx[0])
    
    y_range = np.linspace(ymnx[0], ymnx[1], 100) 
    ax.fill_betweenx(y_range, x1=z_min, x2=z_max, 
                     color='lightgreen', alpha=0.5)
    #zlbl = 'P(z|DM) [95% c.l.]'
    zlbl = '95% c.l. FRB Redshift \n estimated from DM'
    ax.text(np.mean([z_min,z_max]), 20.5, zlbl,
            color='k', 
            size='large', ha='center')

    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$m_r$')
    ax.set_xlim(xmnx)
    ax.set_ylim(ymnx)

    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))

    ax.legend()

    frb_fig_u.set_fontsize(ax, 15.)

    # End
    plt.title(FRB)
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 300
    plt.savefig(outfile, **kwargs)
    plt.close()


''' Run me '''
# Command line execution
if __name__ == '__main__':
    r_vs_dm()