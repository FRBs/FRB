""" Estimate magnitude range for a host FRB given DM and plot it """

import numpy as np

#import healpy as hp

from cycler import cycler

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.coordinates import Distance
from linetools.scripts.utils import coord_arg_to_coord

from frb.frb import FRB, build_table_of_frbs, list_of_frbs
from frb.figures import utils as frb_fig_u
from frb.galaxies import utils as frb_gal_u

from IPython import embed


def r_vs_dm_figure(z_min, z_max, z, PzDM, outfile='fig_r_vs_z.png',
               flipy=True, known_hosts = False, title=None):
    """
    Plots the intersection of galaxy apparent magnitude evolution with redshift 
    and redshift distribution of the FRB and saves it. 

    Args:
        z_min (float) : min redshift of the confidence interval
        z_max (float) : max redshift of the confidence intervale
        z (array) : redshift array
        PzDM (array) : pdf(z) PDF of redshifts for a given DM
        outfile (str) : name of the output file. Default "fig_r_vs_z.png"
        flipy (bool) : flip order of y-axis (mag) values. Default True
        known_hosts (bool) : wish to overplot m_r - z for known FRB hosts. Default False
    
    """

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


    zvals = np.linspace(0.021, z_max+0.5, 200)
    m_Lstar = f_mL(zvals)
    ax.plot(zvals, m_Lstar, '-r', label='L*')

    m_01Lstar = m_Lstar + 2.5
    ax.plot(zvals, m_01Lstar, '--r', label='0.1 L*')

    m_001Lstar = m_Lstar + 5
    ax.plot(zvals, m_001Lstar, ':r', label='0.01 L*')


    # Add P(z|DM)
    xmnx = (0., z_max+0.5)
    ymnx = (14, max(15, m_001Lstar[-1]+0.5))
    if flipy:
        ymnx = (ymnx[1], ymnx[0])
    
    y_range = np.linspace(ymnx[0], ymnx[1], 300) 
    x_range = np.linspace(z_min,z_max,300)

    pzd = np.interp(x_range,z,PzDM)
    X,Y = np.meshgrid(x_range,y_range)
    Z = X/X*pzd
    #ax.fill_betweenx(y_range, x1=z_min, x2=z_max, 
    #                 color='lightgreen', alpha=0.5)
    #zlbl = 'P(z|DM) [95% c.l.]'
    c=ax.pcolor(X, Y, Z*1000, cmap='Blues')
    zlbl = '95% c.l. FRB Redshift \n estimated from DM'
    text_x = 0.35 * (xmnx[1] - xmnx[0]) + xmnx[0]
    text_y = 0.5 * (ymnx[1] - ymnx[0]) + ymnx[0]
    ax.text(text_x, text_y, zlbl,
            color='k', 
            size='large', ha='center')

    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$m_r$')
    ax.set_xlim(xmnx)
    ax.set_ylim(ymnx)

    if z_max <= 0.5:
        ax.set_xscale('log')
        ax.set_xlim(1e-2,z_max+0.5)
        ax.xaxis.set_major_locator(plt.LogLocator(base=10, numticks=12))
    else:
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    plt.colorbar(c,label='p(z|DM) [a.u.]',ax=ax)

    ax.legend(loc='upper right')
    # set the title of the figure
    ax.set_title(title)
    frb_fig_u.set_fontsize(ax, 15.)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 300
    plt.savefig(outfile, **kwargs)
    plt.close()

