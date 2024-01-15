from typing import IO
import numpy as np

import glob, os, sys, json
import pdb

#import healpy as hp

from cycler import cycler
import matplotlib as mpl
from numpy.core.fromnumeric import mean
import seaborn as sns

import pandas
import aplpy

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import ConnectionPatch

import matplotlib
from matplotlib import rc
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap

from scipy.interpolate import interp1d
from scipy.stats import kstest
from scipy.stats import norm, poisson
from scipy.signal import convolve

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time

from IPython import embed

# Local
sys.path.append(os.path.abspath("../../Analysis/py"))
import alopeke_defs
import alopeke_utils2
import alopeke_analy

# Globals
handletextpad=0.3

def ccolor(r,g,b):
    return((r/255.,g/255.,b/255., 1.0))


def fig_fov(outroot='fov_gmos_alopeke'):
    # images
    #gmos_image = '/media/ntejos/caladan01/projects/FRBs/FRB180916-R3/GMOS/FRB180916_Z_1320s_GaiaDR2ac.fits'
    #alopeke_image = '/media/ntejos/caladan01/data/FRBs/R3/Alopeke/science/N20201023A0015r.fits'  # just as an example
    gmos_image = '/Users/ckilpatrick/Dropbox/Data/FRB/FRB180916/Data/Gemini/frb180916.r.ut190713_0001.fits'
    alopeke_image = '/Users/ckilpatrick/Dropbox/Data/FRB/FRB180916/Data/Alopeke/20220908/N20220908A0017r_reduced.fits'

    # define positions of objects
    # We use the position from Marcote+ paper because the one in the repo is slightly different
    frb_coord = alopeke_defs.frb180916.coord  # in Gaia DR2 frame (as given by Marcote+20; ok in our Repo)
    
    star1_coord = SkyCoord(29.4954470, 65.7189804, unit='deg')  # From GaiaDR2 frame (to match that of the FRB)
    star2_coord = SkyCoord(29.4908603, 65.7128539, unit='deg')

    star1_xy_red = (61.515341,197.51952)

    bkg_coord = SkyCoord(frb_coord.ra.value,star1_coord.dec.value, unit='deg')  # place holder check with consuelo
    
    # need to add WCS to Alopeke image (blue)
    hdul_alopeke = fits.open(alopeke_image)

    #update header properly
    alopeke_header = hdul_alopeke[0].header
    alopeke_header['CDELT1']=-np.sqrt(alopeke_header['PC1_1']**2+alopeke_header['PC1_2']**2)
    alopeke_header['CDELT2']=-np.sqrt(alopeke_header['PC2_1']**2+alopeke_header['PC2_2']**2)
    alopeke_header['CD1_1'] = alopeke_header['PC1_1']
    alopeke_header['CD2_2'] = alopeke_header['PC2_2']
    alopeke_header['CD2_1'] = alopeke_header['PC2_1']
    alopeke_header['CD1_2'] = alopeke_header['PC1_2']
    alopeke_header['CRPIX1'] = star1_xy_red[0]
    alopeke_header['CRPIX2'] = star1_xy_red[1]
    alopeke_header['CRVAL1'] = star1_coord.dec.value  # x is declination for red camera
    alopeke_header['CRVAL2'] = star1_coord.ra.value  # y is RA for red camera

    # Write new HDU as image, getitng only the frame 45 as reference
    alopeke_header['NAXIS'] = 2
    alopeke_header.remove('NAXIS3')
    hdul_alopeke[0].data = hdul_alopeke[0].data[45]
    hdul_alopeke[:-1].writeto('data/alopeke_aux.fits', overwrite=True)

    wt = WCS(naxis=2)
    wt.wcs.ctype = ['DEC--TAN', 'RA---TAN']
    wt.wcs.crval = alopeke_header['CRVAL1'], alopeke_header['CRVAL2']
    wt.wcs.crpix = alopeke_header['CRPIX1'], alopeke_header['CRPIX2']
    wt.wcs.cdelt = alopeke_header['CDELT1'], alopeke_header['CDELT2']
    header = wt.to_header()

    hdu = fits.PrimaryHDU(hdul_alopeke[0].data, header=header)
    hdul = fits.HDUList([hdu])
    hdul.writeto('data/alopeke_aux.fits', overwrite=True)

    # new string pointing to this new file, used for Fig.2 of paper
    alopeke_image_example = './data/alopeke_aux.fits'

    # Global things
    lw = 1  # lw of square and connector
    color= 'k' # color square and connector
    ls = 'dashed'
    c_frb = 'r'
    ls_frb='solid'
    radius_frb = 1.2/3600
    c_ref = 'orange'
    ls_ref='solid'
    radius_ref = 1.2/3600
    c2_ref = 'magenta'
    lw_c = 1.5 # lw of circles
    ls_bkg = 'solid'
    c_bkg = 'y'
    radius_bkg = 1.2/3600

    # cmap = 'afmhot_r'
    cmap = 'Blues'

    # relative geometry
    # center_gmos = (29.5017349,65.7145370)
    # center_alopeke = (29.5024181,65.7164417)
    center_alopeke = wt.pixel_to_world(128,128)
    center_alopeke = center_alopeke.ra.value, center_alopeke.dec.value
    center_gmos = center_alopeke[0],center_alopeke[1] + 2/3600.
    size_alopeke = 37/3600.  # 37"
    size_gmos = 60/3600. # 1 arcmin
    vmin_gmos = 0.
    vmax_gmos = 600.
    vmin_alopeke = -50.
    vmax_alopeke= 500.

    # Test for cuts
    if 0:
        vmin,vmax = 500,1000
        f3 = aplpy.FITSFigure(alopeke_image_example)
        f3.show_colorscale(vmin=vmin,vmax=vmax, cmap=cmap)
        # f3.recenter(center_gmos[0],center_gmos[1], radius=size_gmos/2.)
        f3.add_colorbar()
        plt.show()
        stop

    # MAIN FIG
    fig = plt.figure(figsize=(7,4))
    fig.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.95,wspace=0.25,hspace=0.25)

    f1 = aplpy.FITSFigure(gmos_image, figure=fig, subplot=(1,2,1), north=True)
    f2 = aplpy.FITSFigure(alopeke_image_example, figure=fig, subplot=(1,2,2), north=True)

    f1.show_colorscale(vmin=vmin_gmos,vmax=vmax_gmos, cmap=cmap)
    f1.recenter(center_gmos[0],center_gmos[1], radius=size_gmos/2.)
    f2.recenter(center_alopeke[0], center_alopeke[1], radius=size_alopeke/2.)
    f2.show_colorscale(vmin=vmin_alopeke,vmax=vmax_alopeke, cmap=cmap)

    #add rectangle to the fig1 with alopeke Zoom
    xw = center_alopeke[0]
    yw = center_alopeke[1]
    f1.show_rectangles(xw, yw, size_alopeke, size_alopeke, zorder=1, color=color, lw=lw, ls=ls)

    #connectors AX1 to AX2
    axA = f1.ax
    axB = f2.ax
    xyA1 = f1.world2pixel(xw-size_alopeke/2./np.cos(yw*u.deg), yw+size_alopeke/2.)
    xyB1 = (0,1)
    xyA2 = f1.world2pixel(xw-size_alopeke/2./np.cos(yw*u.deg), yw-size_alopeke/2.)
    xyB2 = (0,0)        
    con1 = ConnectionPatch(xyA=xyA1, xyB=xyB1, coordsA="data", coordsB="axes fraction", \
                        axesA=axA, axesB=axB, color=color, lw=lw, ls=ls)
    con2 = ConnectionPatch(xyA=xyA2, xyB=xyB2, coordsA="data", coordsB="axes fraction", \
                        axesA=axA, axesB=axB, color=color, lw=lw, ls=ls)
    axA.add_artist(con1)
    axA.add_artist(con2)

    # f1.add_scalebar((10*u.arcsec).to('deg').value)
    f2.add_scalebar((10*u.arcsec).to('deg').value)
    f2.scalebar.set_corner('bottom left')
    f2.scalebar.set_label('10"')
    f2.scalebar.set_color('k')
    f2.scalebar.set_font_size('large')

    # add FRB position and reference star
    frb_ra = frb_coord.ra.value
    frb_dec = frb_coord.dec.value
    star_ra = star1_coord.ra.value
    star_dec = star1_coord.dec.value
    star2_ra = star2_coord.ra.value
    star2_dec = star2_coord.dec.value
    bkg_ra = bkg_coord.ra.value 
    bkg_dec = bkg_coord.dec.value

    f2.show_circles(frb_ra, frb_dec, radius_frb, zorder=1, color=c_frb, lw=lw_c, ls=ls_frb)
    f1.show_circles(frb_ra, frb_dec, radius_frb, zorder=1, color=c_frb, lw=lw_c, ls=ls_frb)
    f2.show_circles(star_ra, star_dec, radius_ref, zorder=1, color=c_ref, lw=lw_c, ls=ls_ref)
    f1.show_circles(star_ra, star_dec, radius_ref, zorder=1, color=c_ref, lw=lw_c, ls=ls_ref)
    f2.show_circles(star2_ra+0.0005, star2_dec, radius_ref, zorder=1, color=c2_ref, lw=lw_c, ls=ls_ref)
    f1.show_circles(star2_ra, star2_dec, radius_ref, zorder=1, color=c2_ref, lw=lw_c, ls=ls_ref)

    f1.add_label(frb_ra, frb_dec+3*radius_frb, 'FRB pos.', color=c_frb)
    f1.add_label(star_ra+0.5*radius_ref, star_dec-3.5*radius_ref, 'Star-1', color=c_ref)
    f1.add_label(star2_ra+0.5*radius_ref, star2_dec+3.5*radius_ref, 'Star-2', color=c2_ref)

    # ADD labels
    f1.add_label(0.53,0.97,'GMOS r-band', relative='axes', va='top')
    f2.add_label(0.5,0.97,'`Alopeke i-band\n(~10 ms)\n 2020-10-23', relative='axes', va='top')

    #Hide tick labels for f2
    f2.tick_labels.hide()
    f2.axis_labels.hide()

    # write
    fig.savefig(outroot+'.pdf', format='pdf')


def fig_chk_gauss(date, outroot='fig_chk_gauss_', camera='red',
                  show_poisson=False):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()
    
    # Load data
    if camera == 'red':
        clr = 'r'
    elif camera == 'blue':
        clr = 'b'
    else:
        raise IOError(f'Bad camera: {camera}')
    
    outfile = f'{outroot}{camera}_{date}.pdf'
        
    # Load up
    data_dict = alopeke_utils2.load_camera(camera, date, cut=True)
    C_gal = data_dict['C_gal']
    mean_gal, std_gal = data_dict['Gauss_gal'] # Gaussian
    print("mean={}, std={}".format(mean_gal, std_gal))

    poiss_gal = np.mean(C_gal)

    # Plot
    psz = 13.
    plt.figure(figsize=(6, 9))
    gs = gridspec.GridSpec(3,1)

    # Main Gaussian
    ax0 = plt.subplot(gs[0])
    ax0.hist(C_gal, bins=50, density=True, color=clr)
    xmin, xmax = C_gal.min(), C_gal.max()
    x = np.linspace(xmin, xmax, 1000)
    y = norm.pdf(x, mean_gal, std_gal)
    ax0.plot(x, y, 'k')
    if show_poisson:
        x2 = np.arange(int(xmin), int(xmax)+1)
        y2 = poisson.pmf(x2, poiss_gal)
        ax0.plot(x2, y2, 'k--')

    # Label me
    xlbl = r'$C_{\rm FRB}^{\rm Tot} \; (\rm e^-/exposure)$'
    ax0.set_xlabel(xlbl)
    ax0.set_ylabel('PDF')
    set_fontsize(ax0, psz)

    # High-tail
    isort_gal = np.argsort(C_gal)
    cdf_gal = (np.arange(isort_gal.size)+1) / isort_gal.size

    ax1 = plt.subplot(gs[2])
    xval = np.linspace(mean_gal+std_gal, C_gal.max(), 100000)
    ax1.step(C_gal[isort_gal], 1.-cdf_gal, color=clr, where='mid')
    ax1.plot(xval, 1-norm.cdf(xval, loc=mean_gal, scale=std_gal),
             color='k')
    #
    ax1.set_xlim(mean_gal+std_gal, np.max(C_gal))
    ax1.set_ylim(1e-6, 1e-1)
    #
    ax1.set_yscale('log')
    ax1.set_xlabel(xlbl)
    ax1.set_ylabel('1-CDF')
    set_fontsize(ax1, psz)

    # Low tail
    ax2 = plt.subplot(gs[1])
    xval2 = np.linspace(C_gal.min(), mean_gal, 100000)
    
    # Data
    ax2.step(C_gal[isort_gal], cdf_gal, color=clr, where='mid')
    ax2.plot(xval2, norm.cdf(xval2, loc=mean_gal, scale=std_gal),
             color='k')
    #
    ax2.set_xlim(C_gal.min(), mean_gal-std_gal)
    ax2.set_ylim(1e-6, 1e-1)

    ax2.set_yscale('log')
    ax2.set_xlabel(xlbl)
    ax2.set_ylabel('CDF')
    set_fontsize(ax2, psz)
    
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    elif 'pdf' in outfile:
        kwargs['format'] = 'pdf'
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_upper_limit(date, outroot='fig_upper_limit_', camera='red',
                    step=200, nsamp=100):
    if camera == 'red':
        clr = 'r'
        ylim = (1e-4, 1e-1)
    elif camera == 'blue':
        clr = 'b'
        ylim = (1e-4, 3e-1)
    elif camera == 'both':
        ylim = (1e-4, 1e-1)
    else:
        raise IOError(f'Bad camera: {camera}')

    outfile = f'{outroot}{camera}_{date}.pdf'

    # Analyze
    if camera == 'both':
        C_FRB_val_blue, p_exceed_blue, upper_FRB_blue = alopeke_analy.upper_limit(date, camera='blue')
        C_FRB_val_red, p_exceed_red, upper_FRB_red = alopeke_analy.upper_limit(date, camera='red')
    else:
        C_FRB_val, p_exceed, upper_FRB = alopeke_analy.upper_limit(date, camera)

    
    date_str=alopeke_defs.FRB_time[date].datetime.strftime('%Y-%m-%d')

    # Figure time
    psz = 13.
    plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(1,1)

    ax = plt.subplot(gs[0])

    if camera == 'both':
        ax.plot(C_FRB_val_blue, 1-p_exceed_blue, color='blue')
        ax.plot(C_FRB_val_red, 1-p_exceed_red, color='red')
    else:
        ax.plot(C_FRB_val, 1-p_exceed, color=clr)

    ax.axhline(0.0027, color='k', ls='--')

    ax.set_xlabel(r'$\mu_{\rm FRB} \; (\rm e^-)$')
    ax.set_ylabel('Fraction of Events')

    ax.set_ylim(ylim)
    ax.set_yscale('log')

    set_fontsize(ax, psz)

    ax.text(0.8, 0.9, date_str, color='k', transform=ax.transAxes)
    
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    elif 'pdf' in outfile:
        kwargs['format'] = 'pdf'
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_gal_star_corr(date, outroot='fig_gal_star_corr_', camera='red',
                      cut=True):
    if camera == 'red':
        clr = 'r'
    elif camera == 'blue':
        clr = 'b'
    else:
        raise IOError(f'Bad camera: {camera}')

    outfile = f'{outroot}{camera}_{date}.pdf'

    # Load up
    data_dict = alopeke_utils2.load_camera(camera, date, cut=cut)
    C_gal = data_dict['C_gal']
    #C_star = data_dict['C_star1']
    C_star = data_dict['C_star2']

    # Fit Gaussian
    mean_gal, std_gal = norm.fit(C_gal)
    mean_star, std_star = norm.fit(C_star)

    # Convert to sig
    nsig_gal = (C_gal-mean_gal)/std_gal
    nsig_star = (C_star-mean_star)/std_star

    df = pandas.DataFrame(dict(nsig_gal=nsig_gal,
                               nsig_star=nsig_star))
    
    # Figure time
    psz = 13.
    plt.figure(figsize=(6, 5))
    
    fg = sns.displot(df, x='nsig_gal', y='nsig_star',
                    color=clr)
    fg.ax.set_xlabel(r'$\Delta C_{\rm FRB}^{\rm Tot} \; (\sigma_{\rm FRB})$')
    fg.ax.set_ylabel(r'$\Delta C_{\rm 1}^{\rm Tot} \; (\sigma_{\rm 1})$')

    fg.ax.set_aspect('equal')

    fg.ax.set_xlim([-6,6])
    fg.ax.set_ylim([-6,6])

    set_fontsize(fg.ax, psz)
    
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    elif 'pdf' in outfile:
        kwargs['format'] = 'pdf'
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_frb_counts(date, outfile='fig_frb_counts', camera='red',
                    step=200, nsamp=100):
    set_mplrc()


    # Figure time
    psz = 13.
    plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(2,1)

    for ss, clr, camera in zip(range(2), ['r', 'b'], ['red', 'blue']):
        # Load up
        data_dict = alopeke_utils2.load_camera(camera, date, cut=True)

        ax = plt.subplot(gs[ss])

        date_str=alopeke_defs.FRB_time[date].datetime.strftime('%Y-%m-%d')

        # All points
        ax.scatter((data_dict['MJD_gal']-alopeke_defs.FRB_time[date].mjd)*24*3600,
                data_dict['C_gal'], color='k')
        # In time window
        ax.scatter((data_dict['MJD_FRB']-alopeke_defs.FRB_time[date].mjd)*24*3600,
                data_dict['C_FRB'], color=clr)

        # Limits
        ax.set_xlim(-1., 1.)
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
        ax.minorticks_on()

        print(data_dict['MJD_FRB']-alopeke_defs.FRB_time[date].mjd)

        mask = np.abs(data_dict['MJD_FRB']-alopeke_defs.FRB_time[date].mjd)<163.0/2/86400/1000
        print(f'{camera} data:')
        for idx in np.where(mask)[0]:
            mjd = data_dict['MJD_FRB'][idx]
            mjd = '%5.8f'%mjd
            filt = 'r'
            if camera=='r' or camera=='red':
                filt = 'i'
            bkg = '%.2f'%data_dict['C_bkg'][idx]
            star1 = '%.2f'%data_dict['C_star1'][idx]
            star2 = '%.2f'%data_dict['C_star2'][idx]
            frb = '%.2f'%data_dict['C_FRB'][idx]

            # Output measurements for latex table
            print(f'{mjd} & {bkg} & {star1} & {star2} & {frb} \\\\')

        ax.text(0.02, 0.9, date_str, color='k', transform=ax.transAxes)

        ax.set_ylabel(r'$C_{\rm FRB}^{\rm TOT} \, ({\rm e^-})$')

        if ss == 0:
            ax.axes.xaxis.set_ticklabels([])
        else:
            ax.set_xlabel(r'$t - t_{\rm FRB} \;$ (seconds)')

        #ax.set_ylim(1e-4, 1)
        #ax.set_yscale('log')

        set_fontsize(ax, psz)
        # Stats
        print(f"Camera: {camera}")
        print("Nobs = {}".format(len(data_dict['C_FRB'])))
        print("Max = {}".format(max(data_dict['C_FRB'])))
        
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    outfile=outfile+'_'+date+'.pdf'
    if 'png' in outfile:
        kwargs['dpi'] = 700
    elif 'pdf' in outfile:
        kwargs['format'] = 'pdf'
    plt.savefig(outfile, **kwargs)
    plt.close()



def fig_all_counts(date, outfile='fig_all_counts.pdf', zoom=False):
    set_mplrc()

    # Load up
    data_red = alopeke_utils2.load_camera('red', date, cut=True)
    data_blue = alopeke_utils2.load_camera('blue', date, cut=True)

    mjd_day = int(data_red['MJD_gal'][0])

    mjd_FRB = alopeke_defs.FRB_time[date].mjd
    print("MJD FRB = {}".format(mjd_FRB))

    # Figure time
    psz = 13.
    plt.figure(figsize=(12, 5))
    gs = gridspec.GridSpec(2,1)


    for ss, clr, data in zip([0,1], ['b', 'r'], [data_blue, data_red]):
        ax = plt.subplot(gs[ss])
        # FRB
        ax.axvline(mjd_FRB-mjd_day, ls='--', color='gray', zorder=4)
        ax.axvspan(alopeke_defs.mjd_low[date] - mjd_day, alopeke_defs.mjd_high[date] - mjd_day, 0,100, color='gray', alpha=0.5, zorder=4)

        # Data
        ax.scatter(data['MJD_gal']-mjd_day, data['C_gal'], color=clr, s=1, zorder=1)
        ax.scatter(data['MJD_FRB']-mjd_day, data['C_FRB'], color=clr, s=1, zorder=1)
        if ss == 0:
            ax.get_xaxis().set_ticks([])
        set_fontsize(ax, psz)
        ax.set_ylabel(r'$C_{\rm FRB}^{\rm Tot} \, ({\rm e^-/exposure})$')

        # Limits
        if zoom is True:
            xmin = (alopeke_defs.FRB_time - 7*alopeke_defs.ata).mjd - mjd_day
            xmax = (alopeke_defs.FRB_time + 7*alopeke_defs.ata).mjd - mjd_day
            ax.set_xlim(xmin,xmax)

    ax.set_xlabel('MJD - {}'.format(mjd_day))

    
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    if zoom is True:
        outfile = outfile.replace('.','_zoom.')
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 500
    elif 'pdf' in outfile:
        kwargs['format'] = 'pdf'
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_star_gal_counts(date, outfileroot='fig_star_gal_counts.pdf', cam='red',
    star_factors=[6], use_stars=['star1'], ymax=200):
    set_mplrc()

    # Load up
    data = alopeke_utils2.load_camera(cam, date, cut=True)

    mjd_day = int(data['MJD_FRB'][0])

    mjd_FRB = alopeke_defs.FRB_time[date].mjd
    print("MJD FRB = {}".format(mjd_FRB))

    # Figure time
    psz = 13.
    plt.figure(figsize=(9, 3))
    gs = gridspec.GridSpec(1,1)



    #for ss, clr, data in zip([0,1], ['b', 'r'], [data_blue, data_red]):
    if 1:
        clr = cam
        ax = plt.subplot(gs[0])
        # FRB
        ax.axvline(mjd_FRB-mjd_day, ls='--', color='gray', zorder=4)
        ax.axvspan(alopeke_defs.mjd_low[date] - mjd_day, alopeke_defs.mjd_high[date] - mjd_day, 0,100, color='gray', alpha=0.5, zorder=4)

        # Data
        ax.scatter(data['MJD_gal']-mjd_day, data['C_gal'], color=clr, s=1, zorder=1, label='FRB position')
        ax.scatter(data['MJD_FRB']-mjd_day, data['C_FRB'], color=clr, s=1, zorder=1)
        colors = ['orange','purple']
        i=0
        for star,fact in zip(use_stars,star_factors):
            nstar=star.replace('star','')
            ax.scatter(data[f'MJD_{star}_full']-mjd_day,
                data[f'C_{star}_full']/fact, color=colors[i],
                alpha=0.8, s=1, zorder=1,
                label=f'Ref. star {nstar} (scaled by 1/{fact})')
            i=i+1
        
        # label camera and filter
        frb_date = Time(mjd_FRB, format='mjd')
        date_str = frb_date.datetime.strftime('%Y-%m-%d')
        if cam == 'red':
            label = 'Red camera (i-band) \n'+date_str
        elif cam == 'blue':
            label = 'Blue camera (r-band) \n'+date_str
        ax.text(0.05, 0.85, label, color='k', transform=ax.transAxes)

        set_fontsize(ax, psz)
        ax.set_ylabel(r'$ {\rm Total \ counts} \ ({\rm e^-/exposure})$')
        

    ax.set_xlabel('MJD - {}'.format(mjd_day))

    ax.set_ylim(0, ymax)

    # minorticks
    ax.minorticks_on()
    
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    outfile = outfileroot.replace('.pdf', '_{}.pdf'.format(cam+'_'+date))
    print('Writing {:s}'.format(outfile))
    plt.legend(loc='upper right')
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 500
    elif 'pdf' in outfile:
        kwargs['format'] = 'pdf'
    plt.savefig(outfile, **kwargs)

    plt.close()


def fig_model_afterglow_limits(date, outfileroot='fig_model_afterglow_limits.pdf'):

    set_mplrc()

    # Load up
    Evals, nvals, grid = alopeke_analy.afterglow_limits(date)

    # Figure time
    psz = 13.
    plt.figure(figsize=(3.5, 3.5))
    gs = gridspec.GridSpec(1,1)

    ax = plt.subplot(gs[0])

    ncols = 60
    colarray=np.array([ccolor(l,l,l) for l in np.flip(np.arange(255-ncols,255))])
    cm = ListedColormap(colarray)

    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    levels = np.linspace(0, 1, ncols)

    Erange = 10**Evals*1e43
    nrange = 10**nvals

    cf=ax.contourf(Erange, nrange, grid, cmap=cm, norm=norm, levels=levels)
    C=ax.contour(Erange, nrange, grid, linewidths=2, levels=levels,
        zorder=5, colors=('red'))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([np.min(Erange), np.max(Erange)])
    ax.set_ylim([np.min(nrange), np.max(nrange)])

    set_fontsize(ax, psz)

    ax.set_xlabel('Burst Energy (erg)')
    ax.set_ylabel(r'Circumburst Density (cm$^{-3}$)')

    # minorticks
    ax.minorticks_on()

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    outfile = outfileroot.replace('.pdf', '_{}.pdf'.format(date))
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 500
    elif 'pdf' in outfile:
        kwargs['format'] = 'pdf'
    plt.savefig(outfile, **kwargs)

    plt.close()


def set_mplrc():
    mpl.rcParams['mathtext.default'] = 'it'
    mpl.rcParams['font.size'] = 12
    mpl.rc('font',family='Times New Roman')
    mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
    mpl.rc('text', usetex=True)


def set_fontsize(ax,fsz):
    '''
    Parameters
    ----------
    ax : Matplotlib ax class
    fsz : float
      Font size
    '''
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)

def set_mplrc():
    mpl.rcParams['mathtext.default'] = 'it'
    mpl.rcParams['font.size'] = 12
    #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #mpl.rc('font',family='Times New Roman')
    #mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
    mpl.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
    mpl.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
    mpl.rc('text', usetex=True)


#### ########################## #########################
def main(flg_fig, date, cameras=['blue','red']):

    if flg_fig == 'all':
        flg_fig = np.sum( np.array( [2**ii for ii in range(25)] ))
    else:
        flg_fig = int(flg_fig)

    # Check Gaussian
    if flg_fig & (2**0):
        for camera in cameras:
            fig_chk_gauss(date, camera=camera)

    # Upper limit
    if flg_fig & (2**1):
        if 'blue' in cameras and 'red' in cameras:
            fig_upper_limit(date, camera='both')
        else:
            for camera in cameras:
                fig_upper_limit(date, camera=camera)

    # Correlation?
    if flg_fig & (2**2):
        for camera in cameras:
            fig_gal_star_corr(date, camera=camera)

    # FRB counts
    if flg_fig & (2**3):
        if 'blue' in cameras and 'red' in cameras:
            fig_frb_counts(date)

    # All counts: FRB counts for both cameras
    if flg_fig & (2**4):
        if 'blue' in cameras and 'red' in cameras:
            fig_all_counts(date)
        
    # counts for FRB + star in the same panel for each camera
    if flg_fig & (2**5):
        for camera in cameras:
            if camera == 'blue':
                ymax = 120
                use_stars = ['star1','star2']
                star_factors = [34, 20]
            elif camera == 'red':
                ymax = 120
                use_stars = ['star1','star2']
                star_factors = [45, 30]
            fig_star_gal_counts(date, cam=camera, ymax=ymax, star_factors=star_factors,
                use_stars=use_stars)

    # FOV
    if flg_fig & (2**6):
        fig_fov()

    if flg_fig & (2**7):
        fig_model_afterglow_limits(date)
    


# Command line execution
if __name__ == '__main__':

    date = sys.argv[1]

    if len(sys.argv) == 2:
        flg_fig = 0
        flg_fig += 2**0   # Check Gaussianity
        flg_fig += 2**1   # Upper limit
        flg_fig += 2**2   # Star/gal correlation
        flg_fig += 2**3   # FRB counts
        flg_fig += 2**4   # All counts (both cameras)
        flg_fig += 2**5   # FRB counts + Star counts (each camera)
        flg_fig += 2**6   # FOV
    else:
        flg_fig = sys.argv[2]

    main(flg_fig, date, cameras=['blue','red'])
