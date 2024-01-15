from astropy.io import fits
from astropy.table import Table, vstack
import sep
import numpy as np
import UTILS as ut
import os
import sys
import glob
import copy

from astropy.time import Time
from astropy import units as u
from astropy.wcs import WCS, utils
from astropy.coordinates import SkyCoord

from astroquery.mast import Catalogs

from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats

import progressbar

def correct_masked_median_columns(fitsfile, camera):

    image_data = np.median(fitsfile[0].data, axis=0)

    for i in np.arange(fitsfile[0].header['NAXIS2']):
        med = np.median(image_data[:,i])
        fitsfile[0].data[:,:,i] = fitsfile[0].data[:,:,i]-med

    return(fitsfile)


def get_frame_alignment(fitsfile, coord, camera, catalog, save_stack=''):

    image_data = np.median(fitsfile[0].data, axis=0)

    newhdu = fits.PrimaryHDU()
    newhdu.data = image_data

    # Initial WCS guess
    newhdu.header['CTYPE1']='RA---TAN'
    newhdu.header['CTYPE2']='DEC--TAN'
    newhdu.header['CRVAL1']=coord.ra.degree
    newhdu.header['CRVAL2']=coord.dec.degree
    if camera=='r':
        newhdu.header['CRPIX1']=113.92306
        newhdu.header['CRPIX2']=117.73306
        newhdu.header['CD1_1']=0.
        newhdu.header['CD1_2']=-0.00004138888
        newhdu.header['CD2_1']=-0.00004138888
        newhdu.header['CD2_2']=0.
    if camera=='b':
        newhdu.header['CRPIX1']=128.0
        newhdu.header['CRPIX2']=128.0
        newhdu.header['CD1_1']=0.
        newhdu.header['CD1_2']=-0.00004138888
        newhdu.header['CD2_1']=0.00004138888
        newhdu.header['CD2_2']=0.

    wcs = WCS(newhdu.header)

    # Get xy coordinates of all sources in catalog
    coords = SkyCoord(catalog['ra'], catalog['dec'], unit='deg')
    xs, ys = utils.skycoord_to_pixel(coords, wcs, origin=1)

    mask = (xs > 0) & (xs < fitsfile[0].header['NAXIS1']) &\
        (ys > 0) & (ys < fitsfile[0].header['NAXIS2'])

    catalog = catalog[mask]

    n=len(catalog)
    print(f'Got {n} sources for alignment')

    mean, median, std = sigma_clipped_stats(image_data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
    sources = daofind(image_data - median)

    # Do matching from catalog->sources
    catalog.sort('phot_g_mean_mag')

    temp_sources = copy.copy(sources)
    matched_sources = []
    for row in catalog:
        if len(temp_sources)==0: break
        c = SkyCoord(row['ra'], row['dec'], unit='deg')
        x,y=utils.skycoord_to_pixel(c, wcs, origin=1)

        sep = (temp_sources['xcentroid']-x)**2+(temp_sources['ycentroid']-y)**2
        idx = np.argmin(sep.data)

        if sep[idx]>400: continue

        matched_sources.append((row['ra'], row['dec'],
            temp_sources[idx]['xcentroid'], temp_sources[idx]['ycentroid']))

        temp_sources.remove_row(idx)

    n=len(matched_sources)
    print(f'Got {n} matched sources for alignment')

    coords = SkyCoord([m[0] for m in matched_sources],
        [m[1] for m in matched_sources], unit='deg')
    c = coords[0]

    xy = np.array([[m[2] for m in matched_sources],[m[3]
        for m in matched_sources]])

    wcs = utils.fit_wcs_from_points(xy=xy, world_coords=coords, proj_point=c,
        projection='TAN')
    wcs_head = wcs.to_header()
    for key in wcs_head.keys():
        newhdu.header[key] = wcs_head[key]
        fitsfile[0].header[key] = wcs_head[key]

    if save_stack:
        newhdu.writeto(save_stack, overwrite=True, output_verify='silentfix')

    return(fitsfile)

tipo = "reduced"     #change filename

data_path = sys.argv[1]
date = sys.argv[2]
camera = sys.argv[3]

star1_ra = 29.4956219
star1_dec = 65.7191098
star1 = SkyCoord(star1_ra, star1_dec, unit='deg')
star2_ra = 29.4916103
star2_dec = 65.7187927
star2 = SkyCoord(star2_ra, star2_dec, unit='deg')
star3_ra = 29.4911121
star3_dec = 65.7129930
star3 = SkyCoord(star3_ra, star3_dec, unit='deg')

frb_ra = 29.5031258
frb_dec = 65.7167542
frb = SkyCoord(frb_ra, frb_dec, unit='deg')

bkg_ra = 29.5076445
bkg_dec = 65.7183129
bkg_coord = SkyCoord(bkg_ra, bkg_dec, unit='deg')

radio = 2
FWHM = 2*radio
sizebox = 20

master = None

for filename in sorted(glob.glob(os.path.join(data_path, date,
    f'N{date}A*{camera}_{tipo}.fits'))):

    basefile = os.path.basename(filename)
    n = basefile.replace(f'{camera}_{tipo}.fits','')
    n = n.replace(f'N{date}A','')

    origfile = basefile.replace(f'_{tipo}','')
    origfile = os.path.join(data_path, date, 'raw', origfile)

    if not os.path.exists(filename):
        raise Exception(f'{filename} does not exist.  Stop!')
    if not os.path.exists(origfile):
        raise Exception(f'{origfile} does not exist.  Stop!')

    fitsfile = fits.open(filename)
    print(f'Correcting bias rows in {filename}...')
    fitsfile = correct_masked_median_columns(fitsfile, camera)
    print('Done.')
    orighdu = fits.open(origfile)

    print(f'Getting WCS alignment for stacked frame {filename}...')
    coord = SkyCoord(orighdu[0].header['RA'],
        orighdu[0].header['DEC'], unit=(u.hour, u.deg))
    catalog_data = Catalogs.query_region(coord, radius=0.1,
            catalog="Gaia", version=2)
    mask = catalog_data['phot_g_mean_mag'] < 21.0

    basepath = os.path.dirname(filename)
    stackpath = os.path.join(basepath, 'stack')
    if not os.path.exists(stackpath):
        os.makedirs(stackpath)
    save_stack = os.path.join(stackpath,basefile.replace('.fits','.stack.fits'))
    fitsfile = get_frame_alignment(fitsfile, coord, camera, catalog_data,
        save_stack=save_stack)

    fitsfile.writeto(filename, overwrite=True, output_verify='silentfix')

    wcs_head = WCS(fitsfile[0].header).to_header()

    # Start time for all frames
    start_utc = Time(orighdu[0].header['OBSTIME'], format='unix')
    end_utc = Time(orighdu[0].header['EXPENDTM'], format='unix')

    # In principle, the total time per frame should be (start_utc-end_utc)/num_frames
    time_per_frame = (start_utc-end_utc)/orighdu[0].header['NAXIS3']

    # We also want the exposure time, which is slightly different accounting for dead time
    time_per_exposure = orighdu[0].header['EXPTIME']

    data = fitsfile[0].data
    data = data.astype(float)

    ind = []
    flux_1FWHM = []
    flux_2FWHM = []
    fluxerr_1FWHM = []
    fluxerr_2FWHM = []

    flux_star1_1FWHM_recov = []
    flux_star1_2FWHM_recov = []
    x_star1_recov = []
    y_star1_recov = []
    SN_star1_recov = []

    flux_star2_1FWHM_recov = []
    flux_star2_2FWHM_recov = []
    x_star2_recov = []
    y_star2_recov = []
    SN_star2_recov = []

    flux_star3_1FWHM_recov = []
    flux_star3_2FWHM_recov = []
    x_star3_recov = []
    y_star3_recov = []
    SN_star3_recov = []

    flux_bkg_1FWHM_recov = []
    flux_bkg_2FWHM_recov = []
    x_bkg_recov = []
    y_bkg_recov = []

    mjd = []
    utc = []
    exptime = []

    bar = progressbar.ProgressBar(maxval=data.shape[0]).start()

    print(data.shape[0])

    for i in range(data.shape[0]):
        bar.update(i)
        ind += [i+1]
        image = data[i]
        bkg = sep.Background(image)
        image_sub = image - bkg

        mean, median, std = sigma_clipped_stats(image_sub, sigma=3.0)
        daofind = DAOStarFinder(fwhm=6.0, threshold=5.*std)
        sources = daofind(image_sub)
        sources.sort('flux')
        sources.reverse()

        mask = (sources['xcentroid']-wcs_head['CRPIX1'])**2 +\
            (sources['ycentroid']-wcs_head['CRPIX2'])**2 < 100
        sources = sources[mask]

        if len(sources)!=0:
            wcs_copy = copy.copy(wcs_head)
            wcs_copy['CRPIX1']=sources[0]['xcentroid']
            wcs_copy['CRPIX2']=sources[0]['ycentroid']

        wcs = WCS(wcs_copy)

        # Start at frame index=1 (as opposed to 0) because we skip the first
        # blank image frame in the _reduced.fits data
        curr_time = start_utc + (i+1) * time_per_frame

        mjd.append(curr_time.mjd)
        utc.append(curr_time.datetime.strftime('%Y-%m-%d %H:%M:%S.%f'))
        exptime.append(time_per_exposure)

        x_FRB, y_FRB = utils.skycoord_to_pixel(frb, wcs, origin=1)

        flux_2FWHM_, fluxerr_2FWHM_, flag_2FWHM_ = sep.sum_circle(image_sub, [x_FRB], [y_FRB], 2*FWHM, err=bkg.globalrms, gain=1.0)
        flux_1FWHM_, fluxerr_1FWHM_, flag_1FWHM_ = sep.sum_circle(image_sub, [x_FRB], [y_FRB], FWHM, err=bkg.globalrms, gain=1.0)

        flux_1FWHM += [flux_1FWHM_[0]]
        flux_2FWHM += [flux_2FWHM_[0]]
        fluxerr_1FWHM += [fluxerr_1FWHM_[0]]
        fluxerr_2FWHM += [fluxerr_2FWHM_[0]]

        x_star1, y_star1 = utils.skycoord_to_pixel(star1, wcs, origin=1)

        mask1 = ut.create_mask(x_star1, y_star1, sizebox)
        flux_star1_1FWHM_recov_, x_star1_recov_, y_star1_recov_, SN_star1_recov_ = ut.bright_star_parameters(image_sub, mask1, FWHM)
        flux_star1_2FWHM_recov_, x_star1_recov_, y_star1_recov_, SN_star1_recov_ = ut.bright_star_parameters(image_sub, mask1, 2*FWHM)

        flux_star1_1FWHM_recov += [flux_star1_1FWHM_recov_[0]]
        flux_star1_2FWHM_recov += [flux_star1_2FWHM_recov_[0]]
        x_star1_recov += [x_star1_recov_[0]]
        y_star1_recov += [y_star1_recov_[0]]
        SN_star1_recov += [SN_star1_recov_[0]]

        x_star2, y_star2 = utils.skycoord_to_pixel(star2, wcs, origin=1)

        mask2 = ut.create_mask(x_star2, y_star2, sizebox)
        flux_star2_1FWHM_recov_, x_star2_recov_, y_star2_recov_, SN_star2_recov_ = ut.bright_star_parameters(image_sub, mask2, FWHM)
        flux_star2_2FWHM_recov_, x_star2_recov_, y_star2_recov_, SN_star2_recov_ = ut.bright_star_parameters(image_sub, mask2, 2*FWHM)

        flux_star2_1FWHM_recov += [flux_star2_1FWHM_recov_[0]]
        flux_star2_2FWHM_recov += [flux_star2_2FWHM_recov_[0]]
        x_star2_recov += [x_star2_recov_[0]]
        y_star2_recov += [y_star2_recov_[0]]
        SN_star2_recov += [SN_star2_recov_[0]]

        x_star3, y_star3 = utils.skycoord_to_pixel(star3, wcs, origin=1)

        mask3 = ut.create_mask(x_star3, y_star3, sizebox)
        flux_star3_1FWHM_recov_, x_star3_recov_, y_star3_recov_, SN_star3_recov_ = ut.bright_star_parameters(image_sub, mask3, FWHM)
        flux_star3_2FWHM_recov_, x_star3_recov_, y_star3_recov_, SN_star3_recov_ = ut.bright_star_parameters(image_sub, mask3, 2*FWHM)

        flux_star3_1FWHM_recov += [flux_star3_1FWHM_recov_[0]]
        flux_star3_2FWHM_recov += [flux_star3_2FWHM_recov_[0]]
        x_star3_recov += [x_star3_recov_[0]]
        y_star3_recov += [y_star3_recov_[0]]
        SN_star3_recov += [SN_star3_recov_[0]]

        x_bkg, y_bkg = utils.skycoord_to_pixel(bkg_coord, wcs, origin=1)

        flux_bkg_2FWHM_, fluxerr_bkg_2FWHM_, flag_bkg_2FWHM_ = sep.sum_circle(image_sub, [x_bkg], [y_bkg], 2*FWHM, err=bkg.globalrms, gain=1.0)
        flux_bkg_1FWHM_, fluxerr_bkg_1FWHM_, flag_bkg_1FWHM_ = sep.sum_circle(image_sub, [x_bkg], [y_bkg], FWHM, err=bkg.globalrms, gain=1.0)

        flux_bkg_1FWHM_recov += [flux_bkg_1FWHM_[0]]
        flux_bkg_2FWHM_recov += [flux_bkg_2FWHM_[0]]

    bar.finish()

    tab = Table()
    tab['frame'] = ind
    tab['file'] = [os.path.basename(basefile)]*len(flux_1FWHM)
    tab['UTC'] = utc # UTC at start of exposure
    tab['MJD'] = mjd # MJD at start of exposure
    tab['exptime'] = exptime

    tab['flux_1FWHM'] = flux_1FWHM
    tab['fluxerr_1FWHM'] = fluxerr_1FWHM
    tab['flux_2FWHM'] = flux_2FWHM
    tab['fluxerr_2FWHM'] = fluxerr_2FWHM

    tab['flux_star1_1FWHM'] = flux_star1_1FWHM_recov
    tab['flux_star1_2FWHM'] = flux_star1_2FWHM_recov
    tab['x_star1'] = x_star1_recov
    tab['y_star1'] = y_star1_recov
    tab['SN_star1'] = SN_star1_recov

    tab['flux_star2_1FWHM'] = flux_star2_1FWHM_recov
    tab['flux_star2_2FWHM'] = flux_star2_2FWHM_recov
    tab['x_star2'] = x_star2_recov
    tab['y_star2'] = y_star2_recov
    tab['SN_star2'] = SN_star2_recov

    tab['flux_star3_1FWHM'] = flux_star3_1FWHM_recov
    tab['flux_star3_2FWHM'] = flux_star3_2FWHM_recov
    tab['x_star3'] = x_star3_recov
    tab['y_star3'] = y_star3_recov
    tab['SN_star3'] = SN_star3_recov

    tab['flux_bkg_1FWHM'] = flux_bkg_1FWHM_recov
    tab['flux_bkg_2FWHM'] = flux_bkg_2FWHM_recov

    if master is None:
        master = tab
    else:
        master = vstack([master, tab])

    tabfile = f'count_table_{date}_{n}_{camera}.fits'
    tabfile = os.path.join('..','..','Data',tabfile)
    print(f'Writing out {tabfile}')
    tab.write(tabfile, overwrite=True)

masterfile = os.path.join('..','..','Data',f'master_table_{date}_{camera}.fits')
master.write(masterfile, overwrite=True)
