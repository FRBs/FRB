""" Methods related to galaxy photometry """

import os
import warnings
import numpy as np

from pkg_resources import resource_filename

from IPython import embed

from astropy.io import fits
from astropy.table import Table, hstack, vstack, join
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units
from astropy.cosmology import Planck15 as cosmo
from astropy.wcs import utils as wcs_utils
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy import stats

from photutils import SkyCircularAperture
from photutils import aperture_photometry

from frb.galaxies import defs

try:
    import extinction
except ImportError:
    print("extinction package not loaded.  Extinction corrections will fail")

# Photometry globals
table_format = 'ascii.fixed_width'
fill_values_list = [('-999', '0'), ('-999.0', '0')]
fill_value = -999.

def merge_photom_tables(new_tbl, old_file, tol=1*units.arcsec, debug=False):
    """
    Merge photometry tables

    Args:
        new_tbl (astropy.table.Table):
            New table of photometry
        old_file (str or Table):
            Path to the old table

    Returns:
        astropy.table.Table:
            Merged tables

    """
    fill_value = -999.
    # File or tbl?
    if isinstance(old_file, str):
        # New file?
        if not os.path.isfile(old_file):
            return new_tbl
        # Load me
        old_tbl = Table.read(old_file, format=table_format)
    elif isinstance(old_file, Table):
        old_tbl = old_file
    else:
        embed(header='42 of photom')
    # Coords
    new_coords = SkyCoord(ra=new_tbl['ra'], dec=new_tbl['dec'], unit='deg')
    old_coords = SkyCoord(ra=old_tbl['ra'], dec=old_tbl['dec'], unit='deg')
    idx, d2d, _ = match_coordinates_sky(new_coords, old_coords, nthneighbor=1)
    match = d2d < tol

    # Match?
    if np.sum(match) == len(new_coords):
        # Insist on the same RA, DEC
        new_tbl['ra'] = old_tbl['ra'][idx[0]]
        new_tbl['dec'] = old_tbl['dec'][idx[0]]
        # Join
        merge_tbl = join(old_tbl.filled(-999.), new_tbl, join_type='left').filled(-999.)
    elif np.sum(match) == 0:
        merge_tbl = vstack([old_tbl, new_tbl]).filled(-999.)
    else:
        embed(header='50 of photom')  # Best to avoid!!  Use photom_by_name
    # Return
    return merge_tbl

def photom_by_name(name, filelist):
    """
    Generate a Table for a given galaxy from a list of photom files

    Warning:  Order matters!  Use best data last

    Args:
        name (str):
        filelist (list):

    Returns:
        astropy.table.Table:

    """
    # Loop on tables
    final_tbl = None
    for ifile in filelist:
        # Load an insure it is a masked Table
        tbl = Table(Table.read(ifile, format=table_format, fill_values=fill_values_list), masked=True)
        idx = tbl['Name'] == name
        if np.sum(idx) == 1:
            sub_tbl = tbl[idx]
            if final_tbl is None:
                final_tbl = sub_tbl
            else:
                for key in sub_tbl.keys():
                    if sub_tbl[key].mask != True:  # Cannot use "is"
                        final_tbl[key] = sub_tbl[key]
    # Return
    return final_tbl.filled(fill_value)


def extinction_correction(filter, EBV, RV=3.1, max_wave=None, required=True):
    """
    calculate MW extinction correction for given filter

    Uses the Fitzpatrick & Massa (2007) extinction law

    Args:
        filter (str):
            filter name (name of file without .dat extension)
        EBV (float):
            E(B-V) (can get from frb.galaxies.nebular.get_ebv which uses IRSA Dust extinction query
        RV:
            from gbrammer/threedhst eazyPy.py -- characterizes MW dust
        max_wave (float, optional):
            If set, cut off the calculation at this maximum wavelength.
            A bit of a hack for the near-IR, in large part because the
            MW extinction curve ends at 1.4 microns.
        required (bool, optional):
            Crash out if the transmission curve is not present

    Returns:
             float: linear extinction correction

    """
    # Read in filter in Table
    path_to_filters = os.path.join(resource_filename('frb', 'data'), 'analysis', 'CIGALE')
    # Hack for LRIS which does not differentiate between cameras
    if 'LRIS' in filter:
        _filter = 'LRIS_{}'.format(filter[-1])
    else:
        _filter = filter
    filter_file = os.path.join(path_to_filters, _filter+'.dat')
    if not os.path.isfile(filter_file):
        msg = "Filter {} is not in the Repo.  Add it!!".format(filter_file)
        if required:
            raise IOError(msg)
        else:
            warnings.warn(msg)
            return 1.
    filter_tbl = Table.read(filter_file, format='ascii')

    #get wave and transmission (file should have these headers in first row)
    wave = filter_tbl['col1'].data
    throughput = filter_tbl['col2'].data

    if max_wave:
        warnings.warn("Cutting off the extinction correction calculation at {} Ang".format(max_wave))
        gdwv = wave < max_wave
        wave = wave[gdwv]
        throughput = throughput[gdwv]

    #get MW extinction correction
    AV = EBV * RV
    #AlAV = nebular.load_extinction('MW')
    Alambda = extinction.fm07(wave, AV)
    source_flux = 1.

    #calculate linear correction
    delta = np.trapz(throughput * source_flux * 10 ** (-0.4 * Alambda), wave) / np.trapz(
        throughput * source_flux, wave)

    correction = 1./delta

    return correction


def correct_photom_table(photom, EBV, name, max_wave=None, required=True):
    """
    Correct the input photometry table for Galactic extinction
    Table is modified in place

    If there is SDSS photometry, we look for the extinction values
    provided by the Survey itself.

    Uses extinction_correction()

    Args:
        photom (astropy.table.Table):
        EBV (float):
            E(B-V) (can get from frb.galaxies.nebular.get_ebv which uses IRSA Dust extinction query
        name (str):\
            Name of the object to correct
        required (bool, optional):
            Crash out if the transmission curve is not present

    """
    # Cut the table
    mt_name = photom['Name'] == name
    if not np.any(mt_name):
        print("No matches to input name={}.  Returning".format(name))
        return
    elif np.sum(mt_name) > 1:
        raise ValueError("More than 1 match to input name={}.  Bad idea!!".format(name))
    idx = np.where(mt_name)[0][0]
    cut_photom = photom[idx]  # This is a Row

    # Dust correct
    for key in photom.keys():
        if key in ['Name', 'ra', 'dec', 'extinction', 'SDSS_ID',
                   'run', 'rerun'] or 'err' in key:
            continue
        filter = key
        if filter not in defs.valid_filters:
            print("Assumed filter {} is not in our valid list.  Skipping extinction".format(filter))
            continue
        # -999? -- Not even measured
        try:
            if cut_photom[filter] <= -999.:
                continue
        except:
            embed(header='187')
        # SDSS
        if 'SDSS' in filter:
            if 'extinction_{}'.format(filter[-1]) in photom.keys():
                print("Appying SDSS-provided extinction correction")
                cut_photom[key] -= cut_photom['extinction_{}'.format(filter[-1])]
                continue
        # Hack for LRIS
        if 'LRIS' in filter:
            _filter = 'LRIS_{}'.format(filter[-1])
        else:
            _filter = filter
        # Do it
        dust_correct = extinction_correction(_filter, EBV, max_wave=max_wave, required=required)
        mag_dust = 2.5 * np.log10(1. / dust_correct)
        cut_photom[key] += mag_dust
    # Add it back in
    photom[idx] = cut_photom

def sb_at_frb(host, cut_dat:np.ndarray, cut_err:np.ndarray, wcs:WCS, 
          fwhm=3., physical=False, min_uncert=2):
    """ Measure the surface brightness at an FRB location
    in a host galaxy

    Args:
        host (Host object): host galaxy object from frb repo
        cut_dat (np.ndarray): data (data from astorpy 2D Cutout object)
        cut_err (np.ndarray): inverse variance of data (from astropy 2D Cutout object)
        wcs (WCS): WCS for the cutout
        fwhm (float, optional): FWHM of the PSF of the image in either
            pixels or kpc. Defaults to 3 [pix].
        physical (bool, optional): If True, FWHM is in kpc. Defaults to False.
        min_uncert (int, optional): Minimum localization unceratainty
            for the FRB, in pixels.  Defaults to 2.

    Returns:
        tuple: sb_average, sb_average_err  [counts/sqarcsec]
    """
    # Generate the x,y grid of coordiantes
    x = np.arange(np.shape(cut_dat)[0])
    y = np.arange(np.shape(cut_dat)[1])
    xx, yy = np.meshgrid(x, y)
    coords = wcs_utils.pixel_to_skycoord(xx, yy, wcs)
    xfrb, yfrb = wcs_utils.skycoord_to_pixel(host.frb.coord, wcs)
    plate_scale = coords[0, 0].separation(coords[0, 1]).to('arcsec').value

    # Calculate total a, b uncertainty (FRB frame)
    uncerta, uncertb = host.calc_tot_uncert()

    # Put in pixel space
    uncerta /= plate_scale 
    uncertb /= plate_scale 

    # Set a minimum threshold
    uncerta = max(uncerta, min_uncert)
    uncertb = max(uncertb, min_uncert)
        
    # check if in ellipse -- pixel space!
    theta = host.frb.eellipse['theta']
    in_ellipse = ((xx - xfrb.item()) * np.cos(theta) + 
                  (yy - yfrb.item()) * np.sin(theta)) ** 2 / (uncerta ** 2) + (
                      (xx - xfrb.item()) * np.sin(theta) - (
                          yy - yfrb.item()) * np.cos(theta)) ** 2 / (uncertb ** 2) <= 1
    idx = np.where(in_ellipse)
    xval = xx[idx]
    yval = yy[idx]

    # x, y gal on the tilted grid (same for frb coords)
    xp = yval * np.cos(theta) - xval * np.sin(theta)
    yp = xval * np.cos(theta) + yval * np.sin(theta)

    xpfrb = yfrb.item() * np.cos(theta) - xfrb.item() * np.sin(theta)
    ypfrb = xfrb.item() * np.cos(theta) + yfrb.item() * np.sin(theta)

    # convert fwhm from pixels to arcsec or kpc to arcsec
    if physical:
        fwhm_as = fwhm * units.kpc * cosmo.arcsec_per_kpc_proper(host.z)
    else:
        fwhm_as = fwhm * plate_scale * units.arcsec

    # Aperture photometry at every pixel in the ellipse
    photom = []
    photom_var = []
    for i in np.arange(np.shape(idx)[1]):
        aper = SkyCircularAperture(coords[idx[0][i], idx[1][i]], fwhm_as)
        apermap = aper.to_pixel(wcs)

        # aperture photometry for psf-size within the galaxy
        photo_frb = aperture_photometry(cut_dat, apermap)
        photo_err = aperture_photometry(1 / cut_err, apermap)

        photom.append(photo_frb['aperture_sum'][0])
        photom_var.append(photo_err['aperture_sum'][0])

    # ff prob distribution
    p_ff = np.exp(-(xp - xpfrb) ** 2 / (2 * uncerta ** 2)) * np.exp(
        -(yp - ypfrb) ** 2 / (2 * uncertb ** 2))
    f_weight = (photom / (np.pi * fwhm_as.value ** 2)) * p_ff  # weighted photometry
    fvar_weight = (photom_var / (np.pi * fwhm_as.value ** 2)) * p_ff  # weighted sigma

    weight_avg = np.sum(f_weight) / np.sum(p_ff) # per unit area (arcsec^2)

    # Errors
    weight_var_avg = np.sum(fvar_weight) / np.sum(p_ff)
    weight_err_avg = np.sqrt(weight_var_avg)


    return weight_avg, weight_err_avg


def fractional_flux(cutout, frbdat, hg, nsig=3.):
    """Calculate the fractional flux at the FRB location

    Args:
        cutout (WCS Cutout2D): astropy 2D Cutout of data around host galaxy
        frbdat (frb.FRB): frb object loaded from frb repo
        hg (frb.galaxies.frbgalaxy.FRBHost): host galaxy object loaded from frb repo
        nsig (float, optional): sigma for FRB localization within which the measurement should be made. Defaults to 3.

    Returns:
        tuple: median_ff, sig_ff, ff_weight [no units]
    """

    # get image data from cutout
    cut_data = cutout.data
    frbcoord = frbdat.coord

    # shift the data to above zero (all positive values)
    shift_data = cut_data - np.min(cut_data)

    # make mesh grid
    if np.shape(cut_data)[0] != np.shape(cut_data)[1]:
        cut_data = np.resize(cut_data, (np.shape(cut_data)[1], np.shape(cut_data)[1]))
    x = np.arange(np.shape(cut_data)[0])
    y = np.arange(np.shape(cut_data)[1])
    xx, yy = np.meshgrid(x, y)
    coords = wcs_utils.pixel_to_skycoord(xx, yy, cutout.wcs)
    xfrb, yfrb = wcs_utils.skycoord_to_pixel(frbcoord, cutout.wcs)

    # Calc plate scale
    plate_scale = coords[0, 0].separation(coords[0, 1]).to('arcsec').value

    # get a, b, and theta from frb object -- convert to pixel space
    sig_a, sig_b = hg.calc_tot_uncert()
    # Put in pixel space
    sig_a /= plate_scale 
    sig_b /= plate_scale 

    # sigma
    a = nsig * sig_a
    if a < 1:
        print('a is less than 1!')
        a = 3

    b = nsig * sig_b
    if b < 1:
        print('b is less than 1!')
        b = 3


    # check if in ellipse -- pixel space!
    theta = hg.frb.eellipse['theta']
    in_ellipse = ((xx - xfrb.item()) * np.cos(theta).value + (yy - yfrb.item()) * np.sin(theta).value) ** 2 / (
            a ** 2) + (
                         (xx - xfrb.item()) * np.sin(theta).value - (yy - yfrb.item()) * np.cos(
                     theta).value) ** 2 / (
                         b ** 2) <= 1

    #print(frbdat.FRB, a, b, np.size(cut_data), np.size(cut_data[in_ellipse]))

    idx = np.where(in_ellipse)
    xval = xx[idx]
    yval = yy[idx]

    # x, y gal
    xp = yval * np.cos(theta).value - xval * np.sin(theta).value
    yp = xval * np.cos(theta).value + yval * np.sin(theta).value

    xpfrb = yfrb.item() * np.cos(theta).value - xfrb.item() * np.sin(theta).value
    ypfrb = xfrb.item() * np.cos(theta).value + yfrb.item() * np.sin(theta).value

    # sigma clip data to exclude background
    clipp = stats.sigma_clip(shift_data, sigma=1, maxiters=5)
    mask = np.ma.getmask(clipp)
    masked_dat = shift_data[mask]

    # fractional flux for all values in ellipse
    fprime_inlocal = []
    for dat in shift_data[idx]:
        fprime = np.sum(shift_data[shift_data < dat]) / np.sum(shift_data)
        fprime_inlocal.append(fprime)

    # ff prob distribution
    p_ff = np.exp(-(xp - xpfrb) ** 2 / (2 * a ** 2)) * np.exp(-(yp - ypfrb) ** 2 / (2 * b ** 2))
    f_weight = fprime_inlocal * p_ff  # weighted fractional fluxes

    avg_ff = np.sum(fprime_inlocal * p_ff) / np.sum(p_ff)
    var_ff = np.sum((fprime_inlocal - avg_ff) ** 2 * p_ff) / np.sum(p_ff)
    sig_ff = np.sqrt(var_ff)

    med_ff = np.percentile(f_weight, 50)
    l68, u68 = np.abs(np.percentile(f_weight, (16, 84)))

    # make array into list for writing out
    f_weight = np.array(f_weight).tolist()

    # return med_ff, med_flux, fprime_inlocal
    return med_ff, sig_ff, f_weight