""" Module related to host galaxies of FRBs
Warning: Might get chopped up into pieces sommeday
"""
import numpy as np
import pdb

import os
import pandas

from astropy import units
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
from scipy import interpolate
from scipy.integrate import quad
from pkg_resources import resource_filename

def chance_coincidence(rmag, r_i):
    """
    Calculate the chance probability of a galaxy to an FRB

    Taken from Bloom et al. 2002
        https://ui.adsabs.harvard.edu/abs/2002AJ....123.1111B/abstract

    ..todo.. Expand to allow for other filters

    Args:
        rmag (float):  r-band magnitude
        r_i (Angle or Quantity):
            Effective radius, angular
            Should be the max[2Rhalf, 3 sigma_r0, (R_0^2 + 4 Rhalf^2)^1/2]
            See Bloom et al. 2002

    Returns:
        float:  Probability of a chance association

    """
    # WHERE DOES THIS EQUATION COME FROM?
    sigma = 1. / (3600. ** 2 * 0.334 * np.log(10)) * 10 ** (0.334 * (rmag - 22.963) + 4.320)

    # Do it
    eta = np.pi * r_i.to('arcsec').value ** 2 * sigma
    Pch = 1. - np.exp(-eta)

    # Return
    return Pch


def chance_dx(rmag):
    """
    Returns the angular separation for a secure association (1%)
    as in https://ui.adsabs.harvard.edu/abs/2014MNRAS.437.1495T/abstract

    Args:
        rmag (float):
            r-band magnitude

    Returns:
        Quantity: Angular offset in arcsec

    """
    dx = 1.48 * 10**13 * rmag**(-9.53) * units.arcsec
    #
    return dx

def random_separation(catalog, wcs, npix, trim=1*units.arcmin, ntrial=100):
    """
    Find random offsets to

    Args:
        catalog (astropy.table.Table):
        wcs (astropy.WCS.WCS):
        npix (int):
        trim (astropy.units.Quantity, optional):
        ntrial (int, optional):

    Returns:
        astropy.units.Quantity: angular distances

    """

    # Catalog
    cat_coord = SkyCoord(ra=catalog['ra'], dec=catalog['dec'], unit='deg')

    # Trim
    bottom_corner = wcs.pixel_to_world(0, 0)
    bottom_offset = bottom_corner.directional_offset_by(-45.*units.deg, trim*np.sqrt(2))
    x0,y0 = [float(i) for i in wcs.world_to_pixel(bottom_offset)]

    top_corner = wcs.pixel_to_world(npix-1, npix-1)
    top_offset = top_corner.directional_offset_by(135.*units.deg, trim*np.sqrt(2))
    x1,y1 = [float(i) for i in wcs.world_to_pixel(top_offset)]


    # Generate a uniform grid
    ndim = int(np.sqrt(ntrial))

    xval = np.outer(np.linspace(x0, x1, ndim), np.ones(ndim))
    yval = np.outer(np.ones(ndim), np.linspace(y0, y1, ndim))

    # Coordinates now
    grid_coord = wcs.pixel_to_world(xval.flatten(), yval.flatten())

    # Match
    idx, d2d, d3d = match_coordinates_sky(grid_coord, cat_coord, nthneighbor=1)

    return d2d


def get_R(R_frb, R_0=0.2, R_h=0.25):
    """
    Calculates Radius of localisation region in arcsecond
    Based on Bloom et al 2002 and Eftekhari et al 2017
    
    Args:
        R_frb (float): The 1 sigma localization radius of the FRB
        R_0 (float): Radial angular separation between the FRB position and a presumed host
        R_h (float): Galaxy half light radius
    
    Returns:
        float: radius (in arcseconds)
    """

    return np.max([2*R_frb, np.sqrt(R_0**2 + 4*R_h**2)])


def read_r_mags(data_table_path):
    """
    Reads data used in Driver et al (2016).
    https://iopscience.iop.org/article/10.3847/0004-637X/827/2/108

    Args:
        data_table_path (string): Path to the fits file with data 

    Returns:
        array: r band magnitudes
        array: magnitude bin
        array: cosmic variance

    """
    table = fits.open(data_table_path)
    data = table[1].data
    r_mask = np.concatenate((np.where(data['Filtername'] == 'r')[0],
                             np.where(data['Filtername'] == 'F606W')[0]))
    _magbin = data['MagBinCentre'][r_mask]
    indexes = _magbin.argsort()
    magbin = _magbin[indexes]

    r_band_data = data['N(m)'][r_mask][indexes]
    cv = data['CosmicVariance'][r_mask][indexes]
    mag_uniq = np.unique(magbin)

    cvs = []
    r_dat = []
    for mag in mag_uniq:
        loc = np.where(magbin == mag)[0]
        if len(loc) > 1:
            cv_s = cv[loc]
            min_cv_loc = np.where(cv_s == cv_s.min())
            r_dat.append(r_band_data[loc[min_cv_loc]][0])
            cvs.append(cv_s[min_cv_loc][0])
        else:
            cvs.append(cv[loc][0])
            r_dat.append(r_band_data[loc][0])
    cvs = np.array(cvs)
    r_dat = np.array(r_dat)
    
    return r_dat, mag_uniq, cvs


def prob_eb17(R_frb, m, R_0=0.2, R_h=0.25, ret_numgal=False):
    """
    Calculates chance association probability of a galaxy to an FRB
    Taken from: 
        Section 2 of https://ui.adsabs.harvard.edu/abs/2017ApJ...849..162E/abstract

    Args:
        R_frb (float): The 1 sigma localization radius of the FRB in arcsec
        m (float): r band magnitude of the galaxy
        R_0 (float): Radial angular separation between the FRB position and a presumed host
        R_h (float): Galaxy half light radius
        ret_numgal (bool): to return the number of galaxies along with the chance coincidence probability
    
    Returns:
        float: Probability of chance coincidence

    """
    r_dat, mag_uniq, cvs = read_r_mags(resource_filename('frb',os.path.join('data','Galaxies','driver2016_t3data.fits')))
    spl = interpolate.UnivariateSpline(x=mag_uniq,
                                       y=np.log10(r_dat),
                                       bbox=[-100, 100],
                                       k=3)

    def n_gal(m_r):
        return 10 ** spl(m_r)

    num_dens_gal = quad(n_gal, 0, m)[0]
    R = get_R(R_frb, R_0, R_h)

    deg2arcsec = 60 * 60
    num_gals = np.pi * (R / deg2arcsec) ** 2 * num_dens_gal

    if ret_numgal:
        return 1 - np.exp(-1 * num_gals), num_gals
    else:
        return 1 - np.exp(-1 * num_gals)

def load_host_tbl(hosts_file:str=None, host_tbl:pandas.DataFrame=None):
    """ Generate a simple host table from a CSV, usually
    the public file

    Args:
        hosts_file (str, optional): [description]. Defaults to None.
        host_tbl ([type], optional): [description]. Defaults to None.

    Returns:
        pandas.DataFrame: [description]
    """
    galaxy_path = os.path.join(resource_filename('frb', 'data'), 
                               'Galaxies')
    if host_tbl is None:
        if hosts_file is None:
            hosts_file = os.path.join(galaxy_path, 'public_hosts.csv')
        host_tbl = pandas.read_csv(hosts_file)

    # Reformat a few columns
    sfrbs = [str(ifrb) for ifrb in host_tbl.FRB.values]
    host_tbl['FRB'] = sfrbs

    # Return
    return host_tbl