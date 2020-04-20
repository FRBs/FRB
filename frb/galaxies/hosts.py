""" Module related to host galaxies of FRBs
Warning: Might get chopped up into pieces sommeday
"""
import numpy as np
import pdb

from astropy import units
from astropy.coordinates import SkyCoord, match_coordinates_sky


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


def random_separation(catalog, wcs, npix, trim=1*units.arcmin, ntrial=100):

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
