""" Module related to host galaxies of FRBs
Warning: Might get chopped up into pieces sommeday
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from astropy import units
from astropy.coordinates import SkyCoord, match_coordinates_sky

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
