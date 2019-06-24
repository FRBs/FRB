""" Methods related to galaxy photometry """

import os
import warnings
import numpy as np

from IPython import embed

from astropy.table import Table, hstack, vstack, join
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units

# Photometry globals
table_format = 'ascii.fixed_width'

def merge_photom_tables(new_tbl, old_file, tol=1*units.arcsec):
    """

    Args:
        new_tbl (astropy.table.Table):
        old_file (str):

    Returns:
        astropy.table.Table:

    """
    # New file?
    if not os.path.isfile(old_file):
        return new_tbl
    # Load me
    old_tbl = Table.read(old_file, format=table_format)
    # Coords
    new_coords = SkyCoord(ra=new_tbl['ra'], dec=new_tbl['dec'], unit='deg')
    old_coords = SkyCoord(ra=old_tbl['ra'], dec=new_tbl['dec'], unit='deg')
    idx, d2d, _ = match_coordinates_sky(new_coords, old_coords, nthneighbor=1)
    match = d2d < tol

    embed(header='39 of photom')
    # Match?
    if np.sum(match) == len(new_coords):
        merge_tbl = join(old_tbl, new_tbl)
    elif np.sum(match) == 0:
        merge_tbl = vstack([old_tbl, new_tbl])
    else:
        merge_tbl = hstack([old_tbl, new_tbl[match]])
        merge_tbl = hstack([merge_tbl, new_tbl[~match]])
    # Return
    return merge_tbl


