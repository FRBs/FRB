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
    Merge photometry tables

    Args:
        new_tbl (astropy.table.Table):
            New table of photometry
        old_file (str):
            Path to the old table

    Returns:
        astropy.table.Table:
            Merged tables

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

    # Match?
    if np.sum(match) == len(new_coords):
        merge_tbl = join(old_tbl.filled(-999.), new_tbl, join_type='left')
    elif np.sum(match) == 0:
        embed(header='47 of photom')  # Needs testing
        merge_tbl = vstack([old_tbl, new_tbl])
    else:
        embed(header='50 of photom')  # Needs testing
        merge_tbl = hstack([old_tbl, new_tbl[match]])
        merge_tbl = hstack([merge_tbl, new_tbl[~match]])
    # Return
    return merge_tbl


