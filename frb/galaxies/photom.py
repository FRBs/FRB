""" Methods related to galaxy photometry """

import os
import warnings
import numpy as np

from pkg_resources import resource_filename

from IPython import embed

from astropy.table import Table, hstack, vstack, join
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units

from frb.galaxies import nebular

# Photometry globals
table_format = 'ascii.fixed_width'

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
    old_coords = SkyCoord(ra=old_tbl['ra'], dec=new_tbl['dec'], unit='deg')
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
        embed(header='50 of photom')  # Needs testing
        merge_tbl = hstack([old_tbl, new_tbl[match]])
        merge_tbl = hstack([merge_tbl, new_tbl[~match]])
    # Return
    return merge_tbl


def extinction_correction(filter, EBV, RV=3.1):
    """

    calculate MW extinction correction for given filter

    Args:
        filter (str):
            filter name (name of file without .dat extension)
        EBV (float):
            E(B-V) (can get from frb.galaxies.nebular.get_ebv which uses IRSA Dust extinction query
        RV:
            from gbrammer/threedhst eazyPy.py -- characterizes MW dust

    Returns:
             float: linear extinction correction

    """
    # Read in filter in Table
    path_to_filters = os.path.join(resource_filename('frb', 'data'), 'analysis', 'CIGALE')
    filter_file = os.path.join(path_to_filters, filter+'.dat')
    filter_tbl = Table.read(filter_file, format='ascii')

    #get wave and transmission (file should have these headers in first row)
    wave = filter_tbl['col1'].data
    throughput = filter_tbl['col2'].data

    #get MW extinction correction
    AV = EBV * RV
    AlAV = nebular.load_extinction('MW')
    Alambda = AV * AlAV(wave)
    source_flux = 1.

    #calculate linear correction
    delta = np.trapz(throughput * source_flux * 10 ** (-0.4 * Alambda), wave) / np.trapz(
        throughput * source_flux, wave)

    correction = 1./delta

    return correction
