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
from frb.galaxies import defs

try:
    import extinction
except ImportError:
    print("extinction package not loaded.  Extinction corrections will fail")

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
        embed(header='50 of photom')  # Needs testing
        merge_tbl = hstack([old_tbl, new_tbl[match]])
        merge_tbl = hstack([merge_tbl, new_tbl[~match]])
    # Return
    return merge_tbl


def extinction_correction(filter, EBV, RV=3.1, max_wave=None):
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


def correct_photom_table(photom, EBV, max_wave=None):
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

    """

    # Dust correct
    for key in photom.keys():
        if key in ['Name', 'ra', 'dec', 'extinction', 'SDSS_ID',
                   'run', 'rerun'] or 'err' in key:
            continue
        filter = key
        if filter not in defs.valid_filters:
            print("Assumed filter {} is not in our valid list.  Skipping extinction".format(filter))
            continue
        # SDSS
        if 'SDSS' in filter:
            if 'extinction_{}'.format(filter[-1]) in photom.keys():
                print("Appying SDSS-provided extinction correction")
                photom[key] -= photom['extinction_{}'.format(filter[-1])]
                continue
        # Hack for LRIS
        if 'LRIS' in filter:
            _filter = 'LRIS_{}'.format(filter[-1])
        else:
            _filter = filter
        try:
            dust_correct = extinction_correction(_filter, EBV, max_wave=max_wave)
            mag_dust = 2.5 * np.log10(1. / dust_correct)
            photom[key] += mag_dust
        except:
            embed(header='145 of photom; bad filter?')
