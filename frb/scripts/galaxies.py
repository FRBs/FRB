#!/usr/bin/env python
"""
Script to fuss a bit with FRB galaxies
Requres the specDB
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to fuss with FRB galaxies [v1.1]')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1 or FRB name (FRB180924)")
    parser.add_argument("--rho", default=300., type=float, help="Maximum impact parameter in kpc [default=300.]")
    parser.add_argument("--ang_offset", type=float, help="Maximum offset in arcsec [over-rides --rho if set]")
    parser.add_argument("--cat", default=False, action="store_true", help="Only show data from the catalog (not meta)")
    parser.add_argument("--specdb", type=str, help="specDB file; defaults to $SPECDB/FRB_specdb.hdf5")
    parser.add_argument("-p","--plot", default=False, action="store_true", help="Launch a plotting GUI?")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs

def plot_spec(casbah_gspec, meta, frb=None):
    import numpy as np
    from astropy.coordinates import SkyCoord
    from linetools import utils as ltu

    """ Plot galaxy spectra """
    # Grab spectra
    spec = casbah_gspec.spectra_from_meta(meta)
    # Add labels
    lbls = ['None']*len(meta)
    ugroups = np.unique(meta['GROUP'])
    for ugroup in ugroups:
        rows = np.where(meta['GROUP'] == ugroup)[0]
        coords = SkyCoord(ra=meta['RA_GROUP'][rows], dec=meta['DEC_GROUP'][rows], unit='deg')
        if frb is None:
            for kk, row, imeta in zip(range(len(rows)), rows, meta):
                # JNAME
                jname = ltu.name_from_coord(coords[kk])
                lbls[row]='{:s}_{:s}'.format(jname, meta['INSTR'][row])
        else:
            seps = frb.coord.separation(coords).to('arcsec')
            pas = frb.coord.position_angle(coords).to('deg')
            for row, sep, pa, imeta in zip(rows, seps, pas, meta):
                # Separation and PA
                lbls[row]='{:s}_{:s}_{:d}_{:d}'.format(frb.frb_name, meta['INSTR'][row],
                                                   int(pa.value), int(sep.value))
    spec.labels = lbls
    # Add object type
    spec.stypes = ['Galaxy']*len(lbls)
    # Add redshift
    spec.z = meta['zem_GROUP']
    # Plot
    spec.plot(xspec=True)#, scale=1.5)


def main(pargs):
    """ Run
    """
    import warnings
    import numpy as np

    from astropy import units

    from frb import frb
    from frb.galaxies import utils as gutils

    from linetools.scripts.utils import coord_arg_to_coord

    try:
        from specdb import group_utils
    except ImportError:
        print("You need to first install specdb")
        return

    # Load it up
    specDB = gutils.load_specdb(specdb_file=pargs.specdb)
    if specDB is None:
        return

    if pargs.coord[0:3] == 'FRB':
        frb = frb.FRB.by_name(pargs.coord)
        icoord = frb.coord.ra.value, frb.coord.dec.value
    else:
        icoord = coord_arg_to_coord(pargs.coord)

    # Meta?
    if pargs.cat:
        if pargs.ang_offset is not None:
            _, cat, _ = specDB.qcat.query_position(icoord, pargs.ang_offset*units.arcsec)
        else:
            _, cat, _ = specDB.qcat.query_position(icoord, pargs.rho*units.kpc)
        # Show
        ckeys = ['RA', 'DEC', 'zem', 'ZQ']
        group_utils.show_group_meta(cat, show_all_keys=False, meta_keys=ckeys)
    else: # Meta
        if pargs.ang_offset is not None:
            meta = specDB.meta_from_position(icoord, pargs.ang_offset*units.arcsec)
        else:
            meta = specDB.meta_from_position(icoord, pargs.rho*units.kpc)

        # Keys
        mkeys = ['GROUP', 'RA_GROUP', 'DEC_GROUP', 'zem_GROUP', 'ZQ']
        # Show
        if meta is None:
            print("No source found, try another location or a larger tolerance.")
            return
        else:
            group_utils.show_group_meta(meta, show_all_keys=False,
                                        meta_keys=mkeys)#, max_lines=10000000)

    # Plot
    if pargs.plot:
        if pargs.cat:
            print("Cannot mix --plot with --cat.  Try again!")
            return
        plot_spec(specDB, meta, frb=frb)


