""" Module for building specdb files related to FRB spectra
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
from IPython import embed
import h5py

import glob

from pkg_resources import resource_filename

from specdb import defs
from specdb.build import privatedb as pbuild
from specdb.build import utils as spbu

all_instruments = ['SDSS']

def grab_files(all_files, instrument):
    # Setup
    base_files = [os.path.basename(ifile) for ifile in all_files]
    file_subset = []
    # Simple loop 
    for kk, ifile in enumerate(base_files):
        if instrument in ifile:
            file_subset.append(all_files[kk])
    # Return
    return file_subset
    
    
def generate_by_refs(input_refs, outfile):
    # Not elegant but it works
    spectra_path = resource_filename('frb', '../Spectra')
    all_folders = glob.glob(spectra_path+'/*/*')
    all_refs = [os.path.basename(ifolder) for ifolder in all_folders]
    
    # Loop in input refs
    all_spec_files = []
    for ref in input_refs:
        idx = all_refs.index(ref)
        # Grab the list of spectra
        specs = glob.glob(os.path.join(all_folders[idx], '*.fits'))
        # Save
        all_spec_files += specs
    
    # Get it started
    # HDF5 file
    hdf = h5py.File(outfile, 'w')

    # Defs
    zpri = defs.z_priority()

    # Main DB Table
    id_key = 'FRB_ID'
    maindb, tkeys = spbu.start_maindb(id_key)
    gdict = {}

    # Loop on Instruments
    for instr in all_instruments:
        fits_files = grab_files(all_spec_files, instr)
        if len(ifiles) == 0:
            continue
        # Option dicts
        mwargs = {}
        mwargs['toler'] = 1.0 * u.arcsec  # Require an 
        skipz = False
        # Meta
        parse_head, mdict, fname = None, None, True
        if instr == 'SDSS':  # Spectra with Continua only
            mdict = dict(DISPERSER='BOTH', R=2000., TELESCOPE='SDSS 2.5-M', INSTR='SDSS')
            parse_head = {'DATE-OBS': 'DATE-OBS'}
            maxpix = 4000

        # Meta
        full_meta = pbuild.mk_meta(fits_files, ztbl, mdict=mdict, fname=fname,
                                   verbose=True, parse_head=parse_head, skip_badz=skipz,
                                   chkz=True, **mwargs)


        embed()

if __name__ == '__main__':

    # Test
    generate_by_refs(['DR7'], 'tst_specdb.hdf5')
