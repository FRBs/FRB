""" Module for building specdb files related to FRB spectra
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
from IPython import embed
import h5py

import glob

from pkg_resources import resource_filename

from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units
from astropy.table import Table, vstack

from specdb import defs
from specdb.build import privatedb as pbuild
from specdb.build import utils as spbu

from frb.surveys import sdss

# Globals
all_instruments = ['SDSS', 'FORS2', 'MUSE', 'KCWI']
spectra_path = resource_filename('frb', '../DB/Spectra')

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


def load_z_tables(path):
    z_files = glob.glob(os.path.join(path, 'z*'))

    z_tbl = Table()
    for z_file in z_files:
        # Load
        if False:
            pass
        else:
            itbl = Table.read(z_file, format='ascii.fixed_width')
        # Add RA, DEC
        if 'RA' in itbl.keys():
            pass
        elif 'JCOORD' in itbl.keys():
            coords = SkyCoord(itbl['JCOORD'], unit=(units.hourangle, units.deg))
            itbl['RA'] = coords.ra.value
            itbl['DEC'] = coords.dec.value
        # Append
        z_tbl = vstack([z_tbl, itbl])
    # Return
    return z_tbl


def sdss_redshifts():
    """
    Enter the SDSS directory and build a redshift table
    based on the spectra present

    Returns:

    """
    #
    all_folders = glob.glob(spectra_path+'/SDSS/*')
    for folder in all_folders:
        Jnames = []
        # Grab the list of spectra files
        spec_files = glob.glob(os.path.join(folder, 'J*.fits'))
        # Generate the name list
        Jnames += [os.path.basename(ifile).split('_')[0] for ifile in spec_files]
        # Coords
        coords = SkyCoord(Jnames, unit=(units.hourangle, units.deg))  # from DES
        # Setup
        done = np.zeros_like(coords.ra.value, dtype=bool)
        zs = np.zeros_like(coords.ra.value)

        # Loop me
        while np.any(~done):
            # Grab the first not done
            i0 = np.where(~done)[0][0]
            coord = coords[i0]
            # Grab the SDSS data
            sdssSurvey = sdss.SDSS_Survey(coord, 10*units.arcmin)
            #
            sub_coords = coords[~done]
            sep = coord.separation(sub_coords)
            doidx = np.where(sep < 10*units.arcmin)[0]
            dothem = coords[doidx]
            # Now match
            catalog = sdssSurvey.get_catalog()
            sdss_coords = SkyCoord(ra=catalog['ra'], dec=catalog['dec'], unit='deg')
            idx, d2d, d3d = match_coordinates_sky(dothem, sdss_coords, nthneighbor=1)
            # Fill
            zs[doidx] = catalog['z_spec'][idx]
            done[np.where(dothem)] = True

        # Write the catalog
        tbl = Table()
        tbl['RA'] = coords.ra.value
        tbl['DEC'] = coords.dec.value
        tbl['ZEM'] = zs
        tbl['ZEM_SOURCE'] = 'SDSS'
        tbl['ZQ'] = 4
        tbl.write(os.path.join(folder, 'z_SDSS.ascii'), overwrite=True, format='ascii.fixed_width')

    
def generate_by_refs(input_refs, outfile, version):
    # Not elegant but it works
    all_folders = glob.glob(spectra_path+'/*/*')
    all_refs = [os.path.basename(ifolder) for ifolder in all_folders]

    # z_tbl
    allz_tbl = Table()
    
    # Loop in input refs
    all_spec_files = []
    for ref in input_refs:
        idx = all_refs.index(ref)
        # Grab the list of spectra
        specs = glob.glob(os.path.join(all_folders[idx], '*.fits'))
        # Save
        all_spec_files += specs
        # Redshift tables
        z_tbl = load_z_tables(all_folders[idx])
        allz_tbl = vstack([allz_tbl, z_tbl])
    
    # Get it started
    # HDF5 file
    hdf = h5py.File(outfile, 'w')

    # Defs
    zpri = defs.z_priority()

    # Main DB Table
    id_key = 'FRB_ID'
    maindb, tkeys = spbu.start_maindb(id_key)
    tkeys += ['ZQ']
    gdict = {}

    # Loop on Instruments
    pair_groups = []
    badf = None
    for instr in all_instruments:
        fits_files = grab_files(all_spec_files, instr)
        if len(fits_files) == 0:
            continue
        # Option dicts
        mwargs = {}
        mwargs['toler'] = 1.0 * units.arcsec  # Require an
        skipz = False
        swargs = {}
        # Meta
        parse_head, mdict, fname = None, None, True
        if instr == 'SDSS':
            mdict = dict(DISPERSER='BOTH', R=2000., TELESCOPE='SDSS 2.5-M', INSTR='SDSS')
            parse_head = {'DATE-OBS': 'MJD'}
            maxpix = 4000
            scale = 1e-17
        elif instr == 'FORS2':
            mdict = dict(TELESCOPE='VLT', INSTR='FORS2')
            parse_head = {'DATE-OBS': 'MJD', 'DISPERSER': 'DISPNAME', 'R': True}
            maxpix = 2050
            scale = 1e-17
        elif instr == 'MUSE':
            mdict = dict(TELESCOPE='VLT', R=2000.)
            parse_head = {'DATE-OBS': 'MJD-OBS', 'DISPERSER': 'DISPNAME', 'INSTR': 'INSTRUME'}
            maxpix = 4000
            scale = 1e-20
        elif instr == 'KCWI':
            mdict = dict(TELESCOPE='Keck-2')
            parse_head = {'DATE-OBS': 'MJD', 'DISPERSER': 'DISPNAME', 'INSTR': 'INSTRUME', 'R': True}
            maxpix = 4000
            scale = 1e-17
        else:
            embed(header='172')

        # Meta
        full_meta = pbuild.mk_meta(fits_files, allz_tbl, mdict=mdict, fname=fname,
                                   verbose=True, parse_head=parse_head, skip_badz=skipz,
                                   stype='GAL',
                                   chkz=True, **mwargs)
        # Survey flag
        flag_g = spbu.add_to_group_dict(instr, gdict, skip_for_debug=True)
        # IDs
        maindb = spbu.add_ids(maindb, full_meta, flag_g, tkeys, id_key,
                              first=(flag_g==1), close_pairs=(instr in pair_groups))

        # Ingest --
        pbuild.ingest_spectra(hdf, instr, full_meta, max_npix=maxpix, verbose=False,
                              badf=badf, grab_conti=False, scale=scale, **swargs)

    # Write
    spbu.write_hdf(hdf, str('FRB'), maindb, zpri, gdict, version, Publisher=str('JXP'))
    print("Wrote {:s} DB file".format(outfile))

if __name__ == '__main__':

    # Test
    generate_by_refs(['DR7', 'Prochaska2019', 'Bannister2019'], 'CRAFT_specdb.hdf5', 'v0.1')
    #sdss_redshifts()
