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
all_instruments = ['SDSS', 'FORS2', 'MUSE', 'KCWI', 'MagE', 'GMOS-S', 'LRISr', 
                   'DEIMOS', 'XSHOOTER', 'GMOS-N']
db_path = os.getenv('FRB_GDB')


def grab_files(all_files, refs_list, instrument):
    """
    Simple method to parse a set of files for a given an instrument

    Args:
        all_files (list):
            Complete list of files
        refs_list (list):
            List of references of these files
        instrument (str):
            Instrument name to parse on

    Returns:
        list, list: List of files and their references matching the input instrument

    """
    # Setup
    base_files = [os.path.basename(ifile) for ifile in all_files]
    file_subset = []
    ref_subset = []
    # Simple loop
    for kk, ifile in enumerate(base_files):
        if instrument in ifile:
            file_subset.append(all_files[kk])
            ref_subset.append(refs_list[kk])
    # Return
    return file_subset, ref_subset


def load_z_tables(path):
    """
    Load up a redshift table from the Galaxy_DB

    Redshift tables are those that begin with 'z'

    Args:
        path (str):
            Path to the folder holding one or more redshift tables.

    Returns:
        astropy.table.Table: Redshift table with RA, DEC, ZEM, ..

    """
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
    Enter the directory and build a redshift table
    based on the spectra present

    Returns:

    """
    embed(header='THIS NEEDS HELP')
    #
    all_folders = glob.glob(db_path+'/SDSS/*')
    for folder in all_folders:
        Jnames = []
        # Grab the list of spectra files
        spec_files = glob.glob(os.path.join(folder, 'J*_spec.fits'))
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
    """
    Build a specDB file according to the input references

    Args:
        input_refs (list):
            List of references from which to build the specDB
        outfile (str):
            Output filename
        version (str):
            Version number
    """
    # Not elegant but it works
    all_folders = glob.glob(db_path+'/*/*')
    all_refs = [os.path.basename(ifolder) for ifolder in all_folders]

    # z_tbl
    allz_tbl = Table()
    
    # Loop in input refs
    all_spec_files = []
    refs_list = []
    for ref in input_refs:
        idx = all_refs.index(ref)
        # Redshift tables
        z_tbl = load_z_tables(all_folders[idx])
        allz_tbl = vstack([allz_tbl, z_tbl])
        # Grab the list of spectra
        specs = glob.glob(os.path.join(all_folders[idx], 'J*_spec.fits'))
        if len(specs) == 0:
            continue
        # Save
        all_spec_files += specs
        refs_list += [ref]*len(specs)

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
    #pair_groups = ['MUSE']
    pair_groups = []
    badf = None
    for instr in all_instruments:
        print("Working on {}".format(instr))
        fits_files, irefs = grab_files(all_spec_files, refs_list, instr)
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
        elif instr == 'MagE':
            parse_head = {'R': True, 'DATE-OBS': 'MJD-OBS', 'TELESCOPE': 'TELESCOP',
                          'INSTR': 'INSTRUME', 'DISPERSER': 'DISPNAME'}
            maxpix = 18000
            scale = 1e-17
        elif instr == 'GMOS-S':
            mdict = dict(TELESCOPE='Gemini-S', INSTR='GMOS-S')
            parse_head = {'R': True, 'DATE-OBS': 'MJD-OBS', 'DISPERSER': 'DISPNAME'}
            maxpix = 3500
            scale = 1e-17
        elif instr == 'GMOS-N':
            mdict = dict(TELESCOPE='Gemini-N', INSTR='GMOS-N')
            parse_head = {'R': True, 'DATE-OBS': 'MJD-OBS', 
                          'DISPERSER': 'DISPNAME'}
            maxpix = 3500
            scale = 1e-17
        elif instr == 'LRISr':
            mdict = dict(TELESCOPE='Keck-1')
            parse_head = {'DATE-OBS': 'MJD', 'DISPERSER': 'DISPNAME', 'INSTR': 'INSTRUME'}
            maxpix = 2050
            scale = 1e-17
        elif instr == 'DEIMOS':
            mdict = dict(TELESCOPE='Keck-2')
            parse_head = {'DATE-OBS': 'MJD', 'DISPERSER': 'DISPNAME', 'INSTR': 'INSTRUME'}
            maxpix = 9000
            scale = 1e-17
        elif instr == 'XSHOOTER':
            mdict = dict(TELESCOPE='VLT')
            parse_head = {'DATE-OBS': 'MJD', 'DISPERSER': 'DISPNAME', 'INSTR': 'INSTRUME'}
            maxpix = 33000
            scale = 1e-17
        else:
            embed(header='172')
        

        # Meta
        full_meta = pbuild.mk_meta(fits_files, allz_tbl, mdict=mdict, fname=fname,
                                   verbose=True, parse_head=parse_head, skip_badz=skipz,
                                   stype='GAL',
                                   chkz=True, **mwargs)
        full_meta['Ref'] = irefs
        # Survey flag
        flag_g = spbu.add_to_group_dict(instr, gdict, skip_for_debug=True)
        # IDs
        if 'MUSE' in instr:
            embed(header='207 of build specdb')
        maindb = spbu.add_ids(maindb, full_meta, flag_g, tkeys, id_key, first=(flag_g==1), close_pairs=(instr in pair_groups))

        # Ingest --
        pbuild.ingest_spectra(hdf, instr, full_meta, max_npix=maxpix, verbose=False,
                              badf=badf, grab_conti=False, scale=scale, **swargs)

    # Write
    spbu.write_hdf(hdf, str('FRB'), maindb, zpri, gdict, version, Publisher=str('JXP'))
    print("Wrote {:s} DB file".format(outfile))
    print("You probably need to move it into SPECDB")


def main(inflg='all'):

    if inflg == 'all':
        flg = np.sum(np.array( [2**ii for ii in range(25)]))
    else:
        flg = int(inflg)

    # Public
    if flg & (2**0):
        generate_by_refs(['Prochaska2019', 'Bannister2019', 'Bhandari2019',
                          'Heintz2020', 'Simha2020', 'Tendulkar2017'],
                         'FRB_specDB_Public.hdf5', 'v0.4')


# Command line execution
if __name__ == '__main__':
    pass
