""" Utilities related to FRB galaxies"""

import os
import glob
from IPython import embed
from pkg_resources import resource_filename
import numpy as np

try:
    from specdb.specdb import SpecDB
except ImportError:
    flg_specdb = False
else:
    flg_specdb = True

from astropy.coordinates import SkyCoord
import pandas as pd

from frb import frb

def load_specdb(specdb_file=None):
    """
    Automatically load the specDB file from $SPECDB/FRB_specDB.hdf5

    Args:
        specdb_file (str, optional):
            Over-ride the default file

    Returns:
        specdb.specdb.SpecDB:

    """
    if not flg_specdb:
        raise IOError("You must install the specdb package first!")
        return
    if specdb_file is None:
        if os.getenv('SPECDB') is None:
            raise IOError("You must set the SPECDB environmental variable")
        specdb_files = glob.glob(os.path.join(os.getenv('SPECDB'), 'FRB_specDB_*.hdf5'))
        if len(specdb_files) > 0:
            specdb_file = specdb_files[0]
            print("Loading spectra from {:s}".format(specdb_file))
        else:
            raise IOError("There are no FRB_specdb.hdf5 files in your SPECDB folder")
    # Load it up
    specDB = SpecDB(db_file=specdb_file)
    # Return
    return specDB


def list_of_hosts():
    """
    Scan through the Repo and generate a list of FRB Host galaxies

    Also returns a list of the FRBs

    Returns:
        list, list:

    """
    # FRB files
    frb_data = resource_filename('frb', 'data')
    frb_files = glob.glob(os.path.join(frb_data, 'FRBs', 'FRB*.json'))
    frb_files.sort()

    hosts = []
    frbs = []
    for ifile in frb_files:
        # Parse
        name = ifile.split('.')[-2]
        ifrb = frb.FRB.by_name(name)
        host = ifrb.grab_host()
        if host is not None:
            hosts.append(host)
            frbs.append(ifrb)
    # Return
    return frbs, hosts


def build_table_of_hosts(hosts=None):
    """
    Generate a Pandas table of FRB Host galaxy data.  These are slurped
    from the 'derived', 'photom', and 'neb_lines' dicts of each host object

    Warning:  As standard, missing values are given NaN in the Pandas table
        Be careful!

    Note:
        RA, DEC are given as RA_host, DEC_host to avoid conflict with the FRB table

    Returns:
        pd.DataFrame, dict:  Table of data on FRB host galaxies,  dict of their units

    """
    if hosts is None:
        _, hosts = list_of_hosts()
    nhosts = len(hosts)

    # Table
    host_tbl = pd.DataFrame({'Host': [host.name for host in hosts]})
    frb_names = ['FRB'+host.frb.frb_name for host in hosts]
    host_tbl['FRB'] = frb_names
    tbl_units = {}

    # Coordinates
    coords = SkyCoord([host.coord for host in hosts])
    # Named to faciliate merging with an FRB table
    host_tbl['RA_host'] = coords.ra.value
    host_tbl['DEC_host'] = coords.dec.value
    tbl_units['RA_host'] = 'deg'
    tbl_units['DEC_host'] = 'deg'

    # Loop on all the main dicts
    for attr in ['derived', 'photom', 'neb_lines','offsets','morphology','redshift']:
        # Load up the dicts
        dicts = [getattr(host, attr) for host in hosts]

        # Photometry
        all_keys = []
        for idict in dicts:
            all_keys += list(idict.keys())
            #all_keys += list(host.photom.keys())
        #
        all_keys = np.array(all_keys)
        uni_keys = np.unique(all_keys)

        # Slurp using Nan's for missing values
        tbl_dict = {}
        for key in uni_keys:
            tbl_dict[key] = np.array([np.nan]*nhosts)
        for ss in range(nhosts): #, host in enumerate(hosts):
            for pkey in dicts[ss].keys(): #host.photom.keys():
                tbl_dict[pkey][ss] = dicts[ss][pkey]
        for key in tbl_dict.keys():
            # Error check
            if key in host_tbl.keys():
                raise IOError("Duplicate items!!")
            # Set
            host_tbl[key] = tbl_dict[key]
            tbl_units[key] = 'See galaxies.defs.py'

    # Return
    return host_tbl, tbl_units


