""" Utilities related to FRB galaxies"""

import os
import glob
from IPython import embed

import importlib_resources
import numpy as np
from scipy.interpolate import interp1d
import warnings

import pandas

try:
    from specdb.specdb import SpecDB
except ImportError:
    flg_specdb = False
else:
    flg_specdb = True

from astropy.coordinates import SkyCoord
from astropy import units

import pandas as pd

import dust_extinction

from frb import frb

def deredden_spec(spectrum, ebv:float):
    """ Deredden the input spectrum using the input EBV value

    Args:
        spectrum (xspectrum1d.XSpectrum1D): Spectrum
        ebv (float): Galactic reddening

    Returns:
        xspectrum1d.XSpectrum1D: De-reddened spectrum
    """
    # linetools WILL BE DEPRECATED
    from linetools.spectra import xspectrum1d

    # Correct for Galactic extinction
    #   Need to replace it 
    AV = ebv * 3.1  # RV
    extmod = dust_extinction.parameter_averages.G23(Rv=3.1)
    AlAV = extmod(spectrum.wavelength)#*units.AA)
    Al = AlAV * AV
    #Al = extinction.ccm89(spectrum.wavelength.value, AV, 3.1)

    # New spec
    new_flux = spectrum.flux * 10**(Al/2.5)
    new_sig = spectrum.sig * 10**(Al/2.5)
    new_spec = xspectrum1d.XSpectrum1D.from_tuple((spectrum.wavelength, new_flux, new_sig))

    # Return
    return new_spec

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
            raise IOError("There are no FRB_specDB_*.hdf5 files in your SPECDB folder")
    # Load it up
    specDB = SpecDB(db_file=specdb_file)
    # Return
    return specDB


def list_of_hosts(skip_bad_hosts=True, verbose:bool=False):
    """
    Scan through the Repo and generate a list of FRB Host galaxies

    Also returns a list of the FRBs

    Args:
        skip_bad_hosts (bool):
        verbose (bool):
            If True, print more to the screen

    Returns:
        list, list:

    """
    # FRB files
    frb_data = importlib_resources.files('frb.data.FRBs')
    frb_files = glob.glob(str(frb_data/'FRB*.json'))
    frb_files.sort()

    hosts = []
    frbs = []
    host_names = []
    for ifile in frb_files:
        # Parse
        name = ifile.split('/')[-1].split('.')[-2]
        ifrb = frb.FRB.by_name(name)
        try:
            host = ifrb.grab_host(verbose=verbose)
        except AssertionError as e:
            if skip_bad_hosts:
                if verbose:
                    print(f"Skipping bad host of FRB {ifrb}")
                continue
            else:
                raise e
        if host is not None:
            hosts.append(host)
            frbs.append(ifrb)
    # Return
    return frbs, hosts


def build_table_of_hosts(attrs:list=None): #PATH_root_file:str='scale0.5.csv'):
    """
    Generate a Pandas table of FRB Host galaxy data.  These are slurped
    from the 'derived', 'photom', and 'neb_lines' dicts of each host object

    Warning:  As standard, missing values are given NaN in the Pandas table
        Be careful!

    Note:
        RA, DEC are given as RA_host, DEC_host to avoid conflict with the FRB table

    Args:

    Returns:
        pd.DataFrame, dict:  Table of data on FRB host galaxies,  dict of their units

    """
    frbs, hosts = list_of_hosts(verbose=False)
    nhosts = len(hosts)

    # Table
    host_tbl = pd.DataFrame({'Host': [host.name for host in hosts]})
    frb_names = [host.frb.frb_name for host in hosts]
    host_tbl['FRBname'] = frb_names
    tbl_units = {}

    # Coordinates
    coords = SkyCoord([host.coord for host in hosts])
    # Named to faciliate merging with an FRB table
    host_tbl['RA_host'] = coords.ra.value
    host_tbl['DEC_host'] = coords.dec.value
    tbl_units['RA_host'] = 'deg'
    tbl_units['DEC_host'] = 'deg'

    # FRBs
    host_tbl['FRBobj'] = frbs

    # Loop on all the main dicts
    if attrs is None:
        attrs = ['derived', 'photom', 'neb_lines','offsets', 'morphology','redshift']
    for attr in attrs:
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
            #tbl_dict[key] = np.array([np.nan]*nhosts)
            tbl_dict[key] = [np.nan]*nhosts

        # Fill in
        for ss in range(nhosts): #, host in enumerate(hosts):
            for pkey in dicts[ss].keys(): #host.photom.keys():
                tbl_dict[pkey][ss] = dicts[ss][pkey]

        # Now build the table
        #embed(header='170 of utils ')
        pd_tbl = pd.DataFrame(tbl_dict)
        host_tbl = pd.concat([host_tbl, pd_tbl], axis=1)

        for key in tbl_dict.keys():
            # Error check
            if key in tbl_units.keys():
                raise IOError("Duplicate items!!")
        #    # Set
        #    host_tbl[key] = tbl_dict[key]
            tbl_units[key] = 'See galaxies.defs.py'

    '''
    # Add PATH values
    path_tbl = load_PATH(PATH_root_file=PATH_root_file)
    path_coords = SkyCoord(ra=path_tbl.RA, dec=path_tbl.Dec, unit='deg')

    host_coords = SkyCoord(ra=host_tbl.RA_host, dec=host_tbl.DEC_host, unit='deg')
    '''

    # PATH 
    host_tbl['P_Ox'] = np.nan
    host_tbl['P_O'] = np.nan
    host_tbl['ang_size'] = np.nan

    hosts_file = importlib_resources.files('frb.data.Galaxies')/'public_hosts.csv'
    hosts_df = pandas.read_csv(hosts_file, index_col=False)

    # Loop
    for index, host_row in hosts_df.iterrows():
        '''
        # Match to table RA, DEC
        sep = host_coords.separation(path_coords[index])
        imin = np.argmin(sep)
        # REDUCE THIS TOL TO 1 arcsec!!
        print(f"Min sep = {sep[imin].to('arcsec')}")
        if sep[imin] < 1.0*units.arcsec:
            for key in ['P_Ox', 'P_O', 'ang_size']:
                host_tbl.loc[imin,key] = path_row[key]
        '''
        if 'HG'+host_row.FRB in host_tbl.Host.values:
            indx = host_tbl[host_tbl.Host == 'HG'+host_row.FRB].index[0]
            # Set PATH values
            host_tbl.loc[indx,'P_Ox'] = host_row.P_Ox
        #embed(header='236 of utils')

    # Return
    return host_tbl, tbl_units

def load_f_mL():
    """ Generate an interpolater from mag to Luminosity as 
    a function of redshift (up to z=4)

    Warning:  this is rather approximate

    Returns:
        scipy.interpolate.interp1d:

    """
    # Grab m(L) table
    data_file = importlib_resources.files('frb.data.Galaxies')/'galLF_vs_z.txt'
    df = pandas.read_table(data_file, index_col=False)

    # Interpolate
    f_mL = interp1d(df.z, df['m_r(L*)'])

    # Return
    return f_mL


def load_PATH(PATH_root_file:str='adopted.csv'):
    """Load up the PATH table

    Args:
        PATH_root_file (str, optional): [description]. Defaults to 'adopted.csv'.

    Returns:
        pandas.DataFrame: Table of galaxy coordinates and PATH results
    """
    path_file = importlib_resources.files('frb.data.Galaxies.PATH')/PATH_root_file
    path_tbl = pd.read_csv(path_file, index_col=False)

    return path_tbl