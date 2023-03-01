""" Top-level module to build or re-build the JSON files for
FRB host galaxies"""

from pkg_resources import resource_filename
import os
import sys

from IPython import embed

import numpy as np

import pandas

from astropy.coordinates import SkyCoord
from astropy import units

from astropath.priors import load_std_priors

from frb.frb import FRB

from frb.associate import frbassociate
from frb.associate import frbs
from frb.galaxies import hosts

from frb import utils


db_path = os.getenv('FRB_GDB')
if db_path is None:
    raise IOError('You need to have GDB!!')


def run(frb_list:list, host_coords:list, prior:dict, 
        override:bool=False):
    """Main method for running PATH analysis for a list of FRBs

    Args:
        frb_list (list): List of FRB names from the database
        host_coords (list): List of host galaxy coords fom the database
        prior (dict):
            Prior for PATH
        override (bool, optional): Attempt to over-ride errors. 
            Mainly for time-outs of public data. Defaults to False.
        tol (float, optional):  Tolerance for a match to the expected host
            in arcsec.

    Raises:
        e: [description]
        ValueError: [description]

    Returns:
        pandas.DataFrame:  Table of PATH values and a bit more
    """
    good_frb, PATH_O, PATH_Ox, RAs, Decs = [], [], [], [], []
    ang_sizes, separations, sep_err = [], [], []
    skipped = []
    for frb, host_coord in zip(frb_list, host_coords):
        frb_name = utils.parse_frb_name(frb, prefix='frb')
        # Config
        if not hasattr(frbs, frb_name.upper()):
            print(f"PATH analysis not possible for {frb_name}")
            skipped.append(frb_name)
            continue
        print(f"Performing PATH on {frb_name}")
        config = getattr(frbs, frb_name.upper())

        # Adjust prior, as needed
        iprior = prior.copy()
        if 'PU' in config.keys():
            iprior['U'] = config['PU']

        # Run me
        frbA = frbassociate.run_individual(
            config, prior=iprior, 
            posterior_method=config['posterior_method'])

        if frbA is None:
            print(f"PATH analysis not possible for {frb_name}")
            skipped.append(frb_name)
            continue

        # Save for table
        good_frb.append(frb_name.upper())
        PATH_Ox.append(frbA.candidates.P_Ox.values[0])
        PATH_O.append(frbA.candidates.P_O.values[0])
        RAs.append(host_coord.ra.deg)
        Decs.append(host_coord.dec.deg)
        ang_sizes.append(frbA.candidates.ang_size.values[0])
        separations.append(frbA.candidates.separation.values[0])

    # Build the table
    df = pandas.DataFrame()
    df['FRB'] = good_frb
    df['RA'] = RAs
    df['Dec'] = Decs
    df['ang_size'] = ang_sizes
    df['P_O'] = PATH_O
    df['P_Ox'] = PATH_Ox
    df['separation'] = separations

    for frb_name in skipped:
        print(f"PATH analysis not possible for {frb_name}")
    # 
    return df

def main(options:str=None, frb:str=None):
    """ Driver of the analysis

    Args:
        options (str, optional): [description]. Defaults to None.
        frb (str, optional): FRB name
    """
    # Read public host table
    if frb is None:
        host_tbl = hosts.load_host_tbl()

        host_coords = [SkyCoord(host_coord, frame='icrs') for host_coord in host_tbl.Coord.values]

        # Generate FRBs for PATH analysis
        frb_list = host_tbl.FRB.values.tolist()
    else:
        ifrb = FRB.by_name(frb)
        host = ifrb.grab_host()
        frb_list = [ifrb.frb_name]
        host_coords = [host.coord]

    # Load prior
    priors = load_std_priors()
    prior = priors['adopted'] # Default

    # Parse optionsd
    if options is not None:
        if 'new_prior' in options:
            theta_new = dict(method='exp', 
                             max=priors['adopted']['theta']['max'], 
                             scale=0.5)
            prior['theta'] = theta_new
            print("Using new prior with scale=0.5")

    results = run(frb_list, host_coords, prior)
    # Write
    outfile = os.path.join(resource_filename('frb', 'data'), 'Galaxies', 
                           'PATH', 'tmp.csv')
    results.to_csv(outfile)
    print(f"PATH analysis written to {outfile}")
    print("Rename it, push to Repo, and edit the PATH/README file accordingly")

    return results