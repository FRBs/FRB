""" Top-level module to run PATH analysis on a list of
FRB host galaxies"""

from importlib.resources import files
import os

import pandas

from astropy.coordinates import SkyCoord

try:
    from astropath.priors import load_std_priors
except ModuleNotFoundError:
    print("astropath not installed; install it to builld PATH")
else:
    from frb.associate import frbassociate
    from frb.associate import frbs

from frb.frb import FRB
from frb.galaxies import hosts
from frb import utils

from IPython import embed

db_path = os.getenv('FRB_GDB')
if db_path is None:
    raise IOError('You need to have GDB!!')


def run(frb_list:list, 
        prior:dict, 
        write:bool=False,
        show:bool=False,
        add_skipped:bool=False,
        override:bool=False):
    """Main method for running PATH analysis for a list of FRBs

    Args:
        frb_list (list): List of FRB names from the database
        prior (dict):
            Prior for PATH
        write (bool, optional): Write the results to a CSV file. Defaults to False.
            Uses the --options write_indiv   option in the call
        show (bool, optional): Show the segmentation image. Defaults to False.
        add_skipped (bool, optional): Add skipped FRBs to the table. Defaults to False.
        override (bool, optional): Attempt to over-ride errors. 
            Mainly for time-outs of public data. Defaults to False.

    Raises:
        e: [description]
        ValueError: [description]

    Returns:
        pandas.DataFrame:  Table of PATH values and a bit more
    """
    good_frb, PATH_O, PATH_Ox, RAs, Decs = [], [], [], [], []
    ang_sizes, separations, sep_err = [], [], []
    skipped = []
    #for frb, host_coord in zip(frb_list, host_coords):
    for frb in frb_list:
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
            config, prior=iprior, show=show,
            posterior_method=config['posterior_method'])

        if frbA is None:
            print(f"PATH analysis not possible for {frb_name}")
            skipped.append(frb_name)
            continue

        # Save for table
        good_frb.append(frb_name.upper())
        PATH_Ox.append(frbA.candidates.P_Ox.values[0])
        PATH_O.append(frbA.candidates.P_O.values[0])
        RAs.append(frbA.candidates.ra.values[0])
        Decs.append(frbA.candidates.dec.values[0])
        #RAs.append(host_coord.ra.deg)
        #Decs.append(host_coord.dec.deg)
        ang_sizes.append(frbA.candidates.ang_size.values[0])
        separations.append(frbA.candidates.separation.values[0])

        # Write the full candidate table for each FRB?
        if write:
            # Add P_U
            frbA.candidates['P_U'] = prior['U']
            # Drop coords
            frbA.candidates.drop(columns=['coords'], inplace=True)
            outfile = os.path.join(files('frb'), 'data', 'Galaxies', 
                                   'PATH', f'{frb_name.upper()}_PATH.csv')
            #embed(header="build_path.py: 99")
            frbA.candidates.to_csv(outfile)
            print(f"PATH analysis written to {outfile}")

    # Add in the skipped??
    if add_skipped:
        good_idx = 0
        for ss, frb in enumerate(frb_list):
            frb_name = utils.parse_frb_name(frb, prefix='frb').upper()
            if frb_name not in good_frb:
                # Insert the values in order of the frb_list
                good_frb.insert(good_idx, frb_name)
                PATH_O.insert(good_idx, 0.)
                PATH_Ox.insert(good_idx, 0.)
                RAs.insert(good_idx, 0.)
                Decs.insert(good_idx, 0.)
                ang_sizes.insert(good_idx, 0.)
                separations.insert(good_idx, 0.)
            else:
                good_idx = good_frb.index(frb_name) + 1

    # Build the table
    df = pandas.DataFrame()
    df['FRB'] = good_frb
    df['RA'] = RAs
    df['Dec'] = Decs
    df['ang_size'] = ang_sizes
    df['P_O'] = PATH_O
    df['P_Ox'] = PATH_Ox
    df['P_U'] = PATH_Ox
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
            If None, will read the public host table and run on those
    """
    # Read public host table
    if frb is None:
        host_tbl = hosts.load_host_tbl()
        #host_coords = [SkyCoord(host_coord, frame='icrs') for host_coord in host_tbl.Coord.values]

        # Generate FRBs for PATH analysis
        frb_list = host_tbl.FRB.values.tolist()
    else:
        # Generate the list
        frb_list = frb.split(',')

    # Load prior
    priors = load_std_priors()
    prior = priors['adopted'] # Default

    # Parse optionsd
    write_indiv = False
    show = False
    add_skipped = False
    if options is not None:
        if 'new_prior' in options:
            theta_new = dict(method='exp', 
                             max=priors['adopted']['theta']['max'], 
                             scale=0.5)
            prior['theta'] = theta_new
            print("Using new prior with scale=0.5")
        if 'write_indiv' in options:
            write_indiv = True
        if 'show' in options:
            show = True
        if 'add_skipped' in options:
            add_skipped = True
        if 'PU=' in options:
            for sopt in options.split(','):
                if 'PU=' in sopt:
                    prior['U'] = float(sopt.split('=')[1])
        #embed(header="build_path.py: 155")
        

    results = run(frb_list, prior, write=write_indiv, show=show,
                  add_skipped=add_skipped)
    # Write
    outfile = os.path.join(files('frb'),'data','Galaxies','PATH','tmp.csv')
    results.to_csv(outfile)
    print(f"PATH analysis written to {outfile}")
    print("Rename it, push to Repo, and edit the PATH/README file accordingly")

    return results