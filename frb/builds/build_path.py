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

from frb.frb import FRB

from frb.associate import frbassociate
from frb.associate import frbs
from frb.galaxies import hosts

from frb import utils
import pandas

db_path = os.getenv('FRB_GDB')
if db_path is None:
    raise IOError('You need to have GDB!!')

def run(frb_list:list, host_coords:list, prior:dict, override:bool=False,
        tol:float=1.0):
    """Main method for generating a Host JSON file

    Args:
        frb_list (list): List of FRB names from the database
        host_coords (list): List of host galaxy coords fom the database
        prior (dict):
            Prior for PATH
        override (bool, optional): Attempt to over-ride errors. 
            Mainly for time-outs of public data. Defaults to False.
        tol (float, optional):  Tolearance for a match to the expected host
            in arcsec.

    Raises:
        e: [description]
        ValueError: [description]

    Returns:
        pandas.DataFrame:  Table of PATH values and a bit more
    """
    good_frb, PATH_O, PATH_Ox, RAs, Decs = [], [], [], [], []
    for frb, host_coord in zip(frb_list, host_coords):
        frb_name = utils.parse_frb_name(frb, prefix='frb')
        # Config
        if not hasattr(frbs, frb_name.lower()):
            print(f"PATH analysis not possible for {frb_name}")
            continue
        print(f"Performing PATH on {frb_name}")
        #
        config = getattr(frbs, frb_name.lower())
        config['skip_bayesian'] = True
        # Could do this outside the prior loop
        frbA = frbassociate.run_individual(config)

        # Setup for PATH

        # We skirt the usual candidate init
        frbA.candidates['mag'] = frbA.candidates[frbA.filter]
        frbA.init_cand_coords()
        # Set priors
        frbA.init_cand_prior('inverse', P_U=prior['U'])
        frbA.init_theta_prior(prior['theta']['method'], 
                                prior['theta']['max'])
        #frbA.calc_priors(prior['U'], method=prior['O'])
        #prior['theta']['ang_size'] = frbA.candidates.half_light.values
        #frbA.set_theta_prior(prior['theta'])

        # Localization
        frbA.init_localization('eellipse', 
                                center_coord=frbA.frb.coord,
                                eellipse=frbA.frb_eellipse)
        
        # Calculate priors
        frbA.calc_priors()                            

        # Calculate p(O_i|x)
        frbA.calc_posteriors('fixed', box_hwidth=frbA.max_radius)

        # Reverse Sort
        frbA.candidates = frbA.candidates.sort_values('P_Ox', ascending=False)

        # Check
        sep = frbA.candidates.coords.values[0].separation(host_coord).to('arcsec')
        # TODO - Remove this before the PR is merged
        try:
            assert sep < tol*units.arcsec, f'sep = {sep}'
        except:
            if frb_name == 'frb20191001':
                print("WE NEED TO FIX THAT IMAGE!!!")
                pass
            else:
                embed(header='106 of build')

        # Save for table
        good_frb.append(frb_name.upper())
        PATH_Ox.append(frbA.candidates.P_Ox.values[0])
        PATH_O.append(frbA.candidates.P_O.values[0])
        RAs.append(host_coord.ra.deg)
        Decs.append(host_coord.dec.deg)

    # Build the table
    df = pandas.DataFrame()
    df['FRB'] = good_frb
    df['RA'] = RAs
    df['Dec'] = Decs
    df['P(O)'] = PATH_O
    df['P(O|x)'] = PATH_Ox

    # 
    return df

def main(options:str=None, override:bool=False):
    """ Driver of the analysis

    Args:
        options (str, optional): [description]. Defaults to None.
        override (bool, optional): [description]. Defaults to False.
    """
    # Parse optionsd
    #if options is not None:
    #    if 'cigale' in options:
    #        build_cigale = True
    # Read public host table
    host_tbl = hosts.load_host_tbl()#hosts_file=hosts_file)

    host_coords = [SkyCoord(host_coord, frame='icrs') for host_coord in host_tbl.Coord.values]

    # Generate FRBs for PATH analysis
    frb_list = host_tbl.FRB.values.tolist()

    # Priors
    nhalf = 10.

    # Theta priors
    theta_max = 6.
    theta_u = dict(method='uniform', max=theta_max)
    theta_c = dict(method='core', max=theta_max)
    theta_e = dict(method='exp', max=theta_max)

    # Combined priors
    conservative = dict(theta=theta_u, O='identical', U=0, name='Conservative', nhalf=nhalf)
    adopted = dict(theta=theta_e, O='inverse', U=0., name='Adopted', nhalf=nhalf)


    results = run(frb_list, host_coords, adopted)

    # Write
    outfile = os.path.join(resource_filename('frb', 'data'), 'Galaxies', 'PATH.csv')
    results.to_csv(outfile)