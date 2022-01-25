""" Top-level module to build or re-build the JSON files for FRBs """

from pkg_resources import resource_filename
import os
import sys
import warnings

from IPython import embed

import numpy as np
import requests

import pandas

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table
from astropy.coordinates import match_coordinates_sky

from frb.frb import FRB, load_frb_data
from frb.galaxies import frbgalaxy, defs, offsets
from frb.galaxies import photom as frbphotom
from frb.surveys import survey_utils
from frb import utils
import pandas


            

def run(frb_input:pandas.core.series.Series, 
        lit_refs:str=None,
        override:bool=False, out_path:str=None,
        outfile:str=None):
    """Main method for generating a Host JSON file

    Args:
        frb_input (pandas.core.series.Series): Row of the CVS file
            providing the frb items
        lit_refs (str, optional): File of literature references. Defaults to None.
        override (bool, optional): Attempt to over-ride errors. 
            Mainly for time-outs of public data. Defaults to False.
        outfile (str, optional): Over-ride default outfile [not recommended; mainly for testing]
        out_path (str, optional): Over-ride default outfile [not recommended; mainly for testing]


    Raises:
        e: [description]
        ValueError: [description]
    """

    print("--------------------------------------")
    print(f"Building FRB JSON file for {frb_input.Name}")

    # Instantiate
    ifrb = FRB(frb_input.Name, 
              (frb_input.ra, frb_input.dec),
              frb_input.DM*units.pc/units.cm**3,
              z_frb=frb_input.z if np.isfinite(frb_input.z) else None, 
              repeater=frb_input.repeater)

    # DM_err
    if np.isfinite(frb_input['DM_err']):
        ifrb.DM_err = frb_input.DM_err * units.pc / units.cm**3

    # RM
    for key in ['RM', 'RM_err']:
        if np.isfinite(frb_input[key]):
            setattr(ifrb, key, frb_input[key] * units.rad / units.m**2)

    # Fluence
    for key in ['fluence', 'fluence_err']:
        if np.isfinite(frb_input[key]):
            setattr(ifrb, key, frb_input[key] * units.Jy * units.ms)

    # Error ellipse
    ifrb.set_ee(a=frb_input.ee_a, b=frb_input.ee_b, 
                theta=frb_input.ee_theta, cl=68.) 
    if np.isfinite(frb_input.ee_a_sys):
        ifrb.set_ee(a=frb_input.ee_a_sys, b=frb_input.ee_b_sys, 
                theta=frb_input.ee_theta, cl=68., stat=False) 

    # Add DM_ISM from NE2001
    ifrb.set_DMISM()

    # Refs
    ifrb.refs = frb_input.refs.split(',')

    # Pulses
    path = os.path.join(resource_filename('frb', 'data'), 'FRBs')
    tbl_file = os.path.join(path, 'FRB_pulses.csv')
    frb_pulses = pandas.read_csv(tbl_file)

    idx = np.where(frb_pulses.Name == frb_input.Name)[0]
    if len(idx) == 1:
        frb_pulse = frb_pulses.iloc[idx[0]]
        # Pulse properties
        pulse_dict = {}
        # Width and scattering
        for key in ['Wi', 'Wi_err', 'tscatt', 'tscatt_err']:
            if np.isfinite(frb_pulse[key]):
                pulse_dict[key] = frb_pulse[key] * units.ms
        ifrb.set_pulse(frb_pulse.freq*units.GHz, **pulse_dict)
        # References
        prefs = frb_pulse.refs.split(',')
        for pref in prefs:
            if pref not in ifrb.refs:
                ifrb.refs.append(pref)

    # Write
    if out_path is None:
        out_path = os.path.join(resource_filename('frb', 'data'), 'FRBs')
    ifrb.write_to_json(path=out_path) #'/home/xavier/Projects/FRB_Software/FRB/frb/tests/files')


def main(frbs:list, options:str=None, data_file:str=None, lit_refs:str=None,
         override:bool=False, outfile:str=None, out_path:str=None):
    """ Driver of the analysis

    Args:
        frbs (list): [description]
        options (str, optional): [description]. Defaults to None.
        data_file (str, optional): Alternate table than default
            for building FRBs. Defaults to None.
        lit_refs (str, optional): [description]. Defaults to None.
        override (bool, optional): [description]. Defaults to False.
        outfile (str, optional): [description]. Defaults to None.
            Here for testing
        out_path (str, optional): [description]. Defaults to None.
            Here for testing
    """
    '''
    # Options
    build_cigale, build_ppxf = False, False
    if options is not None:
        if 'cigale' in options:
            build_cigale = True
        if 'ppxf' in options:
            build_ppxf = True
    '''

    # Read public FRB table
    frb_tbl = load_frb_data(tbl_file=data_file)

    # Loop me
    if frbs[0] == 'all':
        frbs = frb_tbl.Name.values
    elif isinstance(frbs, list):
        pass

    for frb in frbs:
        # Grab the name
        frb_name = utils.parse_frb_name(frb, prefix='FRB')
        mt_idx = frb_tbl.Name == frb_name
        idx = np.where(mt_idx)[0].tolist()
        # Do it!
        for ii in idx:
            run(frb_tbl.iloc[ii], 
                lit_refs=lit_refs, override=override,
                outfile=outfile, out_path=out_path)

    # 
    print("All done!")

# Run em all
#  frb_build FRBs --frb 20121102,20171020,20180301,20180916,20180924,20181112