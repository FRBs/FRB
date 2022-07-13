""" Module for an FRB event
"""
import inspect

from pkg_resources import resource_filename
import os
import glob
import copy

import numpy as np

import pandas as pd

from astropy.coordinates import SkyCoord
from astropy import units

from linetools import utils as ltu

from frb import utils
from frb import mw
from frb import defs
from frb.galaxies import frbgalaxy

from IPython import embed


class GenericFRB(object):
    """
    Parent object for FRBs


    Args:
        S : Quantity
          Source density of the burst
        nu_c : Quantity
          Centre frequency
        DM : Quantity
        coord (astropy.coordinates.SkyCoord): multi-format, optional
          RA/DEC in one of many formats (see utils.radec_to_coord)
        cosmo:

    Attributes:
        fluence (Quantity):
            Fluence
        fluence_err (Quantity):
        DM (Quantity):
            Dispersion Measure
        DM_err (Quantity):
        RM (Quantity):
            Rotation Measure
        RM_err (Quantity):
        lpol (float):
            Linear Polarization (%)
        lpol_err (Quantity):
        refs (list):
            List of str, reference names
        z (float):
            Redshift
        z_err (float):
            Uncertainty in the redshift
        repeater (bool):
            Marks the FRB as being a Repeater

    """
    @classmethod
    def from_dict(cls, idict, **kwargs):
        """
        Instantiate from a dict

        Args:
            idict (dict):
            **kwargs: Passed to the __init__ call

        Returns:

        """
        # Init
        slf = cls(idict['S'], idict['nu_c'], idict['DM'], **kwargs)
        for key in ['S','nu_c','DM']:
            idict.pop(key)

        # FRB coord
        if 'ra' in idict.keys():
            slf.coord = SkyCoord(ra=idict['ra'], 
                                 dec=idict['dec'], 
                                 unit='deg')

        # Check cosmology
        if slf.cosmo.name != idict['cosmo']:
            raise AssertionError("Your cosmology does not match the expected.  Gotta deal..")

        # dicts
        for ndict in slf.main_dict:
            if ndict in idict.keys():
                setattr(slf,ndict,idict[ndict])
                idict.pop(ndict)

        # Remainder
        for key in idict.keys():
            setattr(slf,key,idict[key])

        # Return
        return slf

    @classmethod
    def from_json(cls, json_file, **kwargs):
        """
        Instantiate from a JSON file
          A simple wrapper to the from_dict method

        Args:
            json_file (str):
            **kwargs: Passed to from_dict()

        Returns:
            slf

        """
        idict = utils.loadjson(json_file)
        slf = cls.from_dict(idict, **kwargs)
        return slf

    def __init__(self, S, nu_c, DM, coord=None, cosmo=None, repeater=None):
        """
        """
        self.S = S
        self.nu_c = nu_c
        # NE2001 (for speed)
        self.DMISM = None
        self.DMISM_err = None
        # Repeater?
        self.repeater = repeater
        # Coord
        if coord is not None:
            self.coord = utils.radec_to_coord(coord)
        else:
            self.coord = None
        # Cosmology
        if cosmo is None:
            self.cosmo = defs.frb_cosmo
        else:
            self.cosmo = cosmo

        # Attributes
        self.z = None
        self.frb_name = None

        self.fluence = None
        self.fluence_err = None
        self.DM = DM
        self.DM_err = None
        self.RM = None
        self.RM_err = None
        self.lpol = None
        self.lpol_err = None

        self.refs = []

        # dicts of attributes to be read/written
        self.eellipse = {}
        self.pulse = {}
        self.main_dict = ['eellipse', 'pulse']

    def set_DMISM(self):
        if self.coord is None:
            print("Need to set coord first!")
        self.DMISM = mw.ismDM(self.coord)

    def set_ee(self, a, b, theta, cl, stat=True):
        """
        Set an error ellipse for the FRB position

        Args:
            a (float): major axis; Arcsec
            b (float):  minor axis; Arcsec
            theta (float): rotation of the major axis E from N (deg)
            cl (float): confidence level
            stat (bool, optional):
                If True, fill in statistical error
                if False, fill in systematic
        """
        if a < b:
            raise IOError("For the ellipse, a must be greater than or equal to b")
        if stat:
            self.eellipse['a'] = a
            self.eellipse['b'] = b
            self.eellipse['theta'] = theta
            self.eellipse['cl'] = cl
        else:
            self.eellipse['a_sys'] = a
            self.eellipse['b_sys'] = b
            self.eellipse['theta_sys'] = theta
            self.eellipse['cl_sys'] = cl
        #
        return

    @property
    def sig_a(self):
        """
        Combined semi-major axis error

        Returns:
            float:

        """

        if len(self.eellipse) == 0:
            return None
        siga = self.eellipse['a']  # arcsec
        if 'a_sys' in self.eellipse.keys():
            siga = np.sqrt(self.eellipse['a_sys']**2 + siga**2)
        return siga

    @property
    def sig_b(self):
        """
        Combined semi-minor axis error

        Returns:
            float:

        """
        if len(self.eellipse) == 0:
            return None
        sigb = self.eellipse['b']  # arcsec
        if 'b_sys' in self.eellipse.keys():
            sigb = np.sqrt(self.eellipse['b_sys']**2 + sigb**2)
        return sigb

    def set_pulse(self, freq,
                  time_res=None, t0=None, Wi=None, Wi_err=None,
                  tscatt=None, tscatt_err=None, scatt_index=None,
                  scatt_index_err=None, DM_smear=None):
        """
        Args:
            freq (Quantity):
                Frequency at which the pulse was analyzed
            time_res (Quantity):
                Time resolution of the telescope/instrument
            t0 (Quantity):
                Pulse arrival time (MJD) at top band frequency
            Wi (Quantity):
                Intrinsic width
            Wi_err (Quantity):
                Error in intrinsic width
            tscatt (Quantity):
                Scattering broadening time
            tscatt_err (Quantity):
                Error in Scattering broadening time
            scatt_index (float):
                Scattering index
            scatt_index_err (float):
                Error in scattering index
            DM_smear (float):
                Dispersion smearing generated observed width
        """
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        self.pulse = dict([(k,values[k]) for k in args[1:]])

    def make_outfile(self):
        """
        Simple method for naming the output file

        Returns:
            str

        """
        if self.frb_name is None:
            outfile = 'Generic_FRB.json'
        else:
            outfile = '{:s}.json'.format(self.frb_name)
        #
        return outfile

    def write_to_json(self, outfile=None, path='./', overwrite=True):
        """
        Write key aspects of the class to disk in a JSON file

        Args:
            outfile (str, optional): Output filename
              If not provided, one will be generated with make_outfile()
            path (str, optional): Path for the output file
            overwrite (bool, optional): Overwrite?

        Returns:

        """
        if outfile is None:
            outfile = self.make_outfile()
        # Build the dict
        frb_dict = {}

        # Basics
        if self.coord is not None:
            frb_dict['ra'] = self.coord.ra.value
            frb_dict['dec'] = self.coord.dec.value
        if self.frb_name is not None:
            frb_dict['FRB'] = self.frb_name
        frb_dict['cosmo'] = self.cosmo.name
        frb_dict['refs'] = self.refs

        if self.repeater is not None:
            frb_dict['repeater'] = self.repeater

        # Measured properties
        for attr in ['S', 'nu_c', 'DM', 'z', 'RM', 'DMISM', 'fluence', 'lpol']:
            # Value
            if getattr(self,attr) is not None:
                frb_dict[attr] = getattr(self, attr)
            # Error
            if hasattr(self, attr+'_err'):
                if getattr(self, attr+'_err') is not None:
                    frb_dict[attr+'_err'] = getattr(self, attr+'_err')

        # Main dicts
        for idict in self.main_dict:
            if getattr(self,idict) is not None and len(getattr(self,idict)) > 0:
                frb_dict[idict] = getattr(self,idict)

        # JSONify
        jdict = utils.jsonify(copy.deepcopy(frb_dict))

        # Write
        utils.savejson(os.path.join(path,outfile), jdict, easy_to_read=True, overwrite=overwrite)
        print("Wrote data to {}".format(os.path.join(path,outfile)))

    def __repr__(self):
        txt = '<{:s}: S={} nu_c={}, DM={}'.format(
                self.__class__.__name__, self.S, self.nu_c, self.DM)
        # Finish
        txt = txt + '>'
        return (txt)


class FRB(GenericFRB):
    """
    FRB class used for actual, observed FRBs


    """

    @classmethod
    def from_dict(cls, idict, **kwargs):
        """
        Instantiate from a dict

        Args:
            idict (dict):
            **kwargs: Passed to the __init__ call

        Returns:

        """
        # Init
        coord = SkyCoord(ra=idict['ra'], dec=idict['dec'], unit='deg')
        DM = units.Quantity(idict['DM']['value'],unit=idict['DM']['unit'])

        slf = cls(idict['FRB'], coord, DM, **kwargs)
        for key in ['ra','dec','DM']:
            idict.pop(key)
        for key in ['DM_err', 'DMISM', 'DMISM_err', 'RM', 'RM_err', 'fluence', 'fluence_err']:
            if key in idict.keys():
                setattr(slf,key,units.Quantity(idict[key]['value'], unit=idict[key]['unit']))
                idict.pop(key)
        # Cosmology
        if slf.cosmo.name != idict['cosmo']:
            raise AssertionError(f"Your cosmology does not match the expected for {idict['FRB']}.  Gotta deal..")
        idict.pop('cosmo')

        # dicts
        for ndict in slf.main_dict:
            if ndict in idict.keys():
                for key, value in idict[ndict].items():
                    if isinstance(value, dict):
                        newvalue = ltu.convert_quantity_in_dict(value)
                    else:
                        newvalue = value
                    idict[ndict][key] = newvalue
                setattr(slf,ndict,idict[ndict])
                # Deal with quantities
                idict.pop(ndict)

        # Remainder
        for key in idict.keys():
            setattr(slf,key,idict[key])

        # Return
        return slf

    @classmethod
    def by_name(cls, frb_name, **kwargs):
        """
        Method to instantiate an FRB by its name

        Args:
            frb_name (str):
              Name of the FRB,
            **kwargs:

        Returns:

        """
        path = os.path.join(resource_filename('frb', 'data'), 'FRBs', frb_name)
        json_file = path + '.json'
        slf = cls.from_json(json_file, **kwargs)
        return slf

    def __init__(self, frb_name, coord, DM, S=None, nu_c=None, z_frb=None, **kwargs):
        """

        Args:
            frb_name (str):
            coord (astropy.coordinates.SkyCoord):
            DM (Quantity):
            S (Quantity):
                Source density
            nu_c:
            z_frb (float):
                Redshift
            **kwargs:
        """
        super(FRB, self).__init__(S, nu_c, DM, coord=coord, **kwargs)

        self.frb_name = frb_name
        self.z = z_frb

    def grab_host(self):
        """
        Returns the FRBHost object for this FRB

        Returns:
            frb.galaxies.frbgalaxy.FRBHost

        """
        frbHost = frbgalaxy.FRBHost.by_frb(self)
        return frbHost

    def __repr__(self):
        txt = '<{:s}: {} J{}{} DM={}'.format(
            self.__class__.__name__, self.frb_name,
            self.coord.icrs.ra.to_string(unit=units.hour, sep='', pad=True),
            self.coord.icrs.dec.to_string(sep='', pad=True, alwayssign=True),
            self.DM)
        if self.z is not None:
            txt += ' z={}'.format(self.z)
        # Finish
        txt = txt + '>'
        return (txt)


def list_of_frbs(require_z=False):
    """
    Generate a list of FRB objects for all the FRBs in the Repo

    Args:
        require_z (bool, optional):
            If True, require z be set

    Returns:
        list:

    """
    # Grab the files
    frb_files = glob.glob(os.path.join(resource_filename('frb', 'data'), 'FRBs', 'FRB*json'))
    frb_files.sort()
    # Load up the FRBs
    frbs = []
    for frb_file in frb_files:
        frb_name = os.path.basename(frb_file).split('.')[0]
        frb = FRB.by_name(frb_name)
        if require_z and frb.z is None:
            continue
        frbs.append(frb)
    # Return
    return frbs


def build_table_of_frbs(frbs=None, fattrs=None):
    """
    Generate a Pandas table of FRB data

    Warning:  As standard, missing values are given NaN in the Pandas table
        Be careful!

    Args:
        fattrs (list, optional):
            Float attributes for the Table
            The code also, by default, looks for accompanying _err attributes

    Returns:
        pd.DataFrame, dict:  Table of data on FRBs,  dict of their units

    """
    if fattrs is None:
        fattrs = ['DM', 'fluence', 'RM', 'lpol', 'z', 'DMISM']
    # Load up the FRBs
    if frbs is None:
        frbs = list_of_frbs()

    # Table
    frb_tbl = pd.DataFrame({'FRB': [ifrb.frb_name for ifrb in frbs]})
    tbl_units = {}
    tbl_units['FRB'] = None

    # Coordinates
    coords = SkyCoord([ifrb.coord for ifrb in frbs])
    frb_tbl['RA'] = coords.ra.value
    frb_tbl['DEC'] = coords.dec.value
    tbl_units['RA'] = 'deg'
    tbl_units['DEC'] = 'deg'

    # Error ellipses
    ee_attrs = ['a', 'b', 'a_sys', 'b_sys', 'theta']
    ee_units = ['arcsec', 'arcsec', 'arcsec', 'arcsec', 'deg']
    for ss, ee_attr in enumerate(ee_attrs):
        alist = [ifrb.eellipse[ee_attr] if ee_attr in ifrb.eellipse.keys() else np.nan for ifrb in frbs]
        frb_tbl['ee_'+ee_attr] = alist
        tbl_units['ee_'+ee_attr] = ee_units[ss]

    # Pulse
    pulse_attrs = ['Wi', 'tscatt']
    pulse_errors = [ipulse+'_err' for ipulse in pulse_attrs]
    pulse_error_units = ['ms']*len(pulse_errors)
    pulse_attrs += pulse_errors
    pulse_units = ['ms', 'ms'] + pulse_error_units
    for ss, pulse_attr in enumerate(pulse_attrs):
        alist = [ifrb.pulse[pulse_attr] if pulse_attr in ifrb.pulse.keys() else np.nan for ifrb in frbs]
        frb_tbl['pulse_'+pulse_attr] = alist
        tbl_units['pulse_'+pulse_attr] = pulse_units[ss]

    # A few others
    for other in ['repeater']:
        alist = [getattr(ifrb, other) if hasattr(ifrb, other) else np.nan for ifrb in frbs]
        frb_tbl[other] = alist

    # Refs
    alist = [','.join(ifrb.refs) for ifrb in frbs]
    frb_tbl['refs'] = alist

    # Float Attributes on an Object
    for fattr in fattrs:
        values = []
        # Error
        errors = []
        has_error = False
        # Now loop me
        for ss, ifrb in enumerate(frbs):
            if hasattr(ifrb, fattr) and getattr(ifrb, fattr) is not None:
                utils.assign_value(ifrb, fattr, values, tbl_units)
            else:
                values.append(np.nan)
            # Try error
            eattr = fattr+'_err'
            if hasattr(ifrb, eattr) and getattr(ifrb, eattr) is not None:
                has_error = True
                utils.assign_value(ifrb, eattr, errors, tbl_units)
            else:
                errors.append(np.nan)
        # Add to Table
        frb_tbl[fattr] = values
        if has_error:
            frb_tbl[eattr] = errors


    # Return
    return frb_tbl, tbl_units

def load_frb_data(tbl_file:str=None):
    if tbl_file is None:
        path = os.path.join(resource_filename('frb', 'data'), 'FRBs')
        tbl_file = os.path.join(path, 'FRBs_base.csv')

    # Load
    frb_tbl = pd.read_csv(tbl_file)

    # Return
    return frb_tbl
