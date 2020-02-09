""" Module for an FRB event
"""

from __future__ import print_function, absolute_import, division, unicode_literals

from pkg_resources import resource_filename
import os
import glob

import numpy as np

import pandas as pd

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.cosmology import Planck15

from frb import utils
from frb import mw
from frb.galaxies import frbgalaxy


class GenericFRB(object):
    """
    Parent object for FRBs


    Args:
        S : Quantity
          Source density of the burst
        nu_c : Quantity
          Centre frequency
        DM : Quantity
        coord : multi-format, optional
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
            slf.coord = SkyCoord(ra=idict['ra'], dec=idict['dec'], unit='deg')

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

    def __init__(self, S, nu_c, DM, coord=None, cosmo=None):
        """
        """
        self.S = S
        self.nu_c = nu_c
        # NE2001 (for speed)
        self.DMISM = None
        self.DMISM_err = None
        # Coord
        if coord is not None:
            self.coord = utils.radec_to_coord(coord)
        else:
            self.coord = None
        # Cosmology
        if cosmo is None:
            self.cosmo = Planck15
        else:
            self.cosmo = cosmo

        # Attributes
        self.eellipse = {}
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
        self.main_dict = ['eellipse']

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

    def set_width(self, wtype, value, overwrite=False):
        """ Set a Width value

        Parameters
        ----------
        wtype : str
          Indicates the type of width
        value : Quantity
        overwrite : bool, optional
        """
        # Restrict types
        assert wtype in ['Wi', 'Wb', 'Wobs']  # Intrinsic, broadened, observed
        # Check
        if hasattr(self, wtype) and (not overwrite):
            raise IOError("Width type {:s} is already set!".format(wtype))
        # Vette
        try:
            value.to('s')
        except:
            raise IOError("Bad Quantity for value")
        # Set
        setattr(self, wtype, value)

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
        jdict = utils.jsonify(frb_dict)

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
            raise AssertionError("Your cosmology does not match the expected.  Gotta deal..")
        idict.pop('cosmo')

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
    def by_name(cls, frb, **kwargs):
        """
        Method to instantiate an FRB by its name

        Args:
            frb (str):
              Name of the FRB,
            **kwargs:

        Returns:

        """
        path = os.path.join(resource_filename('frb', 'data/FRBs/'), frb)
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
        frbHost = frbgalaxy.FRBHost.by_name(self.frb_name[3:])
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


def build_table_of_frbs(fattrs=None):
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
        fattrs = ['DM', 'fluence', 'RM', 'lpol', 'z']
    # Grab the files
    frb_files = glob.glob(os.path.join(resource_filename('frb', 'data'), 'FRBs', 'FRB*json'))
    frb_files.sort()
    # Load up the FRBs
    frbs = []
    for frb_file in frb_files:
        frb_name = os.path.basename(frb_file).split('.')[0]
        frbs.append(FRB.by_name(frb_name))

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

