""" Module for an FRB event
"""

from __future__ import print_function, absolute_import, division, unicode_literals

from pkg_resources import resource_filename
import os

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.cosmology import Planck15

from frb import utils


class generic_FRB(object):
    """
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
        idict = utils.loadjson(json_file)
        slf = cls.from_dict(idict, **kwargs)
        return slf

    def __init__(self, S, nu_c, DM, coord=None, cosmo=None):
        """
        Parameters
        ----------
        S : Quantity
          Source density of the burst
        width : Quantity
          Width
        nu_c : Quantity
          Centre frequency
        DM : Quantity
        coord : multi-format, optional
          RA/DEC in one of many formats (see utils.radec_to_coord)
        """
        self.S = S
        self.nu_c = nu_c
        self.DM = DM
        self.DM_err = None
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

        # dicts of attributes to be read/written
        self.main_dict = ['eellipse']

    def set_ee(self, a, b, theta, cl):
        """
        Set an error ellipse

        Args:
            a (float): major axis; Arcsec
            b (float):  minor axis; Arcsec
            theta (float): rotation of the major axis E from N (deg)
            cl (float): confidence level
        """
        self.eellipse['a'] = a
        self.eellipse['b'] = b
        self.eellipse['theta'] = theta
        self.eellipse['cl'] = cl
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

        # Measured properties
        for attr in ['S', 'nu_c', 'DM', 'DM_err', 'z']:
            if getattr(self,attr) is not None:
                frb_dict[attr] = getattr(self, attr)

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


class FRB(generic_FRB):

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
        if 'DM_err' in idict.keys():
            slf.DM_err = idict['DM_err']['value'] * units.pc / units.cm**3
            idict.pop('DM_err')

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
        path = os.path.join(resource_filename('frb', 'data/FRBs/'), frb)
        json_file = path + '.json'
        slf = cls.from_json(json_file, **kwargs)
        return slf

    def __init__(self, frb_name, coord, DM, S=None, nu_c=None, z_frb=None, **kwargs):
        # Instantiate
        super(FRB, self).__init__(S, nu_c, DM, coord=coord, **kwargs)

        self.frb_name = frb_name
        self.z = z_frb

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

