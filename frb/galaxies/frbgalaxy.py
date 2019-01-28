""" Module for galaxies related to FRBs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

import warnings

from pkg_resources import resource_filename

from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.special import hyp2f1
from scipy.interpolate import interp1d

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.cosmology import Planck15
from astropy.cosmology import z_at_value
from astropy import constants
from astropy.table import Table

from linetools import utils as ltu

from frb.galaxies import defs


class FRBGalaxy(object):
    """

    """
    def __init__(self, ra, dec, frb, cosmo=None):

        # Init
        self.coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        self.frb = frb

        # Cosmology
        if cosmo is None:
            self.cosmo = Planck15
        else:
            self.cosmo = cosmo

        # Main attributes
        self.redshift = {}
        self.photom = {}
        self.morphology = {}
        self.lines = {}
        self.derived = {}
        self.kinematics = {}

    def parse_photom(self, phot_tbl, max_off=1*units.arcsec, overwrite=True):
        phot_coord = SkyCoord(ra=phot_tbl['ra'], dec=phot_tbl['dec'], unit='deg')
        sep = self.coord.separation(phot_coord)
        row = np.argmin(sep)
        # Satisfy minimum offset?
        if sep[row] > max_off:
            print("No photometric sources within {} of the galaxy".format(max_off))
            return
        # Fill
        for filter in defs.valid_filters:
            # Value
            if filter in phot_tbl.keys():
                if (filter in self.photom.keys()) and (not overwrite):
                    pass
                else:
                    self.photom[filter] = phot_tbl[filter][row]
                    # Try error
                    if filter+'_err' in phot_tbl.keys():
                        self.photom[filter+'_err'] = phot_tbl[filter+'_err'][row]

    def parse_cigale(self, cigale_file, overwrite=True):
        # Read
        cigale_tbl = Table.read(cigale_file)

        # Derived quantities
        cigale_translate = [ # Internal key,  CIGALE key
            ('Mstar', 'bayes.stellar.m_star'),
            ('f_AGN', 'bayes.agn.fracAGN_dale2014'),
            ('u-r', 'bayes.param.restframe_u_prime-r_prime'),
            ('Lnu_r', 'bayes.param.restframe_Lnu(r_prime)'),
            ('SFR_photom', 'bayes.sfh.sfr'),
            ('EBV_photom', 'bayes.attenuation.E_BVs.stellar.old'),
            ('Z_photom', 'bayes.stellar.metallicity')
        ]
        # Do it
        cigale = {}
        for item in cigale_translate:
            if not(item[0] in defs.valid_derived_photom):
                raise AssertionError("{} not valid!!".format(item[0]))
            if item[1] in cigale_tbl.keys():
                cigale[item[0]] = cigale_tbl[item[1]][0]          # Solar masses, linear
                cigale[item[0]+'_err'] = cigale_tbl[item[1]+'_err'][0]          # Solar masses, linear

        # Absolute magnitude (astronomers...)
        if 'Lnu_r' in cigale.keys():
            cigale['M_r'] = -2.5*np.log10(cigale['Lnu_r']) + 34.1
            cigale['M_r_err'] = 2.5*(cigale['Lnu_r_err']/cigale['Lnu_r']) / np.log(10)

        # Fill Derived
        for key, item in cigale.items():
            if (key not in self.derived.keys()) or (overwrite):
                self.derived[key] = item

    def vette_dict(self, dict, valid_defs):
        pass

    def write_to_json(self, outfile=None):
        if outfile is None:
            jname = ltu.name_from_coord(self.coord)
            outfile = jname+'_FRB{}.json'.format(self.frb)
        # Build the dict
        frbgal_dict = {}

        # Basics
        frbgal_dict['ra'] = self.coord.ra.value
        frbgal_dict['dec'] = self.coord.dec.value
        frbgal_dict['FRB'] = self.frb
        frbgal_dict['cosmo'] = self.cosmo.name

        # Photometry
        if len(self.photom) > 0:
            frbgal_dict['photom'] = self.photom

        # Derived quantities
        if len(self.derived) > 0:
            frbgal_dict['derived'] = self.derived

        # JSONify
        jdict = ltu.jsonify(frbgal_dict)

        # Write
        ltu.savejson(outfile, jdict, easy_to_read=True, overwrite=True)
        print("Wrote data to {}".format(outfile))

    def __repr__(self):
        txt = '<{:s}: {:s} {:s}, FRB={:s}'.format(
            self.__class__.__name__, self.coord.icrs.ra.to_string(unit=units.hour,sep=':', pad=True),
            self.coord.icrs.dec.to_string(sep=':',pad=True,alwayssign=True), self.frb)
        # Finish
        txt = txt + '>'
        return (txt)


class FRBHost(FRBGalaxy):

    def __init__(self, ra, dec, frb, z_frb=None, **kwargs):
        # Instantiate
        super(FRBHost, self).__init__(ra, dec, frb, **kwargs)

        # Load up FRB info from name

        # Optional
        if z_frb is not None:
            self.redshift['z_FRB'] = z_frb
