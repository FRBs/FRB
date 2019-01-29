""" Module for galaxies related to FRBs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import pdb

from pkg_resources import resource_filename

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.cosmology import Planck15
from astropy import constants
from astropy.table import Table

from linetools import utils as ltu

from frb.galaxies import defs
from frb.galaxies import nebular


class FRBGalaxy(object):
    """

    """
    @classmethod
    def from_dict(cls, idict, **kwargs):
        # Instantiate
        slf = cls(idict['ra'], idict['dec'], idict['FRB'], **kwargs)

        if slf.cosmo.name != idict['cosmo']:
            raise AssertionError("Your cosmology does not match the expected.  Gotta deal..")

        # Fill me up
        for attr in slf.main_attr:
            if attr in idict.keys():
                setattr(slf,attr,idict[attr])
        #
        return slf

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
        self.neb_lines = {}
        self.kinematics = {}
        self.derived = {}
        self.main_attr = ('photom', 'redshift', 'morphology', 'neb_lines', 'kinematics', 'derived')

    @property
    def z(self):
        if len(self.redshift) == 0:
            return None
        else:
            return self.redshift['z']

    def calc_nebular_AV(self, method='Ha/Hb', **kwargs):
        # Checks
        assert len(self.neb_lines) > 0
        # Do it
        AV = nebular.calc_dust_extinct(self.neb_lines, method, **kwargs)
        self.derived['AV_nebular'] = AV

    def calc_nebular_SFR(self, method='Ha', **kwargs):
        # Checks
        assert len(self.neb_lines) > 0
        assert len(self.redshift) > 0
        # Dust?
        if 'AV_nebular' in self.derived.keys():
            AV = self.derived['AV_nebular']
            print("Using AV={} for a dust correction of the SFR".format(AV))
        else:
            print("Not making a dust correction of the SFR.  Set AV_nebular to do so or input AV to this method")
            AV = None
        # Calculate
        SFR = nebular.calc_SFR(self.neb_lines, method, self.redshift['z'], self.cosmo, AV=AV)
        self.derived['SFR_nebular'] = SFR.to('Msun/yr').value

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
    
    def gen_cigale_data_in(self, ID=None,filename='data.fits',overwrite=False):
        """
        Generates the input data file for CIGALE
        given the photometric points and redshift
        of a galaxy
        Parameters
        ----------
        self: FRBGalaxy or FRBHost object
            A Galaxy object having redshift estimates
            and photometric data available.
        ID: str, optional
            An ID for the galaxy. If none, "GalaxyA" is assigned.
        filename: str, optional
            Name of fits file (with path if needed) to store data in.
            Default value is 'data.fits'
        overwrite = bool, optional
            If true, previously written fits files will be
            overwritten
        """
        assert (self.photom != {}),"No photometry found. CIGALE cannot be run."
        assert (self.redshift != {}),"No redshift found. CIGALE cannot be run"

        photom = self.photom
        if ID is None:
            ID = "GalaxyA"
        photom['id'] = ID
        photom['redshift'] = self.z
        
        #Convert DES fluxes to mJy
        for band in defs.DES_bands:
            colname = "DES_"+band
            photom[colname] = 3630780.5*10**(photom[colname]/-2.5)
            photom[colname+"_err"] = photom[colname+"_err"]/1.087*photom[colname]
        
        #Convert WISE fluxes to mJy
        wise_fnu0 = [309.54,171.787,31.674,8.363] #http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2flux
        for band,zpt in zip(defs.WISE_bands,wise_fnu0):
            photom[band] = zpt*10**(-photom[band]/2.5)
            errname = band+"_err"
            if photom[errname]!=-999.0:
                photom[errname] =-99.0
            else:
                photom[errname] = photom[errname]/1.087*photom[band]
        
        #Write to file
        photom = Table([photom])
        photom.write(filename,format="fits",overwrite=overwrite)

    def parse_cigale(self, cigale_file, overwrite=True):
        # Read
        cigale_tbl = Table.read(cigale_file)

        # Derived quantities
        cigale_translate = [  # Internal key,  CIGALE key
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

    def parse_galfit(self, galfit_file, plate_scale, overwrite=True):
        # Will only allow for 1 sersic component for now
        lines = [line.rstrip('\n') for line in open(galfit_file)]

        galfit = {}
        # Search for Sersic
        for kk,line in enumerate(lines):
            if ('sersic' in line) and ('Component type' in line):
                for newline in lines[kk:]:
                    if 'R_e' in newline:
                        galfit['reff_ang'] = float(newline[4:11])*plate_scale
                    elif 'Sersic index' in newline:
                        galfit['n'] = float(newline[4:11])
                    elif 'Axis ratio' in newline:
                        galfit['b/a'] = float(newline[4:11])
                    elif 'Position angle' in newline:
                        galfit['PA'] = float(newline[4:11])
                    elif 'Component number' in newline:
                        break
        # Fill morphology
        for key, item in galfit.items():
            if (key not in self.morphology.keys()) or (overwrite):
                self.morphology[key] = item
        # reff kpc?
        if (self.z is not None) and ('reff_ang' in self.morphology.keys()):
            self.morphology['reff_kpc'] = \
                (self.morphology['reff_ang']*units.arcsec * self.cosmo.kpc_proper_per_arcmin(self.z)).to('kpc').value


    def parse_ppxf(self, ppxf_line_file, overwrite=True, format='ascii.ecsv'):

        ppxf_tbl = Table.read(ppxf_line_file, format=format)
        names = ppxf_tbl['name'].data
        ppxf_translate = [  # Internal key,  CIGALE key
            ('Ha', 'Halpha'),
            ('Hb', 'Hbeta'),
            ('Hg', 'Hgamma'),
            ('[NII] 6583',  '[NII]6583_d'),  # [NII] 6583 flux erg/s/cm^2; pPXF
            ('[OII] 3726',  '[OII]3726'),    # [OII] flux erg/s/cm^2; pPXF
            ('[OII] 3729',  '[OII]3729'),    # [OII] flux erg/s/cm^2; pPXF
            ('[OIII] 5007',  '[OII]5007_d')  # [OII] 5007 flux erg/s/cm^2; pPXF
        ]

        # Fluxes first
        ppxf = {}
        for item in ppxf_translate:
            if not(item[0] in defs.valid_neb_lines):
                raise AssertionError("{} not valid!!".format(item[0]))
            tidx = np.where(names == item[1])[0]
            if len(tidx) > 0:
                tidx = tidx[0]
                ppxf[item[0]] = ppxf_tbl['flux'][tidx]
                ppxf[item[0]+'_err'] = ppxf_tbl['err'][tidx]

        # Fill Nebular
        for line in defs.valid_neb_lines:
            # Value
            if line in ppxf.keys():
                if (line in self.neb_lines.keys()) and (not overwrite):
                    pass
                else:
                    self.neb_lines[line] = ppxf[line]
                    # Try error
                    if line+'_err' in ppxf.keys():
                        self.neb_lines[line+'_err'] = ppxf[line+'_err']

    def set_z(self, z, origin, err=None):
        # Set internal
        if origin == 'spec':
            key = 'z_spec'
        elif origin == 'phot':
            key = 'z_phot'
        self.redshift[key] = z
        if err is not None:
            self.redshift[key+'_err'] = err

        # Preferred?
        if (origin == 'spec') or ((origin == 'phot') and ('z' not in self.redshift.keys())):
            self.redshift['z'] = z
            if err is not None:
                self.redshift['z_err'] = err

    def vette_dict(self, dict, valid_defs):
        pass

    def make_outfile(self):
        jname = ltu.name_from_coord(self.coord)
        outfile = jname+'_FRB{}.json'.format(self.frb)
        return outfile

    def write_to_json(self, outfile=None, path='./'):
        if outfile is None:
            outfile = self.make_outfile()
        # Build the dict
        frbgal_dict = {}

        # Basics
        frbgal_dict['ra'] = self.coord.ra.value
        frbgal_dict['dec'] = self.coord.dec.value
        frbgal_dict['FRB'] = self.frb
        frbgal_dict['cosmo'] = self.cosmo.name

        # Main attributes
        for attr in self.main_attr:
            if len(getattr(self,attr)) > 0:
                frbgal_dict[attr] = getattr(self,attr)

        # JSONify
        jdict = ltu.jsonify(frbgal_dict)

        # Write
        ltu.savejson(os.path.join(path,outfile), jdict, easy_to_read=True, overwrite=True)
        print("Wrote data to {}".format(os.path.join(path,outfile)))

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

    def make_outfile(self):
        outfile = 'FRB{}_host.json'.format(self.frb)
        return outfile
