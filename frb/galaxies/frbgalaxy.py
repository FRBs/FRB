""" Module for galaxies related to FRBs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
from IPython import embed
import warnings

from pkg_resources import resource_filename

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.cosmology import Planck15
from astropy import constants
from astropy.table import Table

from frb.galaxies import defs
from frb.galaxies import nebular
from frb.galaxies import utils as gutils
from frb import utils



class FRBGalaxy(object):
    """
    Parent class for galaxies in FRB fields

    Simple object to hold key observable and derived quantities

    Warning:  Generating hundreds of these objects will likely be slow.
    Especially SkyCoord generation.  A new class will be warranted for that

    Args:
        ra (float): RA in deg
        dec (float): DEC in deg
        frb (str): Nomiker of the FRB, e.g. 121102
        cosmo (astropy.cosmology): Cosmology, e.g. Planck15

    Attributes:
        redshift (dict):
        photom (dict):
        morphology (dict):
        neb_lines (dict):
        kinematics (dict):
        derived (dict):

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
        slf = cls(idict['ra'], idict['dec'], idict['FRB'], **kwargs)

        # FRB coord
        if 'ra_FRB' in idict.keys():
            slf.frb_coord = SkyCoord(ra=idict['ra_FRB'], dec=idict['dec_FRB'], unit='deg')

        # Check cosmology
        if slf.cosmo.name != idict['cosmo']:
            raise AssertionError("Your cosmology does not match the expected.  Gotta deal..")

        # Fill me up
        for attr in slf.main_attr:
            if attr in idict.keys():
                setattr(slf,attr,idict[attr])
        # Return
        return slf

    @classmethod
    def from_json(cls, json_file, **kwargs):
        idict = utils.loadjson(json_file)
        slf = cls.from_dict(idict, **kwargs)
        return slf

    def __init__(self, ra, dec, frb, cosmo=None):

        # Init
        self.coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        self.frb = frb # Name, not coord

        self.frb_coord = None
        #
        self.name = ''

        # Cosmology
        if cosmo is None:
            self.cosmo = Planck15
        else:
            self.cosmo = cosmo

        # Main attributes
        self.eellipse = {}  # Error ellipse
        self.redshift = {}
        self.photom = {}
        self.morphology = {}
        self.neb_lines = {}
        self.kinematics = {}
        self.derived = {}
        self.main_attr = ('eellipse', 'photom', 'redshift', 'morphology', 'neb_lines', 'kinematics', 'derived')

    @property
    def z(self):
        """
        Return the redshift of the galaxy

        Returns:
            float or None: redshift or nadda

        """
        if len(self.redshift) == 0:
            return None
        else:
            return self.redshift['z']
    @property
    def z_err(self):
        """
        Return the redshift error of the galaxy

        Returns:
            float or None: redshift or nadda

        """
        if len(self.redshift) == 0:
            return None
        else:
            return self.redshift['z_err']

    def calc_nebular_lum(self, line):
        """
        Calculate the line luminosity
        Applies dust extinction if self.derived['AV_nebular'] is filled

        Mainly a wrapper to nebular.calc_lum()

        Args:
            line (str):  Name of the line
        """
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

        Lum, Lum_err = nebular.calc_lum(self.neb_lines, line, self.z, self.cosmo, AV=AV)
        return Lum, Lum_err

    def calc_nebular_AV(self, method='Ha/Hb', min_AV=None, **kwargs):
        """
        Calculate an A_V extinction from a pair of Nebular lines

        Mainly a wrapper to nebular.calc_dust_extinct

        self.derived['AV_nebular'] is filled

        Args:
            method (str): Method to use
            min_AV (float): Minimum A_V value allowed;  might set 0. someday
            **kwargs: Passed to nebular.calc_dust_extinct

        Returns:

        """
        # Checks
        assert len(self.neb_lines) > 0
        # Do it
        AV = nebular.calc_dust_extinct(self.neb_lines, method, **kwargs)
        if min_AV is not None:
            AV = max(AV, min_AV)
        # Set
        self.derived['AV_nebular'] = AV

    def calc_nebular_SFR(self, method='Ha', **kwargs):
        """
        Calculate a SFR from a nebular line

        Mainly a wrapper to nebular.calc_SFR

        self.derived['AV_nebular'] is filled with units SFR/yr

        Args:
            method (str):  Method to use, e.g. 'Ha' for Halpha
            **kwargs: passed to nebular.calc_SFR

        Returns:

        """
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
        """
        Parse photometry from an input table

        Fills the self.photom dict

        Args:
            phot_tbl (astropy.table.Table):
            max_off (Angle, optional):
            overwrite (bool, optional):

        Returns:

        """
        try:
            phot_coord = SkyCoord(ra=phot_tbl['ra'], dec=phot_tbl['dec'], unit='deg')
        except:
            embed(header='233 of frbgalaxy')
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
    
    def gen_cigale_data_in(self, ID=None, filename='data.fits', overwrite=False):
        """
        Generates the input data file for CIGALE
        given the photometric points and redshift
        of a galaxy

        Args:
            ID: str, optional
                An ID for the galaxy. If none, "GalaxyA" is assigned.
            filename: str, optional
                Name of fits file (with path if needed) to store data in.
                Default value is 'data.fits'
            overwrite = bool, optional
                If true, previously written fits files will be
                overwritten
        """
        assert (len(self.photom) > 0 ),"No photometry found. CIGALE cannot be run."
        assert (len(self.redshift) > 0),"No redshift found. CIGALE cannot be run"
        new_photom = Table([self.photom])
        if ID is None:
            ID = "GalaxyA"
        new_photom['id'] = ID
        new_photom['redshift'] = self.z
        
        # Convert DES fluxes to mJy
        for band in defs.DES_bands:
            colname = "DES_"+band
            new_photom[colname] = 3630780.5*10**(new_photom[colname]/-2.5)
            new_photom[colname+"_err"] = new_photom[colname]*(10**(new_photom[colname+"_err"]/2.5)-1)
        
        # Convert WISE fluxes to mJy
        wise_fnu0 = [309.54,171.787,31.674,8.363] #http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2flux
        for band,zpt in zip(defs.WISE_bands,wise_fnu0):
            # Data exists??
            if band not in self.photom.keys():
                continue
            #
            new_photom[band] = zpt*10**(-new_photom[band]/2.5)
            errname = band+"_err"
            if new_photom[errname] != -999.0:
                new_photom[errname] = -99.0
            else:
                new_photom[errname] = new_photom[errname]/1.087*new_photom[band]
        
        # Write to file
        try:
            new_photom.write(filename, format="fits", overwrite=overwrite)
        except OSError:
            warnings.warn("File exists;  use overwrite=True if you wish")

    def get_metaspec(self, instr=None, return_all=False, specdb_file=None):

        specDB = gutils.load_specdb(specdb_file=specdb_file)
        if specDB is None:
            return

        # Grab the spectra
        xspec, meta = specDB.spectra_from_coord(self.coord)

        # Return all?
        if return_all:
            return meta, xspec

        # Cut down
        if instr is None:
            if len(meta) > 1:
                warnings.warn("Multiple spectra returned for this galaxy.  Taking the first, but you may wish to specify your instrument")
                xspec = xspec[0]
                meta = meta[0]
        else:
            idx = meta['GROUP'] == instr
            if np.sum(idx) == 0:
                warnings.warn("No spectrum with instrument = {}".format(instr))
                return
            elif np.sum(idx) > 1:
                warnings.warn("Multiple spectra returned for this galaxy.  Taking the first, but you may wish to specify your instrument")
            xspec = xspec[np.where(idx)[0][0]]
            meta = meta[np.where(idx)[0][0]]
        # Return
        return meta, xspec


    def parse_cigale(self, cigale_file, overwrite=True):
        """
        Parse the output file from CIGALE

        Read into self.derived

        Args:
            cigale_file (str): Name of the CIGALE file
            overwrite (bool, optional):  Over-write any previous values

        Returns:

        """
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
        """
        Parse an output GALFIT file

        Loaded into self.morphology

        Args:
            galfit_file (str):
            plate_scale (float):  Plate scale in arcsec/pixel
            overwrite (bool, optional):

        Returns:

        """
        # Will only allow for 1 sersic component for now
        lines = [line.rstrip('\n') for line in open(galfit_file)]

        galfit = {}
        # Search for Sersic
        for kk,line in enumerate(lines):
            if line[1:7] == 'sersic':
                # Values
                prss = line.split(' ')
                keepp = [obj for obj in prss if obj != '']  # Remove white spaces
                galfit['PA'] = float(keepp[-1])
                galfit['b/a'] = float(keepp[-2])
                galfit['n'] = float(keepp[-3])
                galfit['reff_ang'] = float(keepp[-4])*plate_scale
                # Error
                prss = lines[kk+1].split(' ')
                keepp = [obj for obj in prss if obj != '']  # Remove white spaces
                galfit['PA_err'] = float(keepp[-1])
                galfit['b/a_err'] = float(keepp[-2])
                galfit['n_err'] = float(keepp[-3])
                galfit['reff_ang_err'] = float(keepp[-4])*plate_scale
        # Fill morphology
        for key, item in galfit.items():
            if (key not in self.morphology.keys()) or (overwrite):
                self.morphology[key] = item
        # reff kpc?
        if (self.z is not None) and ('reff_ang' in self.morphology.keys()):
            self.morphology['reff_kpc'] = \
                (self.morphology['reff_ang']*units.arcsec * self.cosmo.kpc_proper_per_arcmin(self.z)).to('kpc').value
            self.morphology['reff_kpc_err'] = \
                (self.morphology['reff_ang_err']*units.arcsec * self.cosmo.kpc_proper_per_arcmin(self.z)).to('kpc').value

    def parse_ppxf(self, ppxf_file, overwrite=True, format='ascii.ecsv'):
        """
        Parse an output pPXF file generated by our custom run

        Loaded into self.lines

        Args:
            ppxf_file (str): pPXF results file
            overwrite (bool, optional):
            format (str, optional):  Format of the table

        Returns:

        """

        ppxf_tbl = Table.read(ppxf_file, format=format)
        names = ppxf_tbl['name'].data
        ppxf_translate = [  # Internal key,  CIGALE key
            ('Halpha', 'Halpha'),
            ('Hbeta', 'Hbeta'),
            ('Hgamma', 'Hgamma'),
            ('[NII] 6584',  '[NII]6583_d'),  # [NII] 6583 flux erg/s/cm^2; pPXF
            ('[OII] 3726',  '[OII]3726'),    # [OII] flux erg/s/cm^2; pPXF
            ('[OII] 3729',  '[OII]3729'),    # [OII] flux erg/s/cm^2; pPXF
            ('[OIII] 5007',  '[OIII]5007_d')  # [OII] 5007 flux erg/s/cm^2; pPXF
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
        
        # Fitted quantities
        self.derived['EBV_spec'] = ppxf_tbl.meta['EBV']
        self.derived['Z_spec'] = ppxf_tbl.meta['METALS']
        self.derived['Mstar_spec'] = 10.**ppxf_tbl.meta['LOGMSTAR']


    def set_z(self, z, origin, err=None):
        """
        Set the redshift value(s) in self.redshift

        Args:
            z (float): Redshift value
            origin (str):  Origin
               'spec' for spectroscopic
               'phot' for photometric
            err (float, optional): Error in the redshift

        Returns:

        """
        # Set internal
        if origin == 'spec':
            key = 'z_spec'
        elif origin == 'phot':
            key = 'z_phot'
        else:
            raise IOError("Bad option for origin")
        #
        self.redshift[key] = z
        if err is not None:
            self.redshift[key+'_err'] = err

        # Preferred?
        if (origin == 'spec') or ((origin == 'phot') and ('z' not in self.redshift.keys())):
            self.redshift['z'] = z
            if err is not None:
                self.redshift['z_err'] = err

    def vet_one(self, attr):
        """
        Vette one of the main_attr

        Returns:
            bool: True = passed

        """
        vet = True
        # Check
        assert attr in self.main_attr

        # Setup
        if attr == 'neb_lines':
            defs_list = defs.valid_neb_lines
        elif attr == 'morphology':
            defs_list = defs.valid_morphology
        elif attr == 'photom':
            defs_list = defs.valid_photom
        elif attr == 'derived':
            defs_list = defs.valid_derived
        elif attr == 'redshift':
            defs_list = defs.valid_z
        elif attr == 'eellipse':
            defs_list = defs.valid_e
        else:
            return True
        # Vet
        for key in getattr(self, attr).keys():
            # Skip error
            if '_err' in key:
                continue
            if key not in defs_list:
                vet = False
                warnings.warn("{} in {} is not valid!".format(key,attr))
        # Return
        return vet

    def vet_all(self):
        """
        Vette all of the main dicts

        Args:
            dict:
            valid_defs:

        Returns:
            bool: True = passed

        """
        vet = True
        # Loop me
        for attr in self.main_attr:
            vet &= self.vet_one(attr)
        # Return
        return vet

    def make_outfile(self):
        """
        Auto-generate an output name for the class

        Returns:
            str: Output filename

        """
        jname = utils.name_from_coord(self.coord)
        outfile = jname+'_FRB{}.json'.format(self.frb)
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
        # Generate path as needed
        if not os.path.isdir(path):
            os.mkdir(path)
        if outfile is None:
            outfile = self.make_outfile()
        # Build the dict
        frbgal_dict = {}

        # Basics
        frbgal_dict['ra'] = self.coord.ra.value
        frbgal_dict['dec'] = self.coord.dec.value
        frbgal_dict['FRB'] = self.frb
        if self.frb_coord is not None:
            frbgal_dict['ra_FRB'] = self.frb_coord.ra.value
            frbgal_dict['dec_FRB'] = self.frb_coord.dec.value
        frbgal_dict['cosmo'] = self.cosmo.name

        # Main attributes
        for attr in self.main_attr:
            if len(getattr(self,attr)) > 0:
                frbgal_dict[attr] = getattr(self,attr)

        # JSONify
        jdict = utils.jsonify(frbgal_dict)

        # Write
        utils.savejson(os.path.join(path,outfile), jdict, easy_to_read=True, overwrite=overwrite)
        print("Wrote data to {}".format(os.path.join(path,outfile)))

    def __repr__(self):
        txt = '<{:s}: {:s} {:s}, FRB={:s}'.format(
            self.__class__.__name__, self.coord.icrs.ra.to_string(unit=units.hour,sep=':', pad=True),
            self.coord.icrs.dec.to_string(sep=':',pad=True,alwayssign=True), self.frb)
        # Finish
        txt = txt + '>'
        return (txt)


class FRBHost(FRBGalaxy):
    """
    Child of FRBGalaxy specific for an FRB host

    Args:
        ra (float): RA in deg
        dec (float): DEC in deg
        FRB (str): Nomiker of the FRB, e.g. 121102
        z_frb (float, optional):  Redshift of the host, expected to be provided

    """
    @classmethod
    def by_name(cls, frb, **kwargs):
        """
        
        Args:
            frb (str):  FRB name *without* FRB, e.g. 180924, not FRB180924
            **kwargs: 

        Returns:

        """
        path = os.path.join(resource_filename('frb', 'data/Galaxies/'), frb)
        json_file = os.path.join(path, FRBHost._make_outfile(frb))
        slf = cls.from_json(json_file, **kwargs)
        return slf

    def __init__(self, ra, dec, frb, z_frb=None, **kwargs):
        # Instantiate
        super(FRBHost, self).__init__(ra, dec, frb, **kwargs)

        # Load up FRB info from name
        self.name = 'HG{}'.format(self.frb)

        # Optional
        if z_frb is not None:
            self.redshift['z_FRB'] = z_frb

    @staticmethod
    def _make_outfile(frbname):
        """
        Static method to generate outfile based on frbname

        Args:
            frbname (str):  FRB name, e.g. 121102

        Returns:
            str: outfile

        """
        if frbname[0:3] != 'FRB':
            prefix = 'FRB'
        else:
            prefix = ''
        #
        outfile = '{}{}_host.json'.format(prefix, frbname)
        return outfile

    def make_outfile(self):
        """
        Overloads the parent method for Host specific naming

        Naming is FRBXXXXXX_host.json with XXXXXXX supplied by self.frb

        Returns:
            str:  Name of the default outfile

        """
        outfile = self._make_outfile(self.frb)
        return outfile

    def set_z(self, z, origin, err=None):
        """
        Partially overload the main method

        The main change is that the input z also sets z_FRB

        self.redshift is modified in place

        Args:
            z (float): Redshift value
            origin (str):  Origin
               'spec' for spectroscopic
               'phot' for photometric
            err (float, optional): Error in the redshift

        Returns:

        """
        super(FRBHost,self).set_z(z, origin, err=err)
        self.redshift['z_FRB'] = z
        if err is not None:
            self.redshift['z_FRB_err'] = err


class FGGalaxy(FRBGalaxy):

    def __init__(self, ra, dec, frb, **kwargs):
        # Instantiate
        super(FGGalaxy, self).__init__(ra, dec, frb, **kwargs)

        # Load up FRB info from name
        self.name = 'FG{}_{}'.format(self.frb, utils.name_from_coord(self.coord))



