""" Module for galaxies related to FRBs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import warnings
import glob


from pkg_resources import resource_filename

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.cosmology import Planck15
from astropy.table import Table

from frb.galaxies import defs
from frb.galaxies import nebular
from frb.galaxies import utils as gutils
from frb.galaxies import offsets
from frb.surveys.catalog_utils import convert_mags_to_flux
from frb import utils

from scipy.integrate import simps

from IPython import embed

class FRBGalaxy(object):
    """
    Parent class for galaxies in FRB fields

    Simple object to hold key observable and derived quantities

    Warning:  Generating hundreds of these objects will likely be slow.
    Especially SkyCoord generation.  A new class will be warranted for that

    Args:
        ra (float): RA in deg
        dec (float): DEC in deg
        frb (frb.FRB): FRB object
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
    def from_dict(cls, frb, idict, **kwargs):
        """
        Instantiate from a dict

        Args:
            frb (frb.FRB):
            idict (dict):
            **kwargs: Passed to the __init__ call

        Returns:

        """
        # Init
        slf = cls(idict['ra'], idict['dec'], frb, **kwargs)

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

        # Physical Offset -- remove this when these get into JSON files
        if 'z_spec' in slf.redshift.keys() and 'physical' not in slf.offsets.keys():
            slf.offsets['physical'] = (slf.offsets['ang_best']*units.arcsec / slf.cosmo.arcsec_per_kpc_proper(slf.z)).to('kpc').value
            slf.offsets['physical_err'] = (slf.offsets['ang_best_err']*units.arcsec / slf.cosmo.arcsec_per_kpc_proper(slf.z)).to('kpc').value

        # Return
        return slf

    @classmethod
    def from_json(cls, frb, json_file, **kwargs):
        """

        Args:
            frb (frb.FRB):
            json_file:
            **kwargs:

        Returns:
            FRBGalaxy or None

        """
        try:
            idict = utils.loadjson(json_file)
        except FileNotFoundError:
            warnings.warn("File {} not found.  This galaxy probably does not exist yet.".format(json_file))
            return None
        slf = cls.from_dict(frb, idict, **kwargs)
        return slf

    def __init__(self, ra, dec, frb, cosmo=None):
        """

        Args:
            ra (float):
            dec (float):
            frb (frb.FRB) :
            cosmo (astropy.cosmology, optional):
        """

        # Init
        self.coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        self.frb = frb  # FRB object
        #
        self.name = ''

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
        self.offsets = {}
        self.positional_error = {}
        self.main_attr = ('photom', 'redshift', 'morphology', 'neb_lines',
                          'kinematics', 'derived', 'offsets', 'positional_error')

        # Angular offset
        self.offsets['ang_avg'], self.offsets['ang_avg_err'], \
                self.offsets['ang_best'], self.offsets['ang_best_err'] \
                    = offsets.angular_offset(frb, self)

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
    
    def calc_tot_uncert(self):
        """Calculate total uncertainty in arcsec of 
        FRB localization + Host localization in the 
        reference frame of the FRB

        Returns:
            tuple: uncerta, uncertb [arcsec]
        """
            # set to zero, but change if we have astrometric and source errors
        if hasattr(self, 'positional_error'):
            host_ra_sig = np.sqrt(self.positional_error['ra_astrometric']**2 +  
                self.positional_error['ra_source']**2)
            host_dec_sig = np.sqrt(self.positional_error['dec_astrometric']**2 + 
                self.positional_error['dec_source']**2)
        else:
            host_ra_sig, host_dec_sig = 0., 0.

        # Rotate to the FRB frame
        # sigma**2
        # will be zero if no positional errors saved in host json file
        theta = self.frb.eellipse['theta']
        sig2_gal_a = host_dec_sig ** 2 * np.cos(theta) ** 2 + host_ra_sig ** 2 * np.sin(theta) ** 2
        sig2_gal_b = host_ra_sig ** 2 * np.cos(theta) ** 2 + host_dec_sig ** 2 * np.sin(theta) ** 2

        # will only be FRB error if positional errors saved in host json file
        #  Units are pixels
        uncerta = np.sqrt(self.frb.sig_a**2 + sig2_gal_a)
        uncertb = np.sqrt(self.frb.sig_b**2 + sig2_gal_b) 

        # Return
        return uncerta, uncertb


    def parse_photom(self, phot_tbl, max_off=1*units.arcsec, overwrite=True, EBV=None):
        """
        Parse photometry from an input table

        Fills the self.photom dict

        Args:
            phot_tbl (astropy.table.Table):
            max_off (Angle, optional):
            overwrite (bool, optional):
            EBV (float, optional):  Galactic reddening.  If included, the photometry
               has been corrected for this.  If not, who knows?!  :)

        Returns:

        """
        phot_coord = SkyCoord(ra=phot_tbl['ra'], dec=phot_tbl['dec'], unit='deg')
        sep = self.coord.separation(phot_coord)
        row = np.argmin(sep)
        # Satisfy minimum offset?
        if sep[row] > max_off:
            print("No photometric sources within {} of the galaxy".format(max_off))
            return
        # Get a flux table
        flux_tbl = convert_mags_to_flux(phot_tbl, fluxunits='mJy')
        # Fill
        for filter in defs.valid_filters:
            # Value
            if filter in phot_tbl.keys():
                if (filter in self.photom.keys()) and (not overwrite):
                    pass
                else:
                    # -999. is used as empty fill value
                    if phot_tbl[filter][row] < -990:
                        pass
                    else:
                        self.photom[filter] = phot_tbl[filter][row]
                        self.photom[filter+'_flux'] = flux_tbl[filter][row]
                        # Try error
                        if filter+'_err' in phot_tbl.keys():
                            self.photom[filter+'_err'] = phot_tbl[filter+'_err'][row]
                            self.photom[filter+'_flux_err'] = flux_tbl[filter+'_err'][row]


                        # Add entries for corresponding flux values.
        # EBV
        if EBV is not None:
            self.photom['EBV'] = EBV
    
    def run_cigale(self, data_file="cigale_in.fits", config_file="pcigale.ini",
        wait_for_input=False, plot=True, outdir='out', compare_obs_model=False, **kwargs):
        """
        Generates the input data file for CIGALE
        given the photometric points and redshift
        of a galaxy

        Args:
            ID: str, optional
                An ID for the galaxy. If none, "GalaxyA" is assigned.
            data_file (str, optional):
                Root name for the photometry data file generated used as input to CIGALE
            config_file (str, optional):
                Root name for the file where CIGALE's configuration is generated
            wait_for_input (bool, optional):
                If true, waits for the user to finish editing the auto-generated config file
                before running.
            plot (bool, optional):
                Plots the best fit SED if true
            cores (int, optional):
                Number of CPU cores to be used. Defaults
                to all cores on the system.
            outdir (str, optional):
                Path to the many outputs of CIGALE
                If not supplied, the outputs will appear in a folder named out/
            compare_obs_model (bool, optional):
                If True compare the input observed fluxes with the model fluxes
                This writes a Table to outdir named 'photo_observed_model.dat'

        kwargs:  These are passed into gen_cigale_in() and _initialise()
            sed_modules (list of 'str', optional):
                A list of SED modules to be used in the 
                PDF analysis. If this is being input, there
                should be a corresponding correct dict
                for sed_modules_params.
            sed_module_params (dict, optional):
                A dict containing parameter values for
                the input SED modules. Better not use this
                unless you know exactly what you're doing.
        """
        # Adding import statement here in case CIGALE is
        # not installed.
        from .cigale import run

        assert (len(self.photom) > 0 ),"No photometry found. CIGALE cannot be run."
        assert (len(self.redshift) > 0),"No redshift found. CIGALE cannot be run"
        new_photom = Table([self.photom])
        new_photom['redshift'] = self.z
        if self.name != '':
            new_photom['ID'] = self.name
        else:
            new_photom['ID'] = 'FRBGalaxy'


        run(new_photom, 'redshift', data_file, config_file,
        wait_for_input, plot, outdir, compare_obs_model, **kwargs)
        return

    def get_metaspec(self, instr=None, return_all=False, specdb_file=None):
        """
        Return the meta data and spectra for this FRBGalaxy
        from the specDB

        If there is more than one spectrum, the code returns the first
        unless return_all=True

        Args:
            instr (str, optional):
                Restrict to the input Instrument
            return_all (bool, optional):
                Return all of the meta, spectra
            specdb_file (str, optional):
                Path+name of the specDB file to use (over-ride the default)

        Returns:
            astropy.table.Table, linetools.spectra.XSpectrum1D: meta data, spectra

        """

        specDB = gutils.load_specdb(specdb_file=specdb_file)
        if specDB is None:
            return

        # Grab the spectra
        xspec, meta = specDB.spectra_from_coord(self.coord)  # Tolerance is 0.5''

        # Return all?
        if return_all:
            return meta, xspec

        # Cut down
        if instr is None:
            if len(meta) > 1:
                warnings.warn("Multiple spectra returned for this galaxy.  Taking the first, but you may wish to specify your instrument")
                xspec = xspec[0]
                meta = meta[0:1]
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


    def parse_cigale(self, cigale_file, sfh_file=None, overwrite=True):
        """
        Parse the output file from CIGALE

        Read into self.derived

        Args:
            cigale_file (str): Name of the CIGALE results file
            sfh_file (str, optional): Name of the best SFH model file.
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

        # Compute mass weighted age?
        if sfh_file is not None:
            try:
                sfh_tab = Table.read(sfh_file)
            except:
                warnings.warn("Invalid SFH file. Skipping mass-weighted age.")
                return
            mass = simps(sfh_tab['SFR'], sfh_tab['time']) # M_sun/yr *Myr
            # Computed mass weighted age
            t_mass = simps(sfh_tab['SFR']*sfh_tab['time'], sfh_tab['time'])/mass # Myr
            # Store
            if ('age_mass' not in self.derived.keys()) or (overwrite):
                cigale['age_mass'] = t_mass

        # Fill Derived
        for key, item in cigale.items():
            if (key not in self.derived.keys()) or (overwrite):
                self.derived[key] = item
        

    def parse_galfit(self, galfit_file, overwrite=True, twocomponent=False):
        """
        Parse an output GALFIT file

        Loaded into self.morphology

        Args:
            galfit_file (str): processed 'out.fits' file
                produced by frb.galaxies.galfit.run. Contains
                a binary table with fit parameters.
            overwrite (bool, optional): Need to overwrite
                the object's attributes?
            twocomponent (bool, optional): Should the morphology
                attribute generated contain fit parameters of
                two components?

        """
        assert os.path.isfile(galfit_file), "Incorrect file path {:s}".format(galfit_file)
        try:
            fit_tab = Table.read(galfit_file, hdu=4)
        except:
            raise IndexError("The binary table with fit parameters was not found as the 4th hdu in {:s}. Was GALFIT run using the wrapper?".format(galfit_file))
        for key in fit_tab.keys():
            if 'mag' in key:
                continue
            if (key not in self.morphology.keys()) or (overwrite):
                if twocomponent:
                    self.morphology[key] = fit_tab[key].data
                else:
                    self.morphology[key] = fit_tab[key][0]
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

        # Physical offset
        if origin == 'spec':
            self.offsets['physical'] = (self.offsets['ang_best'] * units.arcsec *
                                        self.cosmo.kpc_proper_per_arcmin(self.z)).to('kpc').value
            self.offsets['physical_err'] = (self.offsets['ang_best_err'] * units.arcsec *
                                        self.cosmo.kpc_proper_per_arcmin(self.z)).to('kpc').value

    def vet_one(self, attr):
        """
        Vette one of the main_attr

        Parameters:
            attr (str):

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
            defs_list = defs.valid_photom+defs.valid_flux
        elif attr == 'derived':
            defs_list = defs.valid_derived
        elif attr == 'redshift':
            defs_list = defs.valid_z
        elif attr == 'offsets':
            defs_list = defs.valid_offsets
        elif attr == 'positional_error':
            defs_list = defs.valid_positional_error
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
        frbgal_dict['FRB'] = self.frb.frb_name
        if self.frb.coord is not None:
            frbgal_dict['ra_FRB'] = self.frb.coord.ra.value
            frbgal_dict['dec_FRB'] = self.frb.coord.dec.value
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
            self.coord.icrs.dec.to_string(sep=':',pad=True,alwayssign=True), self.frb.frb_name)
        # Finish
        txt = txt + '>'
        return (txt)


class FRBHost(FRBGalaxy):
    """
    Child of FRBGalaxy specific for an FRB host

    Args:
        ra (float): RA in deg
        dec (float): DEC in deg
        FRB (frb.FRB):
        z_frb (float, optional):
            Redshift of the host, expected to be provided

    """
    @classmethod
    def by_frb(cls, frb, **kwargs):
        """
        
        Args:
            frb (FRB):  FRB object
            **kwargs: 

        Returns:
            FRBHost:

        """
        # Strip off FRB
        if frb.frb_name[0:3] == 'FRB':
            name = frb.frb_name[3:]
        else:
            name = frb.frb_name
        #
        path = os.path.join(resource_filename('frb', 'data/Galaxies/'), name)
        json_file = os.path.join(path, FRBHost._make_outfile(name))
        slf = cls.from_json(frb, json_file, **kwargs)
        return slf

    def __init__(self, ra, dec, frb, z_frb=None, **kwargs):
        # Instantiate
        super(FRBHost, self).__init__(ra, dec, frb, **kwargs)

        # Name
        if frb.frb_name[0:3] == 'FRB':
            name = frb.frb_name[3:]
        else:
            name = frb.frb_name
        self.name = 'HG{}'.format(name)

        # Optional
        if z_frb is not None:
            self.redshift['z_FRB'] = z_frb

    @staticmethod
    def _make_outfile(frbname):
        """
        Static method to generate outfile based on frbname

        Args:
            frbname (str):  FRB name, e.g. 121102 or FRB121102

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
        outfile = self._make_outfile(self.frb.frb_name)
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
    """
    Foreground galaxy class (child of FRBGalaxy)
    """

    def __init__(self, ra, dec, frb, **kwargs):
        # Instantiate
        super(FGGalaxy, self).__init__(ra, dec, frb, **kwargs)

        # Load up FRB info from name
        self.name = 'FG{}_{}'.format(self.frb, utils.name_from_coord(self.coord))



