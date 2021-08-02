""" Top-level module to build or re-build the JSON files for
FRB host galaxies"""

from pkg_resources import resource_filename
import os
import sys
import warnings

from IPython import embed

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table
from astropy.coordinates import match_coordinates_sky

from frb.frb import FRB
from frb.galaxies import frbgalaxy, defs, offsets
from frb.galaxies import photom as frbphotom
try:
    from frb.galaxies import ppxf
except:
    print('WARNING:  ppxf not installed')
from frb.galaxies import nebular
from frb.surveys import des
from frb.surveys import sdss
from frb.surveys import wise
from frb.surveys import panstarrs
from frb.surveys import catalog_utils
import pandas

try:
    import extinction
except ImportError:
    print("extinction package not loaded.  Extinction corrections will fail")

try:
    from frb.galaxies import cigale
except ModuleNotFoundError:
    warnings.warn("You haven't installed pcigale and won't be able to do that analysis")

try:
    from frb.galaxies import eazy as frbeazy
except:
    warnings.warn("NOT READY FOR EAZY!")

from linetools.spectra.xspectrum1d import XSpectrum1D

db_path = os.getenv('FRB_GDB')
if db_path is None:
    print("Warning, you need to set $FRB_GDB to build hosts")
    #embed(header='You need to set $FRB_GDB')

ebv_method = 'SandF'


# New astrometry
mannings2021_astrom = pandas.read_csv(os.path.join(resource_filename('frb','data'),
                                          'Galaxies','Additional','Mannings2021', 
                                          'astrometry_v2.csv'))
# Probably will rename this                                        
mannings2021_astrom = mannings2021_astrom[
    (mannings2021_astrom.Filter == 'F160W') | (
        mannings2021_astrom.Filter == 'F110W')].copy()



def build_host_121102(build_photom=False, build_cigale=False, use_orig=False):
    """
    Generate the JSON file for FRB 121102

    Writes to 121102/FRB121102_host.json

    The majority of data comes from Tendulkar et al. 2017

    reff from Bassa+17 https://ui.adsabs.harvard.edu/abs/2017ApJ...843L...8B/abstract


    Args:
        build_photom (bool, optional): Generate the photometry file in the Galaxy_DB

    """
    #FRB_coord = SkyCoord('05h31m58.6980s +33d8m52.671s', frame='icrs') # Updated coords from Bassa+17
    # Eyeball Tendulkar+17 PA
    # UPDATE THIS!
    #gal_coord = FRB_coord.directional_offset_by(-45 * units.deg, 286e-3 * units.arcsec)
    gal_coord = SkyCoord('05h31m58.6980s +33d8m52.671s', frame='icrs') # Updated coords from Bassa+17
    
    # Instantiate
    frb121102 = FRB.by_name('FRB121102')
    host121102 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb121102)

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, host121102)

    # Redshift
    host121102.set_z(0.19273, 'spec', err=0.00008)

    # Photometry
    EBV = nebular.get_ebv(gal_coord, definition=ebv_method)['meanValue']  #
    print("EBV={} for the host {}".format(EBV, host121102.name))

    # photom_file = os.path.join(db_path, 'Repeater', 'Tendulkar2017', 'tendulkar2017_photom.ascii')
    photom_file = os.path.join(db_path, 'Repeater', 'Bassa2017', 'bassa2017_photom.ascii')
    if build_photom:
        photom = Table()
        photom['Name'] = ['HG121102']
        photom['ra'] = [host121102.coord.ra.value]
        photom['dec'] = host121102.coord.dec.value
        #
        # GMOS-N photometry from Bassa+2017, different than Tendulkar+2017, see Bassa+17 for details (AB)
        photom['GMOS_N_r'] = 25.46
        photom['GMOS_N_r_err'] = 0.14
        photom['GMOS_N_i'] = 24.75
        photom['GMOS_N_i_err'] = 0.09
        photom['GMOS_N_z'] = 24.30
        photom['GMOS_N_z_err'] = 0.13
        photom['GMOS_N_g'] = 25.85
        photom['GMOS_N_g_err'] = 0.12

        # HST from Bassa+2017 (AB)
        photom['WFC3_F110W'] = 23.675
        photom['WFC3_F110W_err'] = 0.012
        photom['WFC3_F160W'] = 23.31
        photom['WFC3_F160W_err'] = 0.03

        # Spitzer from Bassa (micro-Jy)
        mag_3p6, err_mag_3p6 = catalog_utils.mag_from_flux(flux=1.03e-3*units.mJy,
                                                           flux_err=0.19e-3*units.mJy)  # in micro-Jy
        photom['Spitzer_3.6'] = mag_3p6
        photom['Spitzer_3.6_err'] = err_mag_3p6
        photom['Spitzer_4.5'] = catalog_utils.mag_from_flux(0.9e-3*units.mJy/2.)[0]   # upper limit (6sigma/2 = ~3sigma) in micro-Jy
        photom['Spitzer_4.5_err'] = 999.  # the flux is un upper limit, note it is 3sigma (estimated by dividing the 6sigma/2)

        # Write
        photom = frbphotom.merge_photom_tables(photom, photom_file)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
    # Read
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correction
    frbphotom.correct_photom_table(photom, EBV, 'HG121102')
    # Parse
    host121102.parse_photom(photom, EBV=EBV)


    # CIGALE
    cigale_file = os.path.join(db_path, 'Repeater', 'Bassa2017', 'HG121102_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        # Run
        cigale.host_run(host121102, cigale_file=cigale_file)

    host121102.parse_cigale(cigale_file, sfh_file=sfh_file)

    # Nebular lines
    neb_lines = {}
    neb_lines['Halpha'] = 0.652e-16
    neb_lines['Halpha_err'] = 0.009e-16
    neb_lines['Halpha_Al'] = 0.622
    #
    neb_lines['Hbeta'] = 0.118e-16
    neb_lines['Hbeta_err'] = 0.011e-16
    neb_lines['Hbeta_Al'] = 0.941
    #
    neb_lines['[OIII] 5007'] = 0.575e-16
    neb_lines['[OIII] 5007_err'] = 0.011e-16
    neb_lines['[OIII] 5007_Al'] = 0.911
    #
    neb_lines['[NII] 6584'] = 0.030e-16  # * units.erg/units.cm**2/units.s      # Upper limit
    neb_lines['[NII] 6584_err'] = -999.  # * units.erg/units.cm**2/units.s      # Upper limit
    neb_lines['[NII] 6584_Al'] = 0.619

    AV = 2.42

    # Extinction correct
    for key in neb_lines.keys():
        if '_err' in key:
            continue
        if 'Al' in key:
            continue
        # Ingest
        host121102.neb_lines[key] = neb_lines[key] * 10 ** (neb_lines[key + '_Al'] * AV / 2.5)
        if neb_lines[key + '_err'] > 0:
            host121102.neb_lines[key + '_err'] = neb_lines[key + '_err'] * 10 ** (neb_lines[key + '_Al'] * AV / 2.5)
        else:
            host121102.neb_lines[key + '_err'] = neb_lines[key + '_err']

    # AV
    host121102.calc_nebular_AV('Ha/Hb', min_AV=0.)

    # SFR
    host121102.calc_nebular_SFR('Ha')

    # Vette
    for key in host121102.neb_lines.keys():
        if '_err' in key:
            continue
        assert key in defs.valid_neb_lines

    '''
    # Morphology : Bassa+2017 half-light
    host121102.morphology['reff_ang'] = 0.20   # arcsec
    host121102.morphology['reff_ang_err'] = 0.01
    host121102.morphology['reff_kpc'] = 0.66   # kpc
    host121102.morphology['reff_kpc_err'] = 0.03
    # Other
    host121102.morphology['n'] = 2.2
    host121102.morphology['n_err'] = 1.5
    #
    host121102.morphology['b/a'] = 0.25
    host121102.morphology['b/a_err'] = 0.13
    '''

    # Galfit -- Mannings+2021
    host121102.parse_galfit(os.path.join(db_path, 'F4', 'mannings2020',
                                   'HG121102_galfit.fits'))

    # Derived quantities
    if use_orig:
        host121102.derived['M_r'] = -17.0  # AB; Tendulkar+17
        host121102.derived['M_r_err'] = 0.2  # Estimated by JXP
        host121102.derived['SFR_nebular'] = 0.23  # from Tendulkar+17
        host121102.derived['Mstar'] = 5.5e7  # Msun; Tendulkar+17
        host121102.derived['Mstar_err'] = 1.5e7  # Msun; Tendulkar+17
    host121102.derived['Z_spec'] = -0.16  # Tendulkar+17 on a scale with Solar O/H = 8.86
    host121102.derived['Z_spec_err'] = -999.  # Tendulkar+17

    # References
    host121102.refs = ['Tendulkar2017, Bassa2017']

    # Vet
    assert host121102.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/121102')
    host121102.write_to_json(path=path, overwrite=True)


def build_host_180924(build_photom=False, build_cigale=False):
    """
    Generate the JSON file for FRB 180924
    
    All data are from Bannister et al. 2019
        https://ui.adsabs.harvard.edu/abs/2019Sci...365..565B/abstract

    Writes to 180924/FRB180924_host.json

    Args:
        build_photom (bool, optional): Generate the photometry file in the Galaxy_DB
    """
    frbname = '180924'
    gal_coord = SkyCoord('J214425.25-405400.81', unit=(units.hourangle, units.deg))

    # Instantiate
    frb180924 = FRB.by_name('FRB180924')
    host180924 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb180924)

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, host180924)

    # Redshift -- JXP measured from multiple data sources
    host180924.set_z(0.3212, 'spec')

    # Photometry

    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'CRAFT', 'Bannister2019', 'bannister2019_photom.ascii')
    if build_photom:
        # DES
        search_r = 2 * units.arcsec
        des_srvy = des.DES_Survey(gal_coord, search_r)
        des_tbl = des_srvy.get_catalog(print_query=True)
        host.parse_photom(des_tbl)  # This is a bit of a hack!

        photom = Table()
        photom['Name'] = ['HG{}'.format(frbname)]
        photom['ra'] = host180924.coord.ra.value
        photom['dec'] = host180924.coord.dec.value
        photom['VLT_FORS2_g'] = 21.38  # No extinction correction
        photom['VLT_FORS2_g_err'] = 0.04
        photom['VLT_FORS2_I'] = 20.10
        photom['VLT_FORS2_I_err'] = 0.02
        # Add in DES
        for key in host.photom.keys():
            photom[key] = host.photom[key]
        # Merge/write
        photom = frbphotom.merge_photom_tables(photom, photom_file)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)

    # Load
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correction
    EBV = nebular.get_ebv(gal_coord, definition=ebv_method)['meanValue']
    frbphotom.correct_photom_table(photom, EBV, 'HG{}'.format(frbname))
    # Parse
    host180924.parse_photom(photom, EBV=EBV)

    # PPXF
    host180924.parse_ppxf(os.path.join(db_path, 'CRAFT', 'Bannister2019', 'HG180924_MUSE_ppxf.ecsv'))

    # Derived quantities

    # AV
    host180924.calc_nebular_AV('Ha/Hb')

    # SFR
    host180924.calc_nebular_SFR('Ha')
    host180924.derived['SFR_nebular_err'] = -999.

    # CIGALE
    cigale_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG180924_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        cut_photom = Table()
        for key in host.photom.keys():
            if 'DES' not in key and 'WISE' not in key:
                continue
            cut_photom[key] = [host.photom[key]]
        cigale.host_run(host, cut_photom=cut_photom, cigale_file=cigale_file)
    # Parse
    host180924.parse_cigale(cigale_file, sfh_file=sfh_file)

    # Galfit
    #host180924.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
    #                               'HG180924_DES_i_galfit.fits'))
    # Galfit -- Mannings+2021
    host180924.parse_galfit(os.path.join(db_path, 'F4', 'mannings2020',
                                   'HG180924_galfit.fits'))

    # Vet all
    host180924.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/180924')
    host180924.write_to_json(path=path)


def build_host_181112(build_photom=False, build_cigale=False):
    """ Build the host galaxy data for FRB 181112

    All of the data comes from Prochaska+2019, Science

    Args:
        build_photom (bool, optional):

    """
    frbname = '181112'
    FRB_coord = SkyCoord('J214923.63-525815.39',
                         unit=(units.hourangle, units.deg))  # Cherie;  2019-04-17 (Slack)
    # Coord from DES
    gal_coord = SkyCoord('J214923.66-525815.28',
                          unit=(units.hourangle, units.deg))  # from DES

    # Instantiate
    frb181112 = FRB.by_name('FRB181112')
    host181112 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb181112)
    host181112.frb_coord = FRB_coord

    # Redshift
    host181112.set_z(0.4755, 'spec', err=7e-5)

    # ############
    # Photometry


    # VLT -- Lochlan 2019-05-02
    # VLT -- Lochlan 2019-06-18
    photom_file = os.path.join(db_path, 'CRAFT', 'Prochaska2019', 'prochaska2019_photom.ascii')
    if build_photom:
        # DES
        # Grab the table (requires internet)
        search_r = 2 * units.arcsec
        des_srvy = des.DES_Survey(Host_coord, search_r)
        des_tbl = des_srvy.get_catalog(print_query=True)
        host181112.parse_photom(des_tbl) # A hack of sorts

        photom = Table()
        photom['Name'] = ['HG{}'.format(frbname)]
        photom['ra'] = host181112.coord.ra.value
        photom['dec'] = host181112.coord.dec.value
        photom['VLT_FORS2_g'] = 22.57  # No extinction correction
        photom['VLT_FORS2_g_err'] = 0.04
        photom['VLT_FORS2_I'] = 21.51
        photom['VLT_FORS2_I_err'] = 0.04
        # Add in DES
        for key in host181112.photom.keys():
            photom[key] = host181112.photom[key]
        # Merge/write
        photom = frbphotom.merge_photom_tables(photom, photom_file)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
    # Load
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correction
    EBV = nebular.get_ebv(host181112.coord, definition=ebv_method)['meanValue']
    frbphotom.correct_photom_table(photom, EBV, 'HG{}'.format(frbname))
    # Parse
    host181112.parse_photom(photom, EBV=EBV)

    # Nebular lines
    host181112.parse_ppxf(os.path.join(db_path, 'CRAFT', 'Prochaska2019', 'HG181112_FORS2_ppxf.ecsv'))

    # Adjust errors on Ha, [NII] because of telluric

    # Derived quantities
    host181112.calc_nebular_AV('Ha/Hb', min_AV=0.)

    # Ha is tough in telluric
    host181112.calc_nebular_SFR('Hb', AV=0.15)  # Photometric
    # This would be an upper limit
    #host.calc_nebular_SFR('Ha')

    # CIGALE
    cigale_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG181112_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        cut_photom = Table()
        for key in host181112.photom.keys():
            if 'DES' not in key: #and 'WISE' not in key:
                continue
            cut_photom[key] = [host181112.photom[key]]
        cigale.host_run(host181112, cut_photom=cut_photom, cigale_file=cigale_file)
    # Parse
    host181112.parse_cigale(cigale_file, sfh_file=sfh_file)

    # Galfit
    host181112.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
                                   'HG181112_VLT_i_galfit.fits'))

    # Write
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host181112.write_to_json(path=path)


def build_host_190102(build_photom=False, build_cigale=False,
                      build_ppxf=False):
    """ Build the host galaxy data for FRB 190102

    All of the data comes from Bhandrari+2020, ApJL, in press

    Args:
        build_photom (bool, optional):
    """
    frbname = '190102'
    # Stuart on June 17, 2019
    #  -- Astrometry.net !
    gal_coord = SkyCoord('J212939.60-792832.4',
                         unit=(units.hourangle, units.deg))  # Cherie;  07-Mar-2019
    # Instantiate
    frb190102 = FRB.by_name('FRB190102')
    host190102 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb190102)

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, host190102)

    # Redshift -- Gaussian fit to [OIII 5007] in MagE
    #  Looks great on the other lines
    #  Ok on FORS2 Halpha
    wv_oiii = 6466.48
    z_OIII = wv_oiii / 5008.239 - 1
    host190102.set_z(z_OIII, 'spec')

    photom_file1 = os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'bhandari2019_photom.ascii')
    photom_file2 = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'heintz2020_photom.ascii')
    # VLT/FORS2 -- Pulled from draft on 2019-06-23
    # VLT/FORS2 -- Pulled from spreadsheet 2019-06-23
    if build_photom:
        '''   Bhandari2020
        photom = Table()
        photom['ra'] = [host190102.coord.ra.value]
        photom['dec'] = host190102.coord.dec.value
        photom['Name'] = host190102.name
        photom['VLT_FORS2_u'] = 23.7   # Not dust corrected
        photom['VLT_FORS2_u_err'] = 0.2  # -999.
        photom['VLT_FORS2_g'] = 22.6
        photom['VLT_FORS2_g_err'] = 0.1
        photom['VLT_FORS2_I'] = 21.1
        photom['VLT_FORS2_I_err'] = 0.05
        photom['VLT_FORS2_z'] = 20.8
        photom['VLT_FORS2_z_err'] = 0.2
        '''
        # Photometry from Mannings et al. 2020-09-02
        photom = Table()
        photom['ra'] = [host190102.coord.ra.value]
        photom['dec'] = host190102.coord.dec.value
        photom['Name'] = host190102.name
        photom['WFC3_F160W'] = 20.550
        photom['WFC3_F160W_err'] = 0.006
        # LCOGT?
        # Write
        photom = frbphotom.merge_photom_tables(photom, photom_file2)
        photom.write(photom_file2, format=frbphotom.table_format, overwrite=True)
    # Load
    photom = frbphotom.photom_by_name(host190102.name, [photom_file1, photom_file2])
    # Dust correction
    EBV = nebular.get_ebv(gal_coord, definition=ebv_method)['meanValue']
    frbphotom.correct_photom_table(photom, EBV, host190102.name)
    # Parse
    host190102.parse_photom(photom, EBV=EBV)

    # PPXF
    if build_ppxf:
        # MagE
        results_file = os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'HG190102_MagE_ppxf.ecsv')
        spec_fit = os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'HG190102_MagE_ppxf.fits')
        meta, spectrum = host190102.get_metaspec(instr='MagE')
        R = meta['R']

        # Correct for Galactic extinction
        ebv = float(nebular.get_ebv(host190102.coord)['meanValue'])
        AV = ebv * 3.1  # RV
        Al = extinction.ccm89(spectrum.wavelength.value, AV, 3.1)
        # New spec
        new_flux = spectrum.flux * 10**(Al/2.5)
        new_sig = spectrum.sig * 10**(Al/2.5)
        new_spec = XSpectrum1D.from_tuple((spectrum.wavelength, new_flux, new_sig))

        # Mask
        atmos = [(7550, 7750)]
        ppxf.run(new_spec, R, host190102.z, results_file=results_file, spec_fit=spec_fit, chk=True, atmos=atmos)
    host190102.parse_ppxf(os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'HG190102_MagE_ppxf.ecsv'))

    # Derived quantities

    # AV
    host190102.calc_nebular_AV('Ha/Hb')

    # SFR
    host190102.calc_nebular_SFR('Ha')

    # CIGALE
    cigale_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG190102_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        cigale.host_run(host190102, cigale_file=cigale_file)

    host190102.parse_cigale(cigale_file, sfh_file=sfh_file)

    # Galfit
    #host190102.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
    #                               'HG190102_VLT_i_galfit.fits'))
    # Galfit -- Mannings+2021
    host190102.parse_galfit(os.path.join(db_path, 'F4', 'mannings2020',
                                   'HG190102_galfit.fits'))

    # Vet all
    host190102.vet_all()

    # Write 
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host190102.write_to_json(path=path)


def build_host_190523(build_photom=False, build_cigale=False):  #:run_ppxf=False, build_photom=False):
    """
    Build the host galaxy data for FRB 190523

    Most of the data is from Ravi+2019
        https://ui.adsabs.harvard.edu/abs/2019Natur.572..352R/abstract

    The exception is that CRAFT (S. Simha) have run CIGALE on the photometry for
    a consistent analysis with the ASKAP hosts.


    Args:
        build_photom:

    Returns:

    """
    frbname = '190523'
    S1_gal_coord = SkyCoord(ra=207.06433, dec=72.470756, unit='deg')  # Pan-STARRs; J134815.4392+722814.7216

    # Instantiate
    frb190523 = FRB.by_name('FRB190523')
    host190523_S1 = frbgalaxy.FRBHost(S1_gal_coord.ra.value, S1_gal_coord.dec.value, frb190523)
    host190523_S1.name += '_S1'

    # Load redshift table
    host190523_S1.set_z(0.660, 'spec')

    # Morphology

    # Photometry
    EBV = nebular.get_ebv(S1_gal_coord, definition=ebv_method)['meanValue']  #
    print("EBV={} for the host of {}".format(EBV, frbname))

    # PanStarrs
    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'DSA', 'Ravi2019', 'ravi2019_photom.ascii')
    if build_photom:
        search_r = 1 * units.arcsec
        ps_srvy = panstarrs.Pan_STARRS_Survey(S1_gal_coord, search_r)
        ps_tbl = ps_srvy.get_catalog(print_query=True)
        ps_tbl['Name'] = host190523_S1.name
        photom = frbphotom.merge_photom_tables(ps_tbl, photom_file)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
    # Read
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correction
    frbphotom.correct_photom_table(photom, EBV, 'HG190523_S1')
    # Parse
    host190523_S1.parse_photom(photom, EBV=EBV)

    # PPXF
    '''
    if run_ppxf:
        results_file = os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'HG190608_SDSS_ppxf.ecsv')
        meta, spectrum = host190608.get_metaspec(instr='SDSS')
        spec_fit = None
        ppxf.run(spectrum, 2000., host190608.z, results_file=results_file, spec_fit=spec_fit, chk=True)
    host190608.parse_ppxf(os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'HG190608_SDSS_ppxf.ecsv'))
    '''

    # CIGALE -- PanStarrs photometry but our own CIGALE analysis
    cigale_file = os.path.join(db_path, 'DSA', 'Ravi2019', 'S1_190523_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        cigale.host_run(host190523_S1, cigale_file=cigale_file)
    # Parse
    host190523_S1.parse_cigale(cigale_file, sfh_file=sfh_file)

    # Nebular flux measured by a hand (Gaussian fit) by JXP on 2020-05-19
    #   Corrected for Galactic extinction but not internal
    neb_lines = {}
    neb_lines['Hbeta'] = 3e-18
    neb_lines['Hbeta_err'] = -999

    host190523_S1.neb_lines = neb_lines

    # SFR
    host190523_S1.calc_nebular_SFR('Hb')

    # Derived quantities
    #host190523_S1.derived['SFR_nebular'] = 1.3
    #host190523_S1.derived['SFR_nebular_err'] = -999.

    # Galfit
    host190523_S1.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
                                   'HG190523_LRIS_r_galfit.fits'))
    # Vet all
    host190523_S1.vet_all()

    # Write 
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host190523_S1.write_to_json(path=path)


def build_host_190608(run_ppxf=False, build_photom=False, build_cigale=False):
    """ Build the host galaxy data for FRB 190608

    All of the data comes from Bhandrari+2020, ApJL, in press

    Args:
        build_photom (bool, optional):
    """
    frbname = '190608'
    gal_coord = SkyCoord('J221604.90-075356.0',
                         unit=(units.hourangle, units.deg))  # Cherie;  07-Mar-2019

    # Instantiate
    frb190608 = FRB.by_name('FRB190608')
    host190608 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb190608)

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, host190608)

    # Load redshift table
    ztbl = Table.read(os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'z_SDSS.ascii'),
                      format='ascii.fixed_width')
    z_coord = SkyCoord(ra=ztbl['RA'], dec=ztbl['DEC'], unit='deg')
    idx, d2d, _ = match_coordinates_sky(gal_coord, z_coord, nthneighbor=1)
    if np.min(d2d) > 0.5*units.arcsec:
        embed(header='190608')
    # Redshift -- SDSS
    host190608.set_z(ztbl[idx]['ZEM'], 'spec')

    # Morphology
    #host190608.parse_galfit(os.path.join(photom_path, 'CRAFT', 'Bannister2019',
    #                               'HG180924_galfit_DES.log'), 0.263)

    # Photometry

    # SDSS
    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'bhandari2019_photom.ascii')
    if build_photom:
        # VLT
        search_r = 1 * units.arcsec
        # SDSS
        sdss_srvy = sdss.SDSS_Survey(gal_coord, search_r)
        sdss_tbl = sdss_srvy.get_catalog(print_query=True)
        sdss_tbl['Name'] = 'HG190608'
        photom = frbphotom.merge_photom_tables(sdss_tbl, photom_file)
        # WISE
        wise_srvy = wise.WISE_Survey(gal_coord, search_r)
        wise_tbl = wise_srvy.get_catalog(print_query=True)
        wise_tbl['Name'] = 'HG190608'
        # Write
        photom = frbphotom.merge_photom_tables(wise_tbl, photom, debug=True)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
    # Load
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correction
    EBV = nebular.get_ebv(gal_coord, definition=ebv_method)['meanValue']
    frbphotom.correct_photom_table(photom, EBV, 'HG{}'.format(frbname))
    # Parse
    host190608.parse_photom(photom, EBV=EBV)

    # PPXF
    results_file = os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'HG190608_SDSS_ppxf.ecsv')
    if run_ppxf:
        meta, spectrum = host190608.get_metaspec(instr='SDSS')
        spec_fit = None
        ppxf.run(spectrum, 2000., host190608.z, results_file=results_file, spec_fit=spec_fit, chk=True)
    host190608.parse_ppxf(results_file)

    # Derived quantities

    # AV
    host190608.calc_nebular_AV('Ha/Hb')

    # SFR
    host190608.calc_nebular_SFR('Ha')
    #host.derived['SFR_nebular_err'] = -999.

    # CIGALE
    cigale_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG190608_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        cigale.host_run(host190608, cigale_file=cigale_file)
    # Parse
    host190608.parse_cigale(cigale_file, sfh_file=sfh_file)

    # Galfit
    #host190608.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
    #                               'HG190608_SDSS_i_galfit.fits'))
    # Galfit -- Mannings+2021
    host190608.parse_galfit(os.path.join(db_path, 'F4', 'mannings2020',
                                   'HG190608_galfit.fits'))
    # Vet all
    host190608.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host190608.write_to_json(path=path)


def build_host_180916(run_ppxf=False, build_photom=False, build_cigale=False):
    """ Build the host galaxy data for FRB 180916

    All of the data currently comes from Marcote et al. 2020
    https://ui.adsabs.harvard.edu/abs/2020Natur.577..190M/abstract

    CIGALE and galfit will be published in Kasper et al. 2020

    Args:
        build_photom (bool, optional):
    """
    frbname = '180916'
    gal_coord = SkyCoord('J015800.28+654253.0',
                         unit=(units.hourangle, units.deg))

    # Instantiate
    frb180916 = FRB.by_name('FRB180916')
    host180916 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb180916)

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, host180916)

    # Redshift
    host180916.set_z(0.0337, 'spec')

    # Morphology
    #host190608.parse_galfit(os.path.join(photom_path, 'CRAFT', 'Bannister2019',
    #                               'HG180924_galfit_DES.log'), 0.263)

    # Photometry
    EBV = nebular.get_ebv(gal_coord, definition=ebv_method)['meanValue']  #
    print("EBV={} for the host of {}".format(EBV, frbname))

    # SDSS
    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'CHIME', 'Marcote2020', 'marcote2020_photom.ascii')
    if build_photom:
        # SDSS
        search_r = 1 * units.arcsec
        sdss_srvy = sdss.SDSS_Survey(gal_coord, search_r)
        sdss_tbl = sdss_srvy.get_catalog(print_query=True)
        sdss_tbl['Name'] = host180916.name
        photom = frbphotom.merge_photom_tables(sdss_tbl, photom_file)
        # WISE
        wise_srvy = wise.WISE_Survey(gal_coord, search_r)
        wise_tbl = wise_srvy.get_catalog(print_query=True)
        wise_tbl['Name'] = host180916.name
        photom = frbphotom.merge_photom_tables(wise_tbl, photom, debug=True)
        # Write
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
    # Read
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correction
    frbphotom.correct_photom_table(photom, EBV, 'HG180916')
    # Parse
    host180916.parse_photom(photom, EBV=EBV)

    # CIGALE
    cigale_file = os.path.join(db_path, 'CHIME', 'Marcote2020', 'HG180916_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        # Run
        cigale.host_run(host180916, cigale_file=cigale_file)

    host180916.parse_cigale(cigale_file, sfh_file=sfh_file)

    # PPXF
    #results_file = os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'HG190608_SDSS_ppxf.ecsv')
    #if run_ppxf:
    #    meta, spectrum = host190608.get_metaspec(instr='SDSS')
    #    spec_fit = None
    #    ppxf.run(spectrum, 2000., host190608.z, results_file=results_file, spec_fit=spec_fit, chk=True)
    #host190608.parse_ppxf(results_file)

    # Nebular lines
    host180916.neb_lines['Halpha'] = 6.57e-16  # Not corrected for Galactic extinction
    host180916.neb_lines['Halpha_err'] = 0.04e-16
    host180916.neb_lines['[OIII] 5007'] = 4.76e-16  # Not corrected for Galactic extinction
    host180916.neb_lines['[OIII] 5007_err'] = 0.04e-16
    host180916.neb_lines['[NII] 6584'] = 2.51e-16  # Not corrected for Galactic extinction
    host180916.neb_lines['[NII] 6584_err'] = 0.04e-16
    host180916.neb_lines['[OIII] 4959'] = 0.38e-16  # Not corrected for Galactic extinction
    host180916.neb_lines['[OIII] 4959_err'] = 0.04e-16

    # Correct for extinction
    RV = 3.1
    AV = EBV * RV
    Alambda = extinction.fm07(np.array([6564*(1+host180916.z)]), AV)[0]
    host180916.neb_lines['Halpha'] *= 10**(Alambda/2.5)
    host180916.neb_lines['Halpha_err'] *= 10**(Alambda/2.5)
    Alambda = extinction.fm07(np.array([5007*(1+host180916.z)]), AV)[0]
    host180916.neb_lines['[OIII] 5007'] *= 10**(Alambda/2.5)
    host180916.neb_lines['[OIII] 5007_err'] *= 10**(Alambda/2.5)
    Alambda = extinction.fm07(np.array([6584*(1+host180916.z)]), AV)[0]
    host180916.neb_lines['[NII] 6584'] *= 10**(Alambda/2.5)
    host180916.neb_lines['[NII] 6584_err'] *= 10**(Alambda/2.5)
    Alambda = extinction.fm07(np.array([4959*(1+host180916.z)]), AV)[0]
    host180916.neb_lines['[OIII] 4959'] *= 10**(Alambda/2.5)
    host180916.neb_lines['[OIII] 4959_err'] *= 10**(Alambda/2.5)

    # Derived quantities
    # AV
    #host190608.calc_nebular_AV('Ha/Hb')


    # SFR
    host180916.calc_nebular_SFR('Ha')

    # Galfit
    #host180916.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
    #                               'HG180916_SDSS_i_galfit.fits'))
    # Galfit -- Mannings+2021
    host180916.parse_galfit(os.path.join(db_path, 'F4', 'mannings2020',
                                   'HG180916_galfit.fits'))

    # Vet all
    assert host180916.vet_all()

    # Write 
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host180916.write_to_json(path=path)


def build_host_190611(run_ppxf=False, build_photom=False, build_cigale=False, source='faint'):
    """ Build the host galaxy data for FRB 190611

    There are 2 sources in play.

    Heintz+2020

    Args:
        build_photom (bool, optional):
    """
    frbname = '190611'
    bright_coord = SkyCoord("21h22m58.277s -79d23m50.09s", frame='icrs')  # Kasper on 2020-04-20
    faint_coord = SkyCoord("21h22m58.973s  -79d23m51.69s", frame='icrs')  # Kasper on 2020-04-20

    if source == 'faint':
        gal_coord = faint_coord
    else:
        gal_coord = bright_coord

    # Instantiate
    frb190611 = FRB.by_name('FRB190611')
    host190611 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb190611)

    # Redshift
    if source == 'faint':
        host190611.set_z(0.999, 'phot')
    else:
        host190611.set_z(0.3778, 'spec')
        host190611.name = 'HG190611b'

    # Morphology
    #host190608.parse_galfit(os.path.join(photom_path, 'CRAFT', 'Bannister2019',
    #                               'HG180924_galfit_DES.log'), 0.263)

    # Photometry
    EBV = nebular.get_ebv(gal_coord)['meanValue']  #
    print("EBV={} for the host of {}".format(EBV, frbname))

    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'heintz2020_photom.ascii')
    if build_photom:
        photom = Table()
        photom['Name'] = [host190611.name]
        photom['ra'] = host190611.coord.ra.value
        photom['dec'] = host190611.coord.dec.value

        # These are observed, not corrected
        if source == 'faint':
            photom['GMOS_S_r'] = 26.30
            photom['GMOS_S_r_err'] = 0.25
            photom['GMOS_S_i'] = 25.70
            photom['GMOS_S_i_err'] = 0.25
        else:
            photom['GMOS_S_r'] = 22.65
            photom['GMOS_S_r_err'] = 0.15
            photom['GMOS_S_i'] = 22.75
            photom['GMOS_S_i_err'] = 0.15
        # Merge/write
        photom = frbphotom.merge_photom_tables(photom, photom_file)
        # Write
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
    # Parse
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correction
    frbphotom.correct_photom_table(photom, EBV, 'HG190611b')
    # Parse
    host190611.parse_photom(photom, EBV=EBV)

    # CIGALE
    cigale_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG190611_CIGALE.fits')
    # NOT ENOUGH PHOTOMETRY


    # PPXF
    ppxf_results_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', '{}_FORS2_ppxf.ecsv'.format(host190611.name))
    spec_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', '{}_FORS2_ppxf.fits'.format(host190611.name))
    if run_ppxf:
        meta, spectrum = host190611.get_metaspec(instr='FORS2')
        R = meta['R']
        print("R = {}".format(R))
        # Correct for Galactic extinction
        AV = EBV * 3.1  # RV
        Al = extinction.fm07(spectrum.wavelength.value, AV)#, 3.1)
        # New spec
        new_flux = spectrum.flux * 10**(Al/2.5)
        new_sig = spectrum.sig * 10**(Al/2.5)
        new_spec = XSpectrum1D.from_tuple((spectrum.wavelength, new_flux, new_sig))
        # Mask
        ppxf.run(new_spec, R, host190611.z, results_file=ppxf_results_file, spec_fit=spec_file,
             atmos=[[0., 6400.], [9300., 20000]], chk=True)
    host190611.parse_ppxf(ppxf_results_file)

    # Derived quantities

    # AV
    host190611.calc_nebular_AV('Ha/Hb')

    # SFR
    host190611.calc_nebular_SFR('Ha')
    #host.derived['SFR_nebular_err'] = -999.

    # Galfit
    host190611.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
                                         'HG190611_GMOS_i_galfit.fits'))

    # Vet all
    assert host190611.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host190611.write_to_json(path=path)


def build_host_190614(build_photom=False, build_cigale=False, run_eazy=False,
            build_A=False, build_B=False, build_C=False):
    """ Build the host galaxy data for FRB 190614
    See Law+2020 https://ui.adsabs.harvard.edu/abs/2020ApJ...899..161L/abstract

    Args:
        build_photom (bool, optional):
        build_cigale (bool, optional):
        run_ppxf (bool, optional):
    """
    frbname = '190614'
    eazy_folder = './eazy'

    frb190614 = FRB.by_name('FRB190614')

    #########################################################
    # A
    #########################################################
    if build_A:
        print("Building Host galaxy A for FRB{}".format(frbname))
        gal_coord = SkyCoord(ra=65.073804, dec=73.706356,  # Slack on 2020 Apr 03
                             unit='deg')  # J042017.713+734222.88

        # Instantiate
        host190614A = frbgalaxy.FRBHost(gal_coord.ra.value, 
                                        gal_coord.dec.value, frb190614)
        host190614A.name = 'G190614_A'

        # Photometry
        EBV = nebular.get_ebv(gal_coord)['meanValue']  # 0.1401

        # Grab the table (requires internet)
        photom_file = os.path.join(db_path, 'Realfast', 'Law2020', 'law2020_photom.ascii')
        gname = 'G{}_A'.format(frbname)
        if build_photom:
            photom = Table()
            photom['Name'] = [gname]
            photom['ra'] = host190614A.coord.ra.value
            photom['dec'] = host190614A.coord.dec.value
            # These are observed
            photom['LRISb_V'] = 25.86
            photom['LRISb_V_err'] = 0.25
            photom['GMOS_S_r'] = 23.61
            photom['GMOS_S_r_err'] = 0.15
            photom['LRISr_I'] = 23.09
            photom['LRISr_I_err'] = 0.1
            photom['NOT_z'] = 23.35
            photom['NOT_z_err'] = 0.3
            photom['NIRI_J'] = 21.75 + 0.91
            photom['NIRI_J_err'] = 0.2

            # Write
            # Merge/write
            photom = frbphotom.merge_photom_tables(photom, photom_file)
            photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
            print("Wrote photometry to: {}".format(photom_file))

        # Load
        photom = Table.read(photom_file, format=frbphotom.table_format)
        # Dust correction
        frbphotom.correct_photom_table(photom, EBV, gname)
        # Parse
        host190614A.parse_photom(photom, EBV=EBV)

        # EAZY
        pub_path = os.path.join(db_path, 'Realfast', 'Law2020')
        outputdir_base = 'EAZY_OUTPUT_{}'.format(host190614A.name)
        if run_eazy:
            frbeazy.eazy_input_files(host190614A.photom, os.path.join(eazy_folder, 'inputs'),
                                     host190614A.name,
                                     '../'+outputdir_base,
                                     templates='br07_default',
                                     prior_filter='GMOS_S_r')
            frbeazy.run_eazy(os.path.join(eazy_folder, 'inputs'),
                             host190614A.name,
                             os.path.join(eazy_folder, outputdir_base, 'logfile'))
            # Rename
            os.system('cp -rp {:s} {:s}'.format(os.path.join(eazy_folder, outputdir_base),
                                                pub_path))
        # Redshift
        zgrid, pzi, prior = frbeazy.getEazyPz(-1, MAIN_OUTPUT_FILE='photz',
                                             OUTPUT_DIRECTORY=os.path.join(pub_path, outputdir_base),
                                             CACHE_FILE='Same', binaries=None, get_prior=True)
        zphot, sig_zphot = frbeazy.eazy_stats(zgrid, pzi)

        host190614A.set_z(zphot, 'phot', err=sig_zphot)

        # CIGALE - Run by hand
        # CIGALE _ build now then run by hand after
        cigale_file = os.path.join(db_path, 'Realfast', 'Law2020', 'G190614_A_CIGALE.fits')
        if build_cigale:
            # Prep
            #embed(header='see 190614 to update this')
            cigale_tbl = photom.copy()

            cigale_tbl = Table(cigale_tbl[0])
            #pdb.set_trace()
            cigale_tbl['z'] = host190614A.z
            cigale_tbl['ID'] = 'G190614_A'
            # Run
            cigale.run(cigale_tbl, 'z', outdir='G190614_A', compare_obs_model=True)
            # Rename/move
            os.system('mv G190614_A/results.fits {:s}'.format(cigale_file))
            model_file = cigale_file.replace('CIGALE', 'CIGALE_model')
            os.system('mv G190614_A/{:s}_best_model.fits {:s}'.format('G190614_A', model_file))
            photo_file = cigale_file.replace('CIGALE.fits', 'CIGALE_photo.dat')
            os.system('mv G190614_A/photo_observed_model_G190614_A.dat {:s}'.format(photo_file))

        # Parse
        host190614A.parse_cigale(cigale_file)


        # Vet all
        assert host190614A.vet_all()

        # Write
        path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
        host190614A.write_to_json(path=path, outfile='G190614_A.json')


    #########################################################
    # B
    #########################################################
    if build_B:
        print("Building Host galaxy B for FRB{}".format(frbname))
        gal_coord = SkyCoord(ra=65.074467, dec=73.706783,  # Kasper on 2020 Apr 03
                             unit='deg')   # J042017.872+734224.42

        # Instantiate
        host190614B = frbgalaxy.FRBHost(gal_coord.ra.value, 
                                        gal_coord.dec.value, frb190614)
        host190614B.name = 'G190614_B'

        # Redshift -- JXP measured from FORS2
        #    Should be refined
        host190614B.set_z(0.66, 'phot', err=0.3)

        # Photometry
        EBV = nebular.get_ebv(gal_coord)['meanValue']  # 0.1401
        # Grab the table (requires internet)
        photom_file = os.path.join(db_path, 'Realfast', 'Law2020', 'law2020_photom.ascii')
        gname = 'G{}_B'.format(frbname)
        if build_photom:

            photom = Table()
            photom['Name'] = [gname]
            photom['ra'] = host190614B.coord.ra.value
            photom['dec'] = host190614B.coord.dec.value
            # These are observed
            photom['GMOS_S_r'] = 24.3
            photom['GMOS_S_r_err'] = 0.24
            photom['NOT_z'] = 22.7          # 3-sigma upper limit
            photom['NOT_z_err'] = 999.
            photom['LRISb_V'] = 25.02
            photom['LRISb_V_err'] = 0.16
            photom['LRISr_I'] = 24.00
            photom['LRISr_I_err'] = 0.18
            photom['NIRI_J'] = 22.85 + 0.91  # 3-sigma upper limit
            photom['NIRI_J_err'] = 999.

            # Write
            # Merge/write
            photom = frbphotom.merge_photom_tables(photom, photom_file)
            photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
            print("Wrote photometry to: {}".format(photom_file))

        # Parse
        photom = Table.read(photom_file, format=frbphotom.table_format)
        # Dust correction
        frbphotom.correct_photom_table(photom, EBV, gname)
        #
        host190614B.parse_photom(photom, EBV=EBV)

        # EAZY
        pub_path = os.path.join(db_path, 'Realfast', 'Law2020')
        outputdir_base = 'EAZY_OUTPUT_{}'.format(host190614B.name)
        if run_eazy:
            frbeazy.eazy_input_files(host190614B.photom, os.path.join(eazy_folder, 'inputs'),
                                     host190614B.name,
                                     '../'+outputdir_base,
                                     templates='br07_default',
                                     prior_filter='GMOS_S_r')
            frbeazy.run_eazy(os.path.join(eazy_folder, 'inputs'),
                             host190614B.name,
                             os.path.join(eazy_folder, outputdir_base, 'logfile'))
            # Rename
            os.system('cp -rp {:s} {:s}'.format(os.path.join(eazy_folder, outputdir_base),
                                                pub_path))
        # Redshift
        zgrid, pzi, prior = frbeazy.getEazyPz(-1, MAIN_OUTPUT_FILE='photz',
                                              OUTPUT_DIRECTORY=os.path.join(pub_path, outputdir_base),
                                              CACHE_FILE='Same', binaries=None, get_prior=True)
        zphot, sig_zphot = frbeazy.eazy_stats(zgrid, pzi)

        host190614B.set_z(zphot, 'phot', err=sig_zphot)

        # CIGALE _ build now then run by hand after
        cigale_file = os.path.join(db_path, 'Realfast', 'Law2020', 'G190614_B_CIGALE.fits')
        if build_cigale:
            # Prep
            #embed(header='see 190614 to update this')
            cigale_tbl = photom.copy()

            cigale_tbl = Table(cigale_tbl[1])
            cigale_tbl['z'] = host190614B.z
            cigale_tbl['ID'] = 'G190614_B'

            # Run
            cigale.run(cigale_tbl, 'z', outdir='G190614_B', compare_obs_model=True)
            # Rename/move
            os.system('mv G190614_B/results.fits {:s}'.format(cigale_file))
            model_file = cigale_file.replace('CIGALE', 'CIGALE_model')
            os.system('mv G190614_B/{:s}_best_model.fits {:s}'.format('G190614_B', model_file))
            photo_file = cigale_file.replace('CIGALE.fits', 'CIGALE_photo.dat')
            os.system('mv G190614_B/photo_observed_model_G190614_B.dat {:s}'.format(photo_file))
        # Parse
        host190614B.parse_cigale(cigale_file)

        # Vet all
        assert host190614B.vet_all()

        # Write
        path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
        host190614B.write_to_json(path=path, outfile='G190614_B.json')



    #########################################################
    # C
    #########################################################
    if build_C:
        print("Building Host galaxy C for FRB{}".format(frbname))
        gal_coord = SkyCoord(ra=65.07440, dec=73.70680,  # DES Survey on 06-Jan-2020
                             unit='deg')

        # Instantiate
        host190614C = frbgalaxy.FRBHost(gal_coord.ra.value, 
                                        gal_coord.dec.value, frb190614)

        # Redshift -- JXP measured from FORS2
        #    Should be refined
        host190614C.set_z(1.041, 'phot')

        # Photometry
        EBV = nebular.get_ebv(gal_coord)['meanValue']  # 0.1401

        # Grab the table (requires internet)
        photom_file = os.path.join(db_path, 'Realfast', 'Law2020', 'law2020_photom.ascii')
        if build_photom:
            photom = Table()
            photom['Name'] = ['G{}_C'.format(frbname)]
            photom['ra'] = host190614C.coord.ra.value
            photom['dec'] = host190614C.coord.dec.value
            # These are observed
            photom['GMOS_S_r'] = 25.1
            photom['GMOS_S_r_err'] = 0.3
            photom['LRISb_V'] = 26.1
            photom['LRISb_V_err'] = 0.3
            photom['LRISr_I'] = 25.4
            photom['LRISr_I_err'] = 0.3

            # Dust correct
            for key in photom.keys():
                if key in ['Name', 'ra', 'dec'] or 'err' in key:
                    continue
                # Grab the correction
                if 'LRIS' in key:
                    filt = 'LRIS_'+key[-1]
                else:
                    filt = key
                dust_correct = frbphotom.extinction_correction(filt, EBV)
                mag_dust = 2.5*np.log10(1./dust_correct)
                photom[key] += mag_dust

            # Merge/write
            photom = frbphotom.merge_photom_tables(photom, photom_file)
            # Write
            photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
            print("Wrote photometry to: {}".format(photom_file))

        # Parse
        photom = Table.read(photom_file, format=frbphotom.table_format)
        host190614C.parse_photom(photom, EBV=EBV)

        # CIGALE - Run by hand

        # Vet all
        assert host190614C.vet_all()

        # Write
        path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
        host190614C.write_to_json(path=path, outfile='FRB190614_C.json')


def build_host_190711(build_ppxf=False, build_photom=False, build_cigale=False):
    """ Build the host galaxy data for FRB 191001

    Heintz+2020

    Args:
        build_photom (bool, optional):
        build_cigale (bool, optional):
    """
    frbname = '190711'
    print("Building Host galaxy for FRB{}".format(frbname))
    gal_coord = SkyCoord('J215740.60-802129.25',   # Kasper on Slack 2020 Jan 27; Gemini r-band image
                         unit=(units.hourangle, units.deg))
    EBV = nebular.get_ebv(gal_coord)['meanValue']  #
    print("EBV={} for the host of {}".format(EBV, frbname))

    # Instantiate
    frb190711 = FRB.by_name('FRB190711')
    host190711 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb190711)

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, host190711)
    '''
    # Load redshift table
    ztbl = Table.read(os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'z_SDSS.ascii'),
                      format='ascii.fixed_width')
    z_coord = SkyCoord(ra=ztbl['RA'], dec=ztbl['DEC'], unit='deg')
    idx, d2d, _ = match_coordinates_sky(gal_coord, z_coord, nthneighbor=1)
    if np.min(d2d) > 0.5*units.arcsec:
        embed(header='190608')
    '''
    # Redshift -- JXP estimated from X-Shooter;  should be improved
    host190711.set_z(0.522, 'spec')

    # Morphology
    #host190711.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
    #                                     'HG190711_galfit_GMOS-S_i.log'), 0.16)

    # Photometry

    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'heintz2020_photom.ascii')
    if build_photom:
        # Gemini/GMOS
        photom = Table()
        photom['Name'] = ['HG{}'.format(frbname)]
        photom['ra'] = host190711.coord.ra.value
        photom['dec'] = host190711.coord.dec.value
        photom['GMOS_S_g'] = 24.0
        photom['GMOS_S_g_err'] = 0.2
        photom['GMOS_S_r'] = 23.85
        photom['GMOS_S_r_err'] = 0.15
        photom['GMOS_S_i'] = 23.20
        photom['GMOS_S_i_err'] = 0.15

        # HST
        #photom['WFC3_F160W'] = 22.877 -- USE THIS ONCE THE TRANSMISSION CURVE IS ADDED
        photom['WFC3_F160W'] = 22.877 - 0.072 # Dust corrected
        photom['WFC3_F160W_err'] = 0.012

        # Write
        # Merge/write
        photom = frbphotom.merge_photom_tables(photom, photom_file)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
        print("Wrote photometry to: {}".format(photom_file))
    # Load
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correct
    frbphotom.correct_photom_table(photom, EBV, 'HG190711', required=True)
    # Parse
    host190711.parse_photom(photom, EBV=EBV)

    # CIGALE
    cigale_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG190711_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        cigale.host_run(host190711, cigale_file=cigale_file)
        # Parse
    host190711.parse_cigale(cigale_file, sfh_file=sfh_file)

    # Galfit
    host190711.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
                                   'HG190711_GMOS_i_galfit.fits'))

    '''
    # PPXF
    ppxf_results_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG191001_GMOS_ppxf.ecsv')
    if run_ppxf:
        meta, spectrum = host191001.get_metaspec(instr='GMOS-S')
        spec_fit = None
        ppxf.run(spectrum, 2000., host191001.z, results_file=ppxf_results_file, spec_fit=spec_fit,
                 atmos=[[7150., 7300.], [7580, 7750.]],
                 gaps=[[6675., 6725.]], chk=True)
    host191001.parse_ppxf(ppxf_results_file)

    # Derived quantities

    # AV
    host191001.calc_nebular_AV('Ha/Hb')
    '''

    # Nebular flux measured by a hand (Gaussian fit) by JXP on 2020-05-19
    #   Corrected for Galactic extinction but not internal
    neb_lines = {}
    neb_lines['Hbeta'] = 2.575e-17
    neb_lines['Hbeta_err'] = 5.34e-18

    host190711.neb_lines = neb_lines

    # Galfit -- Mannings+2021
    host190711.parse_galfit(os.path.join(db_path, 'F4', 'mannings2020',
                                   'HG190711_galfit.fits'))

    # SFR
    host190711.calc_nebular_SFR('Hb')

    # Vet all
    assert host190711.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host190711.write_to_json(path=path)


def build_host_190714(build_ppxf=False, build_photom=False, build_cigale=False):
    """ Build the host galaxy data for FRB 190714

    Heintz+2020

    Args:
        build_photom (bool, optional):
        build_cigale (bool, optional):
        build_ppxf (bool, optional):
    """
    frbname = '190714'
    print("Building Host galaxy for FRB{}".format(frbname))
    gal_coord = SkyCoord(ra=183.97955878, dec=-13.02111222, unit='deg')  # Pan-STARRS ID 92371839795415103
                     # J121555.0941-130116.004

    EBV = nebular.get_ebv(gal_coord)['meanValue']  # 0.061

    # Instantiate
    frb190714 = FRB.by_name('FRB190714')
    host190714 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb190714)

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, host190714)

    # Load redshift table
    ztbl = Table.read(os.path.join(db_path, 'CRAFT', 'Heintz2020', 'z_hand.ascii'),
                      format='ascii.fixed_width')
    z_coord = SkyCoord(ztbl['JCOORD'], unit=(units.hourangle, units.deg))  # from DES
    idx, d2d, _ = match_coordinates_sky(gal_coord, z_coord, nthneighbor=1)
    if np.min(d2d) > 0.5*units.arcsec:
        embed(header='190714')

    # Redshift -- JXP measured from LRISr
    host190714.set_z(ztbl['ZEM'][idx], 'spec')

    # Photometry
    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'heintz2020_photom.ascii')
    if build_photom:
        # Pan_STARRS
        search_r = 1 * units.arcsec
        ps_srvy = panstarrs.Pan_STARRS_Survey(gal_coord, search_r)
        ps_tbl = ps_srvy.get_catalog(print_query=True)
        assert len(ps_tbl) == 1
        ps_tbl['Name'] = 'HG{}'.format(frbname)
        # VISTA -- Grabbed from ESO Catalogs
        vista = Table.read(os.path.join(db_path, 'CRAFT', 'Heintz2020', 'FRB190714_eso_vista_photom.csv'), format='ascii')
        assert len(vista) == 1
        assert np.abs(vista['RA2000'][0]-ps_tbl['ra'][0]) < 1e-4
        vtbl = Table()
        vtbl['ra'] = vista['RA2000']
        vtbl['dec'] = vista['DEC2000']
        vtbl['Name'] = ['HG{}'.format(frbname)]
        vtbl['VISTA_Y'] = vista['YAPERMAG6']
        vtbl['VISTA_J'] = vista['JAPERMAG6']
        vtbl['VISTA_H'] = vista['HAPERMAG6']
        vtbl['VISTA_Ks'] = vista['KSAPERMAG6']
        vtbl['VISTA_Y_err'] = vista['YAPERMAG6ERR']
        vtbl['VISTA_J_err'] = vista['JAPERMAG6ERR']
        vtbl['VISTA_H_err'] = vista['HAPERMAG6ERR']
        vtbl['VISTA_Ks_err'] = vista['KSAPERMAG6ERR']
        # HST
        vtbl['WFC3_F160W'] = 18.911
        vtbl['WFC3_F160W_err'] = 0.002
        # VLT
        vtbl['VLT_FORS2_g'] = 20.70
        vtbl['VLT_FORS2_g_err'] = 0.1
        vtbl['VLT_FORS2_I'] = 19.60
        vtbl['VLT_FORS2_I_err'] = 0.1
        # Merge/write
        photom = frbphotom.merge_photom_tables(ps_tbl, vtbl)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
        print("Wrote photometry to: {}".format(photom_file))
    # Parse
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correct
    frbphotom.correct_photom_table(photom, EBV, 'HG190714')
    # Parse
    host190714.parse_photom(photom, EBV=EBV)

    # CIGALE
    cigale_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG190714_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        # Prep
        cut_photom = Table()
        for key in host190714.photom.keys():
            # Let's stick to Pan-STARRS and VISTA only and HST
            if 'Pan-STARRS' not in key and 'VISTA' not in key and 'WFC3' not in key:
                continue
            # Removing Pan-STARRS_y as it is a bit dodgy
            if 'Pan-STARRS_y' in key:
                continue
            cut_photom[key] = [host190714.photom[key]]
        cigale.host_run(host190714, cut_photom=cut_photom, cigale_file=cigale_file)
    # Parse
    host190714.parse_cigale(cigale_file, sfh_file=sfh_file)

    # PPXF
    ppxf_results_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG190714_LRISr_ppxf.ecsv')
    if build_ppxf:
        meta, spectrum = host190714.get_metaspec(instr='LRISr')
        spec_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG190714_LRISr_ppxf.fits')
        R = 1280.
        ppxf.run(spectrum, R, host190714.z, results_file=ppxf_results_file,
                 spec_fit=spec_file,
                 atmos=[(5450., 5600.), (7190., 7380.), (7580, 7700.), (8460., 9000.)], chk=True)
    host190714.parse_ppxf(ppxf_results_file)

    # Derived quantities

    # Galfit
    #host190714.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
    #                               'HG190714_VLT_i_galfit.fits'))
    # Galfit -- Mannings+2021
    host190714.parse_galfit(os.path.join(db_path, 'F4', 'mannings2020',
                                   'HG190714_galfit.fits'))
    # AV
    host190714.calc_nebular_AV('Ha/Hb')

    # SFR
    host190714.calc_nebular_SFR('Ha')

    # Vet all
    assert host190714.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host190714.write_to_json(path=path)


def build_host_191001(build_ppxf=False, build_photom=False, build_cigale=False):
    """ Build the host galaxy data for FRB 191001

    Heintz+2020

    Args:
        build_photom (bool, optional):
        build_cigale (bool, optional):
        run_ppxf (bool, optional):
    """
    frbname = '191001'
    print("Building Host galaxy for FRB{}".format(frbname))
    gal_coord = SkyCoord(ra=323.351851, dec=-54.748515,  # aka J213324.44-544454.65 DES Survey on 06-Jan-2020
                         unit='deg')

    # Instantiate
    frb191001 = FRB.by_name('FRB191001')
    host191001 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb191001)

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, host191001)

    '''
    # Load redshift table
    ztbl = Table.read(os.path.join(db_path, 'CRAFT', 'Bhandari2019', 'z_SDSS.ascii'),
                      format='ascii.fixed_width')
    z_coord = SkyCoord(ra=ztbl['RA'], dec=ztbl['DEC'], unit='deg')
    idx, d2d, _ = match_coordinates_sky(gal_coord, z_coord, nthneighbor=1)
    if np.min(d2d) > 0.5*units.arcsec:
        embed(header='190608')
    '''
    # Redshift -- JXP measured from FORS2
    #    Should be refined
    host191001.set_z(0.2340, 'spec')

    # Morphology
    #host191001.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
    #                                     'HG191001_galfit_FORS2_I.log'), 0.252)

    # Photometry

    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'heintz2020_photom.ascii')
    if build_photom:
        # DES
        search_r = 1 * units.arcsec
        des_srvy = des.DES_Survey(gal_coord, search_r)
        des_tbl = des_srvy.get_catalog(print_query=True)
        host191001.parse_photom(des_tbl)
        # VLT
        photom = Table()
        photom['Name'] = ['HG{}'.format(frbname)]
        photom['ra'] = host191001.coord.ra.value
        photom['dec'] = host191001.coord.dec.value
        photom['VLT_FORS2_g'] = 19.00
        photom['VLT_FORS2_g_err'] = 0.1
        photom['VLT_FORS2_I'] = 17.89
        photom['VLT_FORS2_I_err'] = 0.1
        # Add in DES
        for key in host191001.photom.keys():
            photom[key] = host191001.photom[key]
        # Write
        # Merge/write
        photom = frbphotom.merge_photom_tables(photom, photom_file)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
        print("Wrote photometry to: {}".format(photom_file))
    # Load
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correct
    EBV = nebular.get_ebv(gal_coord)['meanValue']  # 0.061
    frbphotom.correct_photom_table(photom, EBV, 'HG191001')
    # Parse
    host191001.parse_photom(photom, EBV=EBV)

    # CIGALE
    cigale_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG191001_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        # Prep
        cut_photom = Table()
        # Let's stick to DES only
        for key in host191001.photom.keys():
            if 'DES' not in key:
                continue
            cut_photom[key] = [host191001.photom[key]]
        # Run
        cigale.host_run(host191001, cut_photom=cut_photom, cigale_file=cigale_file)
        # Parse
    host191001.parse_cigale(cigale_file, sfh_file=sfh_file)

    # PPXF
    ppxf_results_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG191001_GMOS_ppxf.ecsv')
    spec_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG191001_GMOS_ppxf.fits')
    if build_ppxf:
        meta, spectrum = host191001.get_metaspec(instr='GMOS-S')
        R = meta['R']
        ppxf.run(spectrum, R, host191001.z, results_file=ppxf_results_file, spec_fit=spec_file,
                 atmos=[[7150., 7300.], [7580, 7750.]],
                 gaps=[[6675., 6725.]], chk=True)
    host191001.parse_ppxf(ppxf_results_file)

    # Derived quantities

    # AV
    host191001.calc_nebular_AV('Ha/Hb')

    # SFR
    host191001.calc_nebular_SFR('Ha')

    # Galfit
    #host191001.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
    #                                     'HG191001_VLT_i_galfit.fits'))
    # Galfit -- Mannings+2021
    host191001.parse_galfit(os.path.join(db_path, 'F4', 'mannings2020',
                                   'HG191001_galfit.fits'))

    # Vet all
    assert host191001.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host191001.write_to_json(path=path)


def build_host_200430(build_ppxf=False, build_photom=False, build_cigale=False, run_eazy=False):
    """ Build the host galaxy data for FRB 200430

    Args:
        build_photom (bool, optional):
        build_cigale (bool, optional):
        run_ppxf (bool, optional):
    """
    eazy_folder = './eazy'
    frbname = '200430'
    print("Building Host galaxy for FRB{}".format(frbname))
    gal_coord = SkyCoord('J151849.52+122235.8', unit=(units.hourangle, units.deg)) # SDSS

    EBV = nebular.get_ebv(gal_coord)['meanValue']  # 0.061

    # Instantiate
    frb200430 = FRB.by_name('FRB200430')
    host200430 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frb200430)

    # Load redshift table
    #ztbl = Table.read(os.path.join(db_path, 'CRAFT', 'Heintz2020', 'z_hand.ascii'),
    #                  format='ascii.fixed_width')
    #z_coord = SkyCoord(ztbl['JCOORD'], unit=(units.hourangle, units.deg))  # from DES
    #idx, d2d, _ = match_coordinates_sky(gal_coord, z_coord, nthneighbor=1)
    #if np.min(d2d) > 0.5*units.arcsec:
    #    embed(header='190714')

    # Redshift -- JXP measured from NOT
    host200430.set_z(0.16, 'spec')

    # Photometry
    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'heintz2020_photom.ascii')
    if build_photom:
        # Pan_STARRS
        search_r = 1 * units.arcsec
        ps_srvy = panstarrs.Pan_STARRS_Survey(gal_coord, search_r)
        ps_tbl = ps_srvy.get_catalog(print_query=True)
        assert len(ps_tbl) == 1
        ps_tbl['Name'] = host200430.name
        # WISE?
        # Merge/write
        photom = frbphotom.merge_photom_tables(ps_tbl, photom_file)
        # Write
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
        print("Wrote photometry to: {}".format(photom_file))
    # Load
    photom = Table.read(photom_file, format=frbphotom.table_format)
    # Dust correction
    frbphotom.correct_photom_table(photom, EBV, 'HG200430')
    # Parse
    host200430.parse_photom(photom, EBV=EBV)

    # EAZY
    pub_path = os.path.join(db_path, 'CRAFT', 'Heintz2020')
    outputdir_base = 'EAZY_OUTPUT_{}'.format(host200430.name)
    if run_eazy:
        frbeazy.eazy_input_files(host200430.photom, os.path.join(eazy_folder, 'inputs'),
                                 host200430.name,
                                 '../' + outputdir_base,
                                 #templates='br07_default',
                                 prior_filter='Pan-STARRS_r')
        frbeazy.run_eazy(os.path.join(eazy_folder, 'inputs'),
                         host200430.name,
                         os.path.join(eazy_folder, outputdir_base, 'logfile'))
        # Rename
        os.system('cp -rp {:s} {:s}'.format(os.path.join(eazy_folder, outputdir_base),
                                            pub_path))
    '''
    # Redshift
    zgrid, pzi, prior = frbeazy.getEazyPz(-1, MAIN_OUTPUT_FILE='photz',
                                          OUTPUT_DIRECTORY=os.path.join(pub_path, outputdir_base),
                                          CACHE_FILE='Same', binaries=None, get_prior=True)
    zphot, sig_zphot = frbeazy.eazy_stats(zgrid, pzi)
    host200430.set_z(zphot, 'phot')
    '''

    # CIGALE
    cigale_file = os.path.join(db_path, 'CRAFT', 'Heintz2020', 'HG200430_CIGALE.fits')
    sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
    if build_cigale:
        cigale.host_run(host200430, cigale_file=cigale_file)

    host200430.parse_cigale(cigale_file, sfh_file=sfh_file)

    # Derived quantities

    # AV
    #host190714.calc_nebular_AV('Ha/Hb')

    # SFR
    #host190714.calc_nebular_SFR('Ha')
    
    # Galfit
    host200430.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
                                         'HG200430_SDSS_i_galfit.fits'))

    # Vet all
    assert host200430.vet_all()

    # Write 
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host200430.write_to_json(path=path)


def main(inflg='all', options=None):
    # Options
    build_photom, build_cigale, build_ppxf = False, False, False
    if options is not None:
        if 'photom' in options:
            build_photom = True
        if 'cigale' in options:
            build_cigale = True
        if 'ppxf' in options:
            build_ppxf = True

    if inflg == 'all':
        flg = np.sum(np.array( [2**ii for ii in range(25)]))
    else:
        flg = int(inflg)

    # 121102
    if flg & (2**0):
        build_host_121102(build_photom=build_photom, build_cigale=build_cigale)  # 1

    # 180924
    if flg & (2**1):
        build_host_180924(build_photom=build_photom, build_cigale=build_cigale)  # 2

    # 181112
    if flg & (2**2): # 4
        build_host_181112(build_photom=build_photom, build_cigale=build_cigale)  # 4

    # 190523
    if flg & (2**3):  # 8
        build_host_190523(build_photom=build_photom, build_cigale=build_cigale)

    # 190608
    if flg & (2**4):  # 16
        build_host_190608(build_photom=build_photom, build_cigale=build_cigale)

    # 190102
    if flg & (2**5):  # 32
        build_host_190102(build_photom=build_photom, build_cigale=build_cigale, build_ppxf=build_ppxf)

    # 180916
    if flg & (2**6):  # 64
        build_host_180916(build_photom=build_photom, build_cigale=build_cigale)

    # 190611
    if flg & (2**7):  # 128
        build_host_190611(build_photom=build_photom, build_cigale=build_cigale,
                          source='bright')

    # 190614
    if flg & (2**8):  # 256
        build_host_190614(build_photom=build_photom, build_cigale=build_cigale,
                          build_A=True, build_B=True)

    # 190711
    if flg & (2**9):  # 512
        build_host_190711(build_photom=build_photom, build_cigale=build_cigale)

    # 190714
    if flg & (2**10):  # 1024
        build_host_190714(build_photom=build_photom, build_cigale=build_cigale, build_ppxf=build_ppxf)

    # 191001
    if flg & (2**11):  # 2048
        build_host_191001(build_photom=build_photom, build_cigale=build_cigale, build_ppxf=build_ppxf)

    # 200430
    if flg & (2**12):  # 4096
        build_host_200430(build_photom=build_photom, build_cigale=build_cigale, build_ppxf=build_ppxf)


# Command line execution
if __name__ == '__main__':
    main()


