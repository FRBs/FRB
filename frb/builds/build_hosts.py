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

from frb.galaxies import frbgalaxy, defs
from frb.galaxies import photom as frbphotom
from frb.galaxies import ppxf
from frb.galaxies import nebular
from frb.surveys import des
from frb.surveys import sdss
from frb.surveys import wise
from frb.surveys import panstarrs
from frb.surveys import catalog_utils

try:
    import extinction
except ImportError:
    print("extinction package not loaded.  Extinction corrections will fail")

try:
    from frb.galaxies import cigale
except ModuleNotFoundError:
    warnings.warn("You haven't installed pcigale and won't be able to do that analysis")

from linetools.spectra.xspectrum1d import XSpectrum1D

db_path = os.getenv('FRB_GDB')
if db_path is None:
    print("Warning, you need to set $FRB_GDB to build hosts")
    #embed(header='You need to set $FRB_GDB')

ebv_method = 'SandF'

def build_host_121102(build_photom=False, build_cigale=False, use_orig=False):
    """
    Generate the JSON file for FRB 121102

    Writes to 121102/FRB121102_host.json

    The majority of data comes from Tendulkar et al. 2017

    reff from Bassa+17 https://ui.adsabs.harvard.edu/abs/2017ApJ...843L...8B/abstract


    Args:
        build_photom (bool, optional): Generate the photometry file in the Galaxy_DB

    """
    FRB_coord = SkyCoord('05h31m58.698s +33d8m52.59s', frame='icrs')
    # Eyeball Tendulkar+17 PA
    gal_coord = FRB_coord.directional_offset_by(-45 * units.deg, 286e-3 * units.arcsec)

    # Instantiate
    host121102 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, '121102')

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

    # Vette
    for key in host121102.neb_lines.keys():
        if '_err' in key:
            continue
        assert key in defs.valid_neb_lines

    # Morphology : Bassa+2017 half-light
    host121102.morphology['reff_ang'] = 0.20   # arcsec
    host121102.morphology['reff_ang_err'] = 0.01
    # Other
    host121102.morphology['n'] = 2.2
    host121102.morphology['n_err'] = 1.5
    #
    host121102.morphology['b/a'] = 0.25
    host121102.morphology['b/a_err'] = 0.13

    # Derived quantities
    if use_orig:
        host121102.derived['M_r'] = -17.0  # AB; Tendulkar+17
        host121102.derived['M_r_err'] = 0.2  # Estimated by JXP
        host121102.derived['SFR_nebular'] = 0.23  # MSun/yr; Tendulkar+17
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
    host = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, '180924')

    # Redshift -- JXP measured from multiple data sources
    host.set_z(0.3212, 'spec')

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
        photom['ra'] = host.coord.ra.value
        photom['dec'] = host.coord.dec.value
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
    host.parse_photom(photom, EBV=EBV)

    # PPXF
    host.parse_ppxf(os.path.join(db_path, 'CRAFT', 'Bannister2019', 'HG180924_MUSE_ppxf.ecsv'))

    # Derived quantities

    # AV
    host.calc_nebular_AV('Ha/Hb')

    # SFR
    host.calc_nebular_SFR('Ha')
    host.derived['SFR_nebular_err'] = -999.

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
    host.parse_cigale(cigale_file, sfh_file=sfh_file)

    # Galfit
    host.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
                                   'HG180924_DES_i_galfit.fits'))

    # Vet all
    host.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/180924')
    host.write_to_json(path=path)


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
    Host_coord = SkyCoord('J214923.66-525815.28',
                          unit=(units.hourangle, units.deg))  # from DES

    # Instantiate
    host181112 = frbgalaxy.FRBHost(Host_coord.ra.value, Host_coord.dec.value, frbname)
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
    host190102 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frbname)

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
        #alAV = nebular.load_extinction('MW')
        #Al = alAV(spectrum.wavelength.value) * AV
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
    host190102.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
                                   'HG190102_VLT_i_galfit.fits'))

    # Vet all
    host190102.vet_all()

    # Write -- BUT DO NOT ADD TO REPO (YET)
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
    host190523_S1 = frbgalaxy.FRBHost(S1_gal_coord.ra.value, S1_gal_coord.dec.value, frbname)
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

    # Write -- BUT DO NOT ADD TO REPO (YET)
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
    host190608 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frbname)

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
    host190608.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
                                   'HG190608_SDSS_i_galfit.fits'))
    # Vet all
    host190608.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host190608.write_to_json(path=path)


def build_host_180916(run_ppxf=False, build_photom=False, build_cigale=False):
    """ Build the host galaxy data for FRB 180916

    All of the data currently comes from Marcote et al. 2020
    https://ui.adsabs.harvard.edu/abs/2020Natur.577..190M/abstract

    CIGALE will be published in Kasper et al. 2020

    Args:
        build_photom (bool, optional):
    """
    frbname = '180916'
    gal_coord = SkyCoord('J015800.28+654253.0',
                         unit=(units.hourangle, units.deg))

    # Instantiate
    host180916 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frbname)

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
    host180916.parse_galfit(os.path.join(db_path, 'CRAFT', 'Heintz2020',
                                   'HG180916_SDSS_i_galfit.fits'))

    # Vet all
    assert host180916.vet_all()

    # Write -- BUT DO NOT ADD TO REPO (YET)
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host180916.write_to_json(path=path)


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
        #build_host_180916(build_photom=build_photom)

# Command line execution
if __name__ == '__main__':
    pass

