""" Top-level module to build or re-build the JSON files for
FRB host galaxies"""

from pkg_resources import resource_filename
import os
import sys

from IPython import embed

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table

from frb.galaxies import frbgalaxy, defs
from frb.galaxies import photom as frbphotom
from frb.surveys import des
from frb.surveys import panstarrs

db_path = os.getenv('FRB_GDB')
if db_path is None:
    embed(header='You need to set $FRB_GDB')


def build_host_121102(build_photom=False):
    """
    Generate the JSON file for FRB 121102

    Writes to 121102/FRB121102_host.json

    All of the data currently comes from Tendulkar et al. 2017

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

    # Photometry -- Tendulkar 2017
    photom_file = os.path.join(db_path, 'Repeater', 'Tendulkar2017', 'tendulkar2017_photom.ascii')
    if build_photom:
        photom = Table()
        #photom['Name'] = ['HG121102']  DO NOT USE str columns!
        photom['ra'] = [host121102.coord.ra.value]
        photom['dec'] = host121102.coord.dec.value
        #
        photom['GMOS_r'] = 25.1
        photom['GMOS_r_err'] = 0.1
        photom['GMOS_i'] = 23.9
        photom['GMOS_i_err'] = 0.1
        # Write
        photom = frbphotom.merge_photom_tables(photom, photom_file)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
    host121102.parse_photom(Table.read(photom_file, format=frbphotom.table_format))

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

    # Morphology
    host121102.morphology['reff_ang'] = 0.41
    host121102.morphology['reff_ang_err'] = 0.06
    #
    host121102.morphology['n'] = 2.2
    host121102.morphology['n_err'] = 1.5
    #
    host121102.morphology['b/a'] = 0.25
    host121102.morphology['b/a_err'] = 0.13

    # Derived quantities
    host121102.derived['M_r'] = -17.0  # AB; Tendulkar+17
    host121102.derived['M_r_err'] = 0.2  # Estimated by JXP
    host121102.derived['SFR_nebular'] = 0.23  # MSun/yr; Tendulkar+17
    host121102.derived['Mstar'] = 5.5e7  # Msun; Tendulkar+17
    host121102.derived['Mstar_err'] = 1.5e7  # Msun; Tendulkar+17
    host121102.derived['Z_spec'] = -0.16  # Tendulkar+17 on a scale with Solar O/H = 8.86
    host121102.derived['Z_spec_err'] = -999.  # Tendulkar+17

    # Vet
    assert host121102.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/121102')
    host121102.write_to_json(path=path, overwrite=True)


def build_host_180924(build_photom=True):
    """
    Generate the JSON file for FRB 180924
    
    All data are from Bannister et al. 2019
        https://ui.adsabs.harvard.edu/abs/2019Sci...365..565B/abstract

    Writes to 180924/FRB180924_host.json

    Args:
        build_photom (bool, optional): Generate the photometry file in the Galaxy_DB
    """
    gal_coord = SkyCoord('J214425.25-405400.81', unit=(units.hourangle, units.deg))

    # Instantiate
    host = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, '180924')

    # Redshift
    host.set_z(0.3212, 'spec')

    # Morphology
    host.parse_galfit(os.path.join(db_path, 'CRAFT', 'Bannister2019',
                                   'HG180924_DES_galfit.log'), 0.263)

    # Photometry

    # DES
    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'CRAFT', 'Bannister2019', 'bannister2019_photom.ascii')
    if build_photom:
        search_r = 2*units.arcsec
        des_srvy = des.DES_Survey(gal_coord, search_r)
        des_tbl = des_srvy.get_catalog(print_query=True)
        # Write
        des_tbl.write(photom_file, format=frbphotom.table_format)
    else:
        des_tbl = Table.read(photom_file, format=frbphotom.table_format)

    # Parse
    host.parse_photom(des_tbl)

    # PPXF
    host.parse_ppxf(os.path.join(db_path, 'CRAFT', 'Bannister2019', 'HG180924_MUSE_ppxf.ecsv'))

    # Derived quantities

    # AV
    host.calc_nebular_AV('Ha/Hb')

    # SFR
    host.calc_nebular_SFR('Ha')
    host.derived['SFR_nebular_err'] = -999.

    # CIGALE
    host.parse_cigale(os.path.join(db_path, 'CRAFT', 'Bannister2019',
                                   'HG180924_CIGALE.fits'), 0.263)

    # Vet all
    host.vet_all()

    # Write
    path = resource_filename('frb', 'data/Galaxies/180924')
    host.write_to_json(path=path)


def build_host_181112(build_photom=False):
    """ Build the host galaxy data for FRB 181112

    All of the data comes from Prochaska+2019, Science, in press

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

    # DES
    # Grab the table (requires internet)
    search_r = 2 * units.arcsec
    des_srvy = des.DES_Survey(Host_coord, search_r)
    des_tbl = des_srvy.get_catalog(print_query=True)

    host181112.parse_photom(des_tbl)

    # VLT -- Lochlan 2019-05-02
    # VLT -- Lochlan 2019-06-18
    photom_file = os.path.join(db_path, 'CRAFT', 'Prochaska2019', 'prochaska2019_photom.ascii')
    if build_photom:
        photom = Table()
        photom['Name'] = ['HG{}'.format(frbname)]
        photom['ra'] = host181112.coord.ra.value
        photom['dec'] = host181112.coord.dec.value
        photom['VLT_g'] = 22.57
        photom['VLT_g_err'] = 0.04
        photom['VLT_I'] = 21.51
        photom['VLT_I_err'] = 0.04
        # Add in DES
        for key in host181112.photom.keys():
            photom[key] = host181112.photom[key]
        # Merge/write
        photom = frbphotom.merge_photom_tables(photom, photom_file)
        embed(header='246 of hosts')
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
    host181112.parse_photom(Table.read(photom_file, format=frbphotom.table_format))

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
    host181112.parse_cigale(os.path.join(db_path, 'CRAFT', 'Prochaska2019', 'HG181112_CIGALE.fits'))

    # Write
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host181112.write_to_json(path=path)


def build_host_190523(build_photom=False):  #:run_ppxf=False, build_photom=False):
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
    gal_coord = SkyCoord(ra=207.06433, dec=72.470756, unit='deg')

    # Instantiate
    host190523 = frbgalaxy.FRBHost(gal_coord.ra.value, gal_coord.dec.value, frbname)

    # Load redshift table
    host190523.set_z(0.660, 'spec')

    # Morphology

    # Photometry

    # PanStarrs
    # Grab the table (requires internet)
    photom_file = os.path.join(db_path, 'DSA', 'Ravi2019', 'ravi2019_photom.ascii')
    if build_photom:
        search_r = 1 * units.arcsec
        ps_srvy = panstarrs.Pan_STARRS_Survey(gal_coord, search_r)
        ps_tbl = ps_srvy.get_catalog(print_query=True)
        photom = frbphotom.merge_photom_tables(ps_tbl, photom_file)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
    # Parse
    host190523.parse_photom(Table.read(photom_file, format=frbphotom.table_format))

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
    host190523.parse_cigale(os.path.join(db_path, 'DSA', 'Ravi2019',
                                         'S1_190523_CIGALE.fits'))

    # Derived quantities
    host190523.derived['SFR_nebular'] = 1.3
    host190523.derived['SFR_nebular_err'] = -999.

    # Vet all
    host190523.vet_all()

    # Write -- BUT DO NOT ADD TO REPO (YET)
    path = resource_filename('frb', 'data/Galaxies/{}'.format(frbname))
    host190523.write_to_json(path=path)

    
def main(inflg='all', build_photom=False):

    if inflg == 'all':
        flg = np.sum(np.array( [2**ii for ii in range(25)]))
    else:
        flg = int(inflg)

    # 121102
    if flg & (2**0):
        build_host_121102(build_photom=build_photom)

    # 180924
    if flg & (2**1):
        build_host_180924(build_photom=build_photom)

    # 181112
    if flg & (2**2):
        build_host_181112(build_photom=build_photom)

    # 181112
    if flg & (2**3):
        build_host_190523(build_photom=build_photom)


# Command line execution
if __name__ == '__main__':
    pass

