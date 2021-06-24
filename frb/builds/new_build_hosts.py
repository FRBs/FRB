""" Top-level module to build or re-build the JSON files for
FRB host galaxies"""

from pkg_resources import resource_filename
import os
import sys
import warnings

from IPython import embed

import numpy as np

import pandas

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
from frb.galaxies import utils as galaxy_utils
from frb.galaxies import hosts
from frb.surveys import des
from frb.surveys import sdss
from frb.surveys import wise
from frb.surveys import panstarrs
from frb.surveys import catalog_utils
from frb.surveys import survey_utils
from frb import utils
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

def assign_z(ztbl_file:str, host:frbgalaxy.FRBHost):
    # Load redshift table
    ztbl = Table.read(ztbl_file, format='ascii.fixed_width')

    z_coord = SkyCoord(ztbl['JCOORD'], unit=(units.hourangle, units.deg))  
    idx, d2d, _ = match_coordinates_sky(host.coord, z_coord, nthneighbor=1)
    if np.min(d2d) > 0.5*units.arcsec:
        raise ValueError("No matching galaxy!")

    # Set redshift 
    host.set_z(ztbl['ZEM'][idx], 'spec')

def search_for_file(projects, references, root,
                    prefix='ref', return_last_file=False):
    found_file = None
    found = False
    for project, ref in zip(projects, references):
        GDB_path = os.path.join(db_path, project, ref)
        if prefix == 'ref':
            filename = os.path.join(GDB_path, ref.lower()+root)
        else:
            filename = os.path.join(GDB_path, prefix+root)
        # Is it there?
        if os.path.isfile(filename):
            found_file = filename
            found = True
    if not found and return_last_file:
        found_file = filename
    #
    return found, found_file
            


def run(host_input:pandas.core.series.Series, 
        build_ppxf=False, build_photom=False, 
        build_cigale=False, is_host=True):

    frbname = utils.parse_frb_name(host_input.FRB)
    print(f"Building Host galaxy for {frbname}")
    gal_coord = SkyCoord(host_input.Coord, frame='icrs')

    # Instantiate
    Frb = FRB.by_name(frbname)
    Host = frbgalaxy.FRBHost(gal_coord.ra.value, 
                             gal_coord.dec.value, 
                             Frb)
                            
    # File root
    file_root = Host.name if is_host else utils.name_from_coord(Host.coord)

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, Host)

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
    Host.set_z(host_input.z, 'spec')


    # ####################
    # Photometry

    # Grab the table (requires internet)
    project_list = host_input.Projects.split(',')
    ref_list = host_input.References.split(',')
    if len(project_list) > 1:
        embed(header='2422 -- not ready for this yet')
    GDB_path = os.path.join(db_path, project_list[0],
        ref_list[0])
    found_photom, photom_file = search_for_file(
        project_list, ref_list, '_photom.ascii', 
        return_last_file=True)
    #                           ref_list[0].lower()+'_photom.ascii')
    if build_photom:
        search_r = 1 * units.arcsec

        # Find all the surveys
        inside = survey_utils.in_which_survey(Frb.coord)

        # Loop
        merge_tbl = None
        for key in inside.keys():
            # Skip? 
            if key in ['NVSS', 'FIRST', 'WENSS'] or not inside[key]:
                continue
            # Skip WISE?
            if key in ['WISE'] and inside['DES']:
                continue
            # Slurp
            survey = survey_utils.load_survey_by_name(key, 
                                                      gal_coord, 
                                                      search_r)
            srvy_tbl = survey.get_catalog(print_query=True)
            if len(srvy_tbl) == 0:
                continue
            # Merge
            if merge_tbl is None:
                merge_tbl = srvy_tbl
            else:
                merge_tbl = frbphotom.merge_photom_tables(srvy_tbl, merge_tbl)

        # Add name
        merge_tbl['Name'] = file_root

        '''
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
        '''
        # Write
        # Merge with table from disk
        photom = frbphotom.merge_photom_tables(merge_tbl, photom_file)
        photom.write(photom_file, format=frbphotom.table_format, overwrite=True)
        print("Wrote photometry to: {}".format(photom_file))
        found_photom = True

    # Load photometry
    if found_photom:
        photom = Table.read(photom_file, format=frbphotom.table_format)
        # Dust correct
        EBV = nebular.get_ebv(gal_coord)['meanValue']  # 0.061
        frbphotom.correct_photom_table(photom, EBV, Host.name)
        # Parse
        Host.parse_photom(photom, EBV=EBV)
    else:
        print(f"No photom file found {file_root}")
        

    # CIGALE
    cigale_file = os.path.join(GDB_path, 
                               file_root+'_CIGALE.fits')
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
    if os.path.isfile(cigale_file):
        Host.parse_cigale(cigale_file, sfh_file=sfh_file)
    else:
        print(f"No CIGALE file to read for {file_root}")

    # PPXF
    ppxf_results_file = os.path.join(GDB_path, 
                               Host.name+f'_{host_input.Spectrum}_ppxf.ecsv')
    spec_file = ppxf_results_file.replace('ecsv', 'fits')

    if build_ppxf:
        meta, spectrum = host191001.get_metaspec(instr='GMOS-S')
        R = meta['R']
        ppxf.run(spectrum, R, host191001.z, results_file=ppxf_results_file, spec_fit=spec_file,
                 atmos=[[7150., 7300.], [7580, 7750.]],
                 gaps=[[6675., 6725.]], chk=True)
    # Load
    if os.path.isfile(ppxf_results_file):
        Host.parse_ppxf(ppxf_results_file)
    else:
        print(f"No pPXF file to read for {file_root}")

    # Derived quantities

    # AV
    Host.calc_nebular_AV('Ha/Hb')

    # SFR
    Host.calc_nebular_SFR('Ha')

    # Galfit
    warnings.warn("ADD GALFIT BACK!!")
    '''
    Host.parse_galfit(os.path.join(db_path, 'F4', 'mannings2020',
                                   'HG191001_galfit.fits'))
    '''

    # Vet all
    assert Host.vet_all()

    # Write
    out_path = os.path.join(resource_filename('frb', 'data'),
        'Galaxies', f'{frbname[3:]}')
    outfile = None if is_host else \
        utils.name_from_coord(Host.coord) + '.json'
        #utils.name_from_coord(Host.coord) + '_{}.json'.format(frbname)
    Host.write_to_json(path=out_path, outfile=outfile)

def main(frbs, options=None):
    # Options
    build_photom, build_cigale, build_ppxf = False, False, False
    if options is not None:
        if 'photom' in options:
            build_photom = True
        if 'cigale' in options:
            build_cigale = True
        if 'ppxf' in options:
            build_ppxf = True

    # Read public host table
    host_tbl = hosts.load_host_tbl()

    # Loop me
    if frbs == 'all':
        embed(header='Generate code to (i) load up the FRB table; (ii) generate a list')
    elif isinstance(frbs, list):
        pass

    for frb in frbs:
        frb_name = utils.parse_frb_name(frb, prefix='')
        mt_idx = host_tbl.FRB == frb_name
        idx = np.where(mt_idx)[0].tolist()
        # Loop on em
        is_host = True
        # Do it!
        for ii in idx:
            run(host_tbl.iloc[ii], build_photom=build_photom, 
                build_cigale=build_cigale, build_ppxf=build_ppxf,
                is_host=is_host)
            # Any additional ones are treated as candidates
            is_host = False
