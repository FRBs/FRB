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
fill_value = -999.

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
            
def read_lit_table(lit_entry, coord=None):
    """ Reade a literature table

    Args:
        lit_entry (pandas row): Row of the overview table
        coord (astropy.coordiantes.SkyCoord, optional): [description]. Defaults to None.

    Raises:
        ValueError: [description]

    Returns:
        astropy.table.Table: table of literature data
    """
    literature_path = os.path.join(resource_filename('frb', 'data'), 
                                       'Galaxies', 'Literature')
    lit_file = os.path.join(literature_path, lit_entry.Table)
    if lit_entry.Format == 'csv':
        lit_tbl = pandas.read_csv(lit_file)
    else:
        lit_tbl = Table.read(lit_file, format=lit_entry.Format)
    
    # 
    if coord is not None:
        tbl_coord = SkyCoord(ra=lit_tbl['ra'], dec=lit_tbl['dec'], unit='deg')
        sep = coord.separation(tbl_coord)
        match = sep < 1*units.arcsec
        nmatch = np.sum(match)
        if nmatch == 0:
            return None
        elif nmatch == 1:
            idx = int(np.where(match)[0])
            return lit_tbl[idx:idx+1]
        else:
            raise ValueError("More than one match in the table!!!")
    else:
        # Return
        return lit_tbl

def run(host_input:pandas.core.series.Series, 
        build_ppxf=False, build_photom=False, 
        lit_refs=None,
        build_cigale=False, is_host=True):

    frbname = utils.parse_frb_name(host_input.FRB)
    print("--------------------------------------")
    print(f"Building Host galaxy for {frbname}")
    gal_coord = SkyCoord(host_input.Coord, frame='icrs')

    # Instantiate
    Frb = FRB.by_name(frbname)
    Host = frbgalaxy.FRBHost(gal_coord.ra.value, 
                             gal_coord.dec.value, 
                             Frb)

    # File root
    file_root = Host.name if is_host else utils.name_from_coord(Host.coord)
    project_list = host_input.Projects.split(',')
    ref_list = host_input.References.split(',')

    # UPDATE RA, DEC, OFFSETS
    offsets.incorporate_hst(mannings2021_astrom, Host)

    '''
    # Load redshift table
    assign_z()
    '''

    # Redshift 
    warnings.warn("We should be using the z Table..")
    Host.set_z(host_input.z, 'spec')


    # ####################
    # Photometry

    search_r = 1 * units.arcsec

    # Survey data
    inside = survey_utils.in_which_survey(Frb.coord)

    # Loop on em
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
        if srvy_tbl is None or len(srvy_tbl) == 0:
            continue
        elif len(srvy_tbl) > 1:
            raise ValueError("You found more than 1 galaxy.  Uh-oh!")
        warnings.warn("We need a way to reference the survey")
        # Merge
        if merge_tbl is None:
            merge_tbl = srvy_tbl
            merge_tbl['Name'] = file_root
        else:
            embed(header='207 of new')
            merge_tbl = frbphotom.merge_photom_tables(srvy_tbl, merge_tbl)
    
    # Literature time
    if lit_refs is None:
        lit_refs = os.path.join(resource_filename('frb', 'data'), 'Galaxies',
            'Literature', 'all_refs.csv')
    lit_tbls = pandas.read_csv(lit_refs)

    for kk in range(len(lit_tbls)):
        lit_entry = lit_tbls.iloc[kk]
        if 'photom' not in lit_entry.Table:
            continue
        # Load table
        sub_tbl = read_lit_table(lit_entry, coord=Host.coord)
        if sub_tbl is not None:
            # Add Ref
            for key in sub_tbl.keys():
                if 'err' in key:
                    newkey = key.replace('err', 'ref')
                    sub_tbl[newkey] = lit_entry.Reference
            # Merge?
            if merge_tbl is not None:
                for key in sub_tbl.keys():
                    if key in merge_tbl.keys():
                        if sub_tbl[key] == fill_value:
                            continue
                        else:
                            merge_tbl[key] = sub_tbl[key]
                    else:
                        if sub_tbl[key] != fill_value:
                            merge_tbl[key] = sub_tbl[key]
                merge_tbl = frbphotom.merge_photom_tables(sub_tbl, merge_tbl)
            else:
                merge_tbl = sub_tbl
                merge_tbl['Name'] = file_root

    # Finish
    if merge_tbl is not None:
        # Dust correct
        EBV = nebular.get_ebv(gal_coord)['meanValue']  # 0.061
        frbphotom.correct_photom_table(merge_tbl, EBV, Host.name)
        # Parse
        Host.parse_photom(merge_tbl, EBV=EBV)
    else:
        print(f"No photometry for {file_root}")
        

    # CIGALE
    found_cigale, cigale_file = search_for_file(
        project_list, ref_list, '_CIGALE.fits',
        prefix=file_root,
        return_last_file=build_cigale)
    #cigale_file = os.path.join(GDB_path, 
    #                           file_root+'_CIGALE.fits')
    if cigale_file is not None:
        sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')

    if build_cigale:
        embed(header='251 -- not ready for this!')
        # Prep
        cut_photom = Table()
        # Cut me!
        for key in Host.photom.keys():
            cut_photom[key] = [host191001.photom[key]]
        # Run
        cigale.host_run(host191001, cut_photom=cut_photom, cigale_file=cigale_file)
        found_cigale = True

    # Parse
    if found_cigale:
        print(f"Slupring in CIGALE outputs from {cigale_file}")
        Host.parse_cigale(cigale_file, sfh_file=sfh_file)
    else:
        print(f"No CIGALE file to read for {file_root}")
    

    # PPXF
    found_ppxf, ppxf_results_file = search_for_file(
        project_list, ref_list, '_ppxf.ecsv',
        prefix=file_root+f'_{host_input.Spectrum}',
        return_last_file=build_ppxf)
    #ppxf_results_file = os.path.join(GDB_path, 
    #                           Host.name+f'_{host_input.Spectrum}_ppxf.ecsv')
    if ppxf_results_file is not None:
        spec_file = ppxf_results_file.replace('ecsv', 'fits')

    if build_ppxf:
        meta, spectrum = Host.get_metaspec(instr=host_input.Spectrum)
        R = meta['R']
        gaps_str = host_input.ppxf_cuts.split(';')
        gaps = []
        for gap in gaps_str:
            gaps.append([float(item) for item in gap.split(',')])

        # Correct for Galactic extinction
        ebv = float(nebular.get_ebv(Host.coord)['meanValue'])
        print(f'Correcting the spectrum for Galactic extinction with reddening E(B-V)={ebv}')
        AV = ebv * 3.1  # RV
        Al = extinction.ccm89(spectrum.wavelength.value, AV, 3.1)
        # New spec
        new_flux = spectrum.flux * 10**(Al/2.5)
        new_sig = spectrum.sig * 10**(Al/2.5)
        new_spec = XSpectrum1D.from_tuple((spectrum.wavelength, new_flux, new_sig))

        #
        embed(header='300 of new')
        ppxf.run(new_spec, R, host_input.z, 
                 results_file=ppxf_results_file, 
                 spec_fit=spec_file,
                 #atmos=[[7150., 7300.], [7580, 7750.]],
                 gaps=gaps, chk=True)
        found_ppxf = True
    # Load
    if found_ppxf:
        print(f"Slupring in pPXF outputs from {ppxf_results_file}")
        Host.parse_ppxf(ppxf_results_file)
    else:
        print(f"No pPXF file to read for {file_root}")

    # Slurp in literature Nebular
    for kk in range(len(lit_tbls)):
        lit_entry = lit_tbls.iloc[kk]
        if 'nebular' not in lit_entry.Table:
            continue
        # Load table
        lit_tbl = read_lit_table(lit_entry, coord=Host.coord)
        if lit_tbl is None:
            continue
        # Fill me in 
        for key in lit_tbl.keys():
            if 'err' in key:
                Host.neb_lines[key] = float(lit_tbl[key].values[0])
                newkey = key.replace('err', 'ref')
                Host.neb_lines[newkey] = lit_entry.Reference
                # Value
                newkey = newkey.replace('_ref', '')
                Host.neb_lines[newkey] = float(lit_tbl[newkey].values[0])

    # AV
    if 'Halpha' in Host.neb_lines.keys() and 'Hbeta' in Host.neb_lines.keys():
        Host.calc_nebular_AV('Ha/Hb')

    # SFR
    if 'Halpha' in Host.neb_lines.keys():
        Host.calc_nebular_SFR('Ha')
    elif 'Hbeta' in Host.neb_lines.keys():
        Host.calc_nebular_SFR('Hb')

    # Galfit
    if isinstance(host_input.Galfit_filter, str):
        found_galfit, galfit_file = search_for_file(
            project_list, ref_list, '_galfit.fits',
            prefix=file_root+'_'+host_input.Galfit_filter)
        if found_galfit:
            Host.parse_galfit(galfit_file)
        else:
            print(f"Galfit file {galfit_file} not found!")
    else:
        print("Galfit analysis not enabled")

    # Vet all
    assert Host.vet_all()

    # Write
    out_path = os.path.join(resource_filename('frb', 'data'),
        'Galaxies', f'{frbname[3:]}')
    outfile = None if is_host else \
        utils.name_from_coord(Host.coord) + '.json'
        #utils.name_from_coord(Host.coord) + '_{}.json'.format(frbname)
    Host.write_to_json(path=out_path, outfile=outfile)

def main(frbs, options=None, hosts_file=None, lit_refs=None):
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
    host_tbl = hosts.load_host_tbl(hosts_file=hosts_file)

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
                is_host=is_host, lit_refs=lit_refs)
            # Any additional ones are treated as candidates
            is_host = False

    # 
    print("All done!")