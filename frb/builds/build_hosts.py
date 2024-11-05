""" Top-level module to build or re-build the JSON files for
FRB host galaxies"""

import importlib_resources
from importlib.resources import files
import os
import sys
import warnings

from IPython import embed

import numpy as np
import requests

import pandas

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table
from astropy.coordinates import match_coordinates_sky

from frb.frb import FRB
from frb.galaxies import frbgalaxy, offsets
from frb.galaxies import photom as frbphotom
try:
    from frb.galaxies import ppxf
except:
    print('WARNING:  ppxf not installed')
from frb.galaxies import nebular
from frb.galaxies import utils as galaxy_utils
from frb.galaxies import hosts
from frb.galaxies import defs as galaxy_defs
from frb.surveys import survey_utils
from frb import utils
import pandas

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
mannings2021_astrom = pandas.read_csv(importlib_resources.files('frb.data.Galaxies.Additional.Mannings2021')/'astrometry_v2.csv')
# Probably will rename this                                        
mannings2021_astrom = mannings2021_astrom[
    (mannings2021_astrom.Filter == 'F160W') | (
        mannings2021_astrom.Filter == 'F110W')].copy()

def chk_fill(value):
    # Masked?
    if isinstance(value,str):
        # Allow for -999 as a string
        try:
            value = float(value)
        except:
            return False
        else:
            return np.isclose(value, fill_value)
    elif isinstance(value,np.ma.core.MaskedConstant):
        return False
    else:
        return np.isclose(value, fill_value)

def assign_z(ztbl_file:str, host:frbgalaxy.FRBHost):
    """Assign a redshift using one of the Galaxy_DB tables

    Args:
        ztbl_file (str): table file
        host (frbgalaxy.FRBHost): host object

    Raises:
        ValueError: [description]
    """
    # Load redshift table
    ztbl = Table.read(ztbl_file, format='ascii.fixed_width')

    z_coord = SkyCoord(ztbl['JCOORD'], unit=(units.hourangle, units.deg))  
    idx, d2d, _ = match_coordinates_sky(host.coord, z_coord, nthneighbor=1)
    if np.min(d2d) > 0.5*units.arcsec:
        raise ValueError("No matching galaxy!")

    # Set redshift 
    host.set_z(ztbl['ZEM'][idx], 'spec')

def search_for_file(projects, references, root:str,
                    prefix='ref', return_last_file=False):
    """ Search for a given data file

    If multiple files are found, the *last* one is returned

    Args:
        projects ([type]): [description]
        references ([type]): [description]
        root (str): [description]
        prefix (str, optional): [description]. Defaults to 'ref'.
        return_last_file (bool, optional): [description]. Defaults to False.

    Returns:
        tuple: bool, str  [file was found?, name of file with path]
    """
    
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
        coord (astropy.coordiantes.SkyCoord, optional): Coordinate
            for the galaxy of interest. Defaults to None.

    Raises:
        ValueError: [description]

    Returns:
        astropy.table.Table: table of literature data
    """
    literature_path = importlib_resources.files('frb.data.Galaxies.Literature')
    lit_file = literature_path/lit_entry.Table
    if lit_entry.Format == 'csv':
        lit_tbl = Table.from_pandas(pandas.read_csv(lit_file))
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
        build_ppxf:bool=False, 
        lit_refs:str=None,
        skip_surveys:bool=False,
        build_cigale:bool=False, is_host:bool=True,
        override:bool=False, out_path:str=None,
        outfile:str=None):
    """Main method for generating a Host JSON file

    Args:
        host_input (pandas.core.series.Series): Row of the CVS file
            providing the host inputs
        build_ppxf (bool, optional): Run pPXF?. Defaults to False.
        lit_refs (str, optional): File of literature references. Defaults to None.
        build_cigale (bool, optional): Run CIGALE?. Defaults to False.
            NOT IMPLEMENTED (yet, but maybe never)
        is_host (bool, optional): Object is a Host, as opposed
            to a neighboring/foreground galaxy. Defaults to True.
        override (bool, optional): Attempt to over-ride errors. 
            Mainly for time-outs of public data. Defaults to False.
        outfile (str, optional): Over-ride default outfile [not recommended; mainly for testing]
        out_path (str, optional): Over-ride default outfile [not recommended; mainly for testing]
        skip_surveys (bool, optional): 
            Skip the survey data.  Useful for testing. Defaults to False.


    Raises:
        e: [description]
        ValueError: [description]
    """

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
    project_list = host_input.Projects.split(',') if isinstance(
        host_input.Projects,str) else []
    ref_list = host_input.References.split(',') if isinstance(
        host_input.References,str) else []
    assert len(project_list) == len(ref_list)

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
    if not skip_surveys:
        try:
            inside = survey_utils.in_which_survey(Frb.coord)
        except (requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout) as e:  # Catches time-out from survey issues
            if override:
                print("Survey timed out.  You should re-run it sometime...")
                inside = {}
            else:
                raise e
    else:
        inside = {}

    # Loop on em
    merge_tbl = None
    for key in inside.keys():
        # Skip? 
        if key in ['NVSS', 'FIRST', 'WENSS'] or not inside[key]:
            continue

        # Slurp
        survey = survey_utils.load_survey_by_name(key, 
                                                    gal_coord, 
                                                    search_r)
        srvy_tbl = survey.get_catalog(print_query=True)

        if srvy_tbl is None or len(srvy_tbl) == 0:
            continue
        elif len(srvy_tbl) > 1:
            if override:
                warnings.warn(f"There was more than 1 galaxy for this FRB for survey: {key}")
                print("Proceeding without using this survey")
                continue
            else:
                #print("You found more than 1 galaxy.  Taking the 2nd one")
                #srvy_tbl = srvy_tbl[1:]
                #srvy_tbl = srvy_tbl[:1]
                raise ValueError("You found more than 1 galaxy.  Uh-oh!")
        warnings.warn("We need a way to reference the survey")
        # Merge
        if merge_tbl is None:
            merge_tbl = srvy_tbl
            merge_tbl['Name'] = file_root
        else:
            merge_tbl = frbphotom.merge_photom_tables(srvy_tbl, merge_tbl)

    # Literature time
    if lit_refs is None:
        lit_refs = importlib_resources.files('frb.data.Galaxies.Literature')/'all_refs.csv'
    lit_tbls = pandas.read_csv(lit_refs, comment='#')

    for kk in range(len(lit_tbls)):
        lit_entry = lit_tbls.iloc[kk]
        if 'photom' not in lit_entry.Table:
            continue
        # Load table
        sub_tbl = read_lit_table(lit_entry, coord=Host.coord)
        if sub_tbl is not None:
            # Add References, unless the value is masked
            for key in sub_tbl.keys():
                if 'err' in key and not chk_fill(sub_tbl[key].data[0]):
                    newkey = key.replace('err', 'ref')
                    sub_tbl[newkey] = lit_entry.Reference
            # Merge?
            if merge_tbl is not None:
                for key in sub_tbl.keys():
                    if key == 'Name':
                        continue
                    if key in merge_tbl.keys() and not chk_fill(sub_tbl[key].data[0]):
                        merge_tbl[key] = sub_tbl[key]
                    else:
                        if not chk_fill(sub_tbl[key].data[0]):
                            merge_tbl[key] = sub_tbl[key]
            else:
                merge_tbl = sub_tbl
                merge_tbl['Name'] = file_root

    # Remove NSC for now
    if merge_tbl is not None:
        for key in merge_tbl.keys():
            if 'NSC' in key:
                merge_tbl.remove_column(key)
                print(f"Removing NSC column: {key}")
    # Finish
    embed(header='324 of build')
    if merge_tbl is not None:
        # Dust correct
        EBV = nebular.get_ebv(gal_coord)['meanValue']
        code = frbphotom.correct_photom_table(merge_tbl, EBV, Host.name)

        if code == -1:
            raise ValueError("Bad extinction correction!")
        # Parse
        Host.parse_photom(merge_tbl, EBV=EBV)
    else:
        print(f"No photometry for {file_root}")

    # CIGALE
    found_cigale, cigale_file = search_for_file(
        project_list, ref_list, '_CIGALE.fits',
        prefix=file_root,
        return_last_file=build_cigale)


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
        new_spec = galaxy_utils.deredden_spec(spectrum, ebv)

        # Run it
        ppxf.run(new_spec, R, host_input.z, 
                 results_file=ppxf_results_file, 
                 spec_fit=spec_file,
                 gaps=gaps, chk=True)
        found_ppxf = True
    # Load
    if found_ppxf:
        print(f"Slurping in pPXF outputs from {ppxf_results_file}")
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
                Host.neb_lines[key] = float(lit_tbl[key].data[0])
                newkey = key.replace('err', 'ref')
                Host.neb_lines[newkey] = lit_entry.Reference
                # Value
                newkey = newkey.replace('_ref', '')
                Host.neb_lines[newkey] = float(lit_tbl[newkey].data[0])

    # Remove bad lines
    if isinstance(host_input.Bad_EM_lines, str):
        lines = host_input.Bad_EM_lines.split(',')
        for line in lines:
            Host.neb_lines.pop(line)
            Host.neb_lines.pop(line+'_err')

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
            print(f"Galfit analysis slurped in via: {galfit_file}")
            Host.parse_galfit(galfit_file)
        else:
            embed(header='435 of build')
            raise IOError(f"Galfit file with filter {host_input.Galfit_filter} not found!")
    else:
        print("Galfit analysis not enabled")

    # Derived from literature
    for kk in range(len(lit_tbls)):
        lit_entry = lit_tbls.iloc[kk]
        if 'derived' not in lit_entry.Table:
            continue
        # Load table
        lit_tbl = read_lit_table(lit_entry, coord=Host.coord)
        if lit_tbl is None:
            continue
        # Fill me in 
        for key in lit_tbl.keys():
            # Handle multiple approaches to errors
            for err_type in galaxy_defs.allowed_errors:
                if err_type in key and '_uperr' not in key: 
                    valkey = key.replace(err_type, '')
                    refkey = valkey+'_ref'
                    # Scrub any existing!
                    if valkey in Host.derived.keys():
                        Host.derived.pop(valkey)
                        for err_type in galaxy_defs.allowed_errors:
                            errkey = valkey+err_type
                            if errkey in Host.derived.keys():
                                Host.derived.pop(errkey)
                    if refkey in Host.derived.keys():
                        Host.derived.pop(refkey)
                    # Fill
                    Host.derived[valkey] = float(lit_tbl[valkey].data[0])
                    Host.derived[refkey] = lit_entry.Reference
                    # Errors
                    Host.derived[key] = float(lit_tbl[key].data[0])
                    if '_loerr' in key:
                        hikey = key.replace('lo', 'up')
                        Host.derived[hikey] = float(lit_tbl[hikey].data[0])

    # Vet all
    assert Host.vet_all()

    # Write
    if out_path is None:
        out_path = str(files('frb').joinpath('data', 'Galaxies', frbname[3:]))
    if outfile is None:
        outfile = None if is_host else \
            utils.name_from_coord(Host.coord) + '.json'
        #utils.name_from_coord(Host.coord) + '_{}.json'.format(frbname)
    Host.write_to_json(path=out_path, outfile=outfile)


def main(frbs:list, options:str=None, hosts_file:str=None, lit_refs:str=None,
         override:bool=False, outfile:str=None, out_path:str=None):
    """ Driver of the analysis

    Args:
        frbs (list): [description]
        options (str, optional): [description]. Defaults to None.
        hosts_file (str, optional): [description]. Defaults to None.
        lit_refs (str, optional): [description]. Defaults to None.
        override (bool, optional): [description]. Defaults to False.
        outfile (str, optional): [description]. Defaults to None.
            Here for testing
        out_path (str, optional): [description]. Defaults to None.
            Here for testing
    """
    # Options
    build_cigale, build_ppxf, skip_surveys = False, False, False
    if options is not None:
        if 'cigale' in options:
            build_cigale = True
        if 'ppxf' in options:
            build_ppxf = True
        if 'skip_surveys' in options:
            skip_surveys = True

    # Read public host table
    host_tbl = hosts.load_host_tbl(hosts_file=hosts_file)

    # Loop me
    if frbs == 'all':
        embed(header='Generate code to (i) load up the FRB table; (ii) generate a list')
    elif isinstance(frbs, list):
        pass

    for frb in frbs:
        frb_name = utils.parse_frb_name(frb, prefix='')
        print(f'Working on {frb_name}')
        mt_idx = host_tbl.FRB == frb_name
        idx = np.where(mt_idx)[0].tolist()
        # Loop on em
        is_host = True
        # Do it!
        for ii in idx:
            run(host_tbl.iloc[ii], 
                build_cigale=build_cigale, build_ppxf=build_ppxf,
                skip_surveys=skip_surveys,
                is_host=is_host, lit_refs=lit_refs, override=override,
                outfile=outfile, out_path=out_path)
            # Any additional ones are treated as candidates
            is_host = False

    # 
    print("All done!")

# Run em all
#  frb_build Hosts --frb 20181112A
# Gordon+23 hosts: frb_build Hosts --frb 20121102A,20180301A,20180916B,20180924B,20181112A,20190102C,20190520B,20190608B,20190611B,20190711A,20190714A,20191001A,20200430A,20200906A,20201124A,20210117A,20210320C,20210410D,20210807D,20211127I,20211203C,20211212A,20220105A
