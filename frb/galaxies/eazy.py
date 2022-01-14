""" Module to faciliate scripting of EAZY analysis"""

import os
import warnings
from pkg_resources import resource_filename
from distutils import spawn
import subprocess

import numpy as np
import pandas

from astropy.table import Table

from frb.surveys import catalog_utils
from frb import defs

from IPython import embed

# Necessary because people might not execute eazy from src
# but might have a copy in bin or some other location.
try:
    _eazy_root = os.environ['EAZYDIR']
except KeyError:
    warnings.warn('Please define the variable EAZYDIR in your environment pointing to the EAZY folder.')

_template_list = ['br07_default','br07_goods','cww+kin','eazy_v1.0','eazy_v1.1_lines','eazy_v1.2_dusty','eazy_v1.3','pegase','pegase13']
_acceptable_priors = ['prior_R_zmax7', 'prior_K_zmax7', 'prior_R_extend', 'prior_K_extend'] # F160W_TAO not included yet.
_default_prior = 'prior_R_zmax7'
_acceptable_combos = [1,2,99,-1,'a']

# This syncs to our custom FILTERS.RES.latest file
frb_to_eazy_filters = {"GMOS_S_r":349,
                           "LRISb_V":346,
                           "LRISr_I":345,
                           "NOT_z":348,
                           "NIRI_J":257,
                           "DECaL_g":294, # Turns out these are legacy transmission curves
                           "DECaL_r":295,
                           "DECaL_z":296,
                           "DES_u":351, # Added DR1 filter curves
                           "DES_g":352,
                           "DES_r":353,
                           "DES_i":354,
                           "DES_z":355,
                           "DES_y":356,
                           "SDSS_u":156,
                           "SDSS_g":157,
                           "SDSS_r":158,
                           "SDSS_i":159,
                           "SDSS_z":160,
                           "WISE_W1":244,
                           "WISE_W2":245,
                           "WISE_W3":246,
                           "WISE_W4":247,
                           "Pan-STARRS_g":334,
                           "Pan-STARRS_r":335,
                           "Pan-STARRS_i":336,
                           "Pan-STARRS_z":337,
                           "Pan-STARRS_y":338,
                           "VISTA_Y":256,
                           "VISTA_J":257,
                           "VISTA_H":258,
                           "VISTA_Ks":259
                            }

def eazy_filenames(input_dir, name):
    """
    Generate names for EAZY files

    Args:
        input_dir (str):
            Path to eazy inputs/ folder  (can be relative)
            This is where eazy will be run
        name (str):
            Name of the source being analzyed

    Returns:
        tuple:  catalog_filename, parameter_filename, translate_file

    """
    if not os.path.isdir(input_dir):
        os.mkdir(input_dir)
    catfile = os.path.join(input_dir, '{}.cat'.format(name))
    param_file = os.path.join(input_dir, 'zphot.param.{}'.format(name))
    translate_file = os.path.join(input_dir, 'zphot.translate.{}'.format(name))
    #
    return catfile, param_file, translate_file


def eazy_setup(input_dir, template_dir=None):
    """
    Setup for EAZY

    Args:
        input_dir (str):
            Path to personal eazy inputs/ folder  (can be relative)
        template_dir(str, optional):
            Path to templates/ folder in EAZY software package.
            If not given, it looks for the folder of `eazy`,
            the executable and navigates from there.

    Returns:

    """
    if template_dir is None:
        template_dir = os.path.join(_eazy_root, "templates")
    if not os.path.isdir(input_dir):
        os.mkdir(input_dir)
    # And link
    os.system('ln -s {:s} {:s}'.format(template_dir, os.path.join(input_dir, 'templates')))
    # And FILTER files
    filter_info = os.path.join(resource_filename('frb', 'data'), 'analysis', 'EAZY', 'FILTER.RES.latest.info')
    os.system('cp -rp {:s} {:s}'.format(filter_info, input_dir))
    filter_latest = os.path.join(resource_filename('frb', 'data'), 'analysis', 'EAZY', 'FILTER.RES.latest')
    os.system('cp -rp {:s} {:s}'.format(filter_latest, input_dir))
    return

def eazy_input_files(photom, input_dir, name, out_dir, prior_filter=None,
                     templates='eazy_v1.3', combo="a", cosmo=defs.frb_cosmo,
                     magnitudes=False, prior=_default_prior,
                     zmin=0.050, zmax=7.000, zstep=0.0010, prior_ABZP=23.9,
                     n_min_col=5, write_full_table=False):
    """
    Write to disk a series of files needed to run EAZY
      - catalog file
      - translation file
      - param file

    Args:
        photom (dict or Table):
            Held by an FRBGalaxy object
        input_dir (str):
            Path to eazy inputs/ folder  (can be relative)
            This is where eazy will be run
        name (str):
            Name of the source being analzyed
        out_dir (str):
            Path to eazy OUTPUT folder *relative* to the input_dir
        id_col (str, optional):
            Column name to be used as the ID. Looks for a
            column with "id" in its name by default.
        prior_filter (str, optional):
            If provided, use the flux in this filter for EAZY's prior
        templates (str, optional):
            Template set name to be used. Should be one of
            'br07_deafult','br07_goods','cww+kin','eazy_v1.0',
            'eazy_v1.1_lines','eazy_v1.2_dusty','eazy_v1.3','pegase',
            'pegase13'.
        combo (int or str, optional):
            Combinations of templates to be used for analysis. Can be
            one of 1,2,99,-1 and 'a'. Read EAZY's zphot.param.default
            file for details
        cosmo (astropy.cosmology, optional):
            Defaults to Repo cosmology
        prior (str, optional):
            Name of the prior file found in the EAZY templates folder.
            Default value is 'prior_R_zmax7'.
        magnitudes (bool, optional):
            True if catalog contains magnitudes as opposed to F_nu values.
        zmin (float, optional):
            Minimum search redshift for EAZY. Default value is 0.05.
        zmax (float, optional):
            Maximum search redshift for EAZY. Be careful about the prior
            file not having information beyond a redshift less than zmax.
        zstep (float, optional):
            Step size of the redshift grid. (z_{i+1} = z_i+zstep*(1+z_i)).
            Default value is 0.001.
        prior_ABZP (float, optional):
            Zero point redshift for the band on which prior will be applied.
            Default value is for DECam r (https://cdcvs.fnal.gov/redmine/projects/des-sci-verification/wiki/Photometry)
        write_full_table (bool, optional):
            Are you trying to use this function for a table of objects instead of
            a single object? If so, set this to True.
    """
    # Output filenames
    catfile, param_file, translate_file = eazy_filenames(input_dir, name)

    # Check output dir
    full_out_dir = os.path.join(input_dir, out_dir)
    if not os.path.isdir(full_out_dir):
        warnings.warn("Output directory {} does not exist, creating it!".format(full_out_dir))
        os.makedirs(full_out_dir)

    # Prior
    if prior_filter is not None:
        assert prior in _acceptable_priors, "Allowed priors are {}".format(_acceptable_priors)
        if prior_filter.split('_')[-1] not in ['r', 'R', 'k', 'K', 'Ks']:
            raise IOError("Not prepared for this type of prior filter")
    
    # Test combo
    assert combo in _acceptable_combos, "Allowed values of 'combo' are {}".format(_acceptable_combos)

    # Create ID column if it's not there.
    if id_col not in photom.colnames:
        photom[id_col] = np.arange(len(photom))
    else:
        assert id_col in photom.colnames, "Could not find {:s} in the photometry table.".format(id_col)

    # Generate the translate file
    filters = []
    codes = []
    for filt in photom.colnames:
        if 'EBV' in filt:
            continue
        if 'err' in filt:
            ifilt = filt[:-4]
            pref = 'E'
        else:
            ifilt = filt
            pref = 'F'
        # Check
        if ifilt not in frb_to_eazy_filters.keys():
            warnings.warn("Filter {} not in our set.  Add it to frb.galaxies.eazy.frb_to_eazy_filters".format(ifilt))
            continue
        # Grab it
        code = '{}{}'.format(pref,frb_to_eazy_filters[ifilt])
        # Append
        filters.append(filt)
        codes.append(code)
    # Do it
    with open(translate_file, 'w') as f:
        f.write(id_col+" id\n")
        for code, filt in zip(codes, filters):
            f.write('{} {} \n'.format(filt, code))
    print("Wrote: {}".format(translate_file))

    # Catalog file
    # Generate a simple table
    phot_tbl = Table()
    if np.isscalar(photom[filters[0]]):
        phot_tbl[filters[0]] = [photom[filters[0]]]
    else:
        phot_tbl[filters[0]] = photom[filters[0]]
    #import pdb;pdb.set_trace()
    for filt in filters[1:]:
        phot_tbl[filt] = photom[filt]
    # Convert --
    fluxtable = catalog_utils.convert_mags_to_flux(photom, fluxunits='uJy')
    # Write
    newfs, newv = [], []
    if write_full_table:
        fluxtable.write(catfile, format="ascii.commented_header", overwrite=True)
    else:
        for key in fluxtable.keys():
            newfs.append(key)
            newv.append(str(fluxtable[key].data[0]))
        with open(catfile, 'w') as f:
            # Filters
            allf = ' '.join(newfs)
            f.write('# {} \n'.format(allf))
            # Values
            f.write(' '.join(newv))
    print("Wrote catalog file: {}".format(catfile))
    base_cat = os.path.basename(catfile)

    # Input file
    default_file = os.path.join(resource_filename('frb', 'data'), 'analysis', 'EAZY', 'zphot.param.default')
    #with open(default_file, 'r') as df:
    #    df_lines = df.readlines()
    in_tab = pandas.read_table(default_file, delim_whitespace=True, comment="#",
                            header=None, names=('eazy_par', 'par_val'))

    # Expect it lands in src/ # Won't work if someone has put eazy in their bin folder.
    #_eazy_path = os.path.abspath(os.path.realpath(spawn.find_executable('eazy')))
    #_eazy_root = _eazy_path[0:_eazy_path.find('src')]

    # Change default parameters to reflect current values
    in_tab.par_val[in_tab.eazy_par == 'FILTERS_RES'] = os.path.join(resource_filename('frb', 'data'), 'analysis', 'EAZY', 'FILTER.RES.latest')
    in_tab.par_val[in_tab.eazy_par == 'TEMPLATES_FILE'] = os.path.join(_eazy_root, 'templates/' + templates + ".spectra.param")
    in_tab.par_val[in_tab.eazy_par == 'TEMP_ERR_FILE'] = os.path.join(_eazy_root,'templates/TEMPLATE_ERROR.eazy_v1.0')
    in_tab.par_val[in_tab.eazy_par == 'TEMPLATE_COMBOS'] = combo
    in_tab.par_val[in_tab.eazy_par == 'WAVELENGTH_FILE'] = os.path.join(_eazy_root,'templates/EAZY_v1.1_lines/lambda_v1.1.def')
    in_tab.par_val[in_tab.eazy_par == 'LAF_FILE'] = os.path.join(_eazy_root,'templates/LAFcoeff.txt')
    in_tab.par_val[in_tab.eazy_par == 'DLA_FILE'] = os.path.join(_eazy_root,'templates/DLAcoeff.txt')
    in_tab.par_val[in_tab.eazy_par == 'CATALOG_FILE'] = base_cat
    if magnitudes:
        in_tab.par_val[in_tab.eazy_par == 'MAGNITUDES'] = 'y'
    else:
        in_tab.par_val[in_tab.eazy_par == 'MAGNITUDES'] = 'n'
    in_tab.par_val[in_tab.eazy_par == 'N_MIN_COLORS'] = n_min_col
    in_tab.par_val[in_tab.eazy_par == 'OUTPUT_DIRECTORY'] = out_dir
    # Prior
    if prior_filter is not None:
        in_tab.par_val[in_tab.eazy_par == 'APPLY_PRIOR'] = 'y'
        in_tab.par_val[in_tab.eazy_par == 'PRIOR_FILE'] = os.path.join(_eazy_root,'templates/' + prior + '.dat')
        in_tab.par_val[in_tab.eazy_par == 'PRIOR_FILTER'] = str(frb_to_eazy_filters[prior_filter])
        in_tab.par_val[in_tab.eazy_par == 'PRIOR_ABZP'] = prior_ABZP
    in_tab.par_val[in_tab.eazy_par == 'Z_MIN'] = zmin
    in_tab.par_val[in_tab.eazy_par == 'Z_MAX'] = zmax
    in_tab.par_val[in_tab.eazy_par == 'Z_STEP'] = zstep
    in_tab.par_val[in_tab.eazy_par == 'H0'] = cosmo.H0.value
    in_tab.par_val[in_tab.eazy_par == 'OMEGA_M'] = cosmo.Om0
    in_tab.par_val[in_tab.eazy_par == 'OMEGA_L'] = cosmo.Ode0

    # Create infile
    in_tab.to_csv(param_file, header=False, index=False, sep="\t")

#    with open(param_file, 'w') as f:
#        for dfline in df_lines:
#            if 'CATALOG_FILE' in dfline:
#                line = dfline.replace('REPLACE.cat', base_cat)
#            elif prior_filter is not None and 'APPLY_PRIOR' in dfline:
#                line = dfline.replace('n', 'y', 1)
#            elif prior_filter is not None and 'PRIOR_FILTER' in dfline:
#                line = dfline.replace('999', str(frb_to_eazy_filters[prior_filter]), 1)
#            elif prior_filter is not None and 'PRIOR_FILE' in dfline:
#                line = dfline  # Deal with this if we do anything other than r
#            elif prior_filter is not None and 'PRIOR_ABZP' in dfline:
#                line = dfline  # Deal with this if we do anything other than r
#            elif 'Directory to put output files in' in dfline:  # Relative to the Input directory
#                line = dfline[0:10]+dfline[10:].replace('OUTPUT', out_dir, -1)
#            else:
#                line = dfline
            # Write
#            f.write(line)
#    print("Wrote param file: {}".format(param_file))

    # Generate a soft link to the templates
    template_folder = os.path.join(_eazy_root, 'templates')
    link = os.path.join(input_dir, '..', 'templates')
    if not os.path.isdir(link):
        os.symlink(template_folder, link)
    return


def run_eazy(input_dir, name, logfile):
    """
    Find and run the EAZY executable on the files

    Args:
        input_dir (str):
            Path to eazy inputs/ folder  (can be relative)
            This is where eazy will be run
        name (str):
            Name of the source being analzyed
        logfile (str):
    """
    _, param_file, translate_file = eazy_filenames(input_dir, name)

    # Find the eazy executable
    path_to_eazy = spawn.find_executable('eazy')
    if path_to_eazy is None:
        raise ValueError("You must have eazy in your Unix path..")
    # Run it!
    command_line = [path_to_eazy, '-p', os.path.basename(param_file),
                    '-t', os.path.basename(translate_file)]
    #
    eazy_out = subprocess.run(command_line, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, cwd=input_dir)
    # Dump stdout to logfile
    with open(logfile, "a") as fstream:
        fstream.write(eazy_out.stdout.decode('utf-8'))

    # Check if the process ran successfully
    if eazy_out.returncode == 0:
        print("EAZY ran successfully!")
    elif eazy_out.returncode == -6:
        print(
            "Try to put your input parameter and translate files in a place with a shorter relative paths. Otherwise, you get this buffer overflow.")
    else:
        # Dump stderr to logfile
        with open(logfile, "a") as fstream:
            fstream.write("ERROR\n------------------------------------\n")
            fstream.write(eazy_out.stderr.decode('utf-8'))
        print("Somethign went wrong. Look at {:s} for details".format(logfile))


def eazy_stats(zgrid, pzi):
    """
    Calculate the 'best' zphot and error

    Args:
        zgrid (np.ndarray):
        pzi (np.ndarray):

    Returns:

    """
    # p(z) weighted
    zphot = np.sum(zgrid * pzi) / np.sum(pzi)
    # Uncertainty
    cum_pzi = np.cumsum(pzi) / np.sum(pzi)
    l68 = zgrid[np.argmin(np.abs(cum_pzi - 0.16))]
    u68 = zgrid[np.argmin(np.abs(cum_pzi - (1 - 0.16)))]
    #
    return zphot, (u68-l68)/2.


####################################################
#  The following are taken from threedhst by Brummer
####################################################

def readEazyBinary(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):
    """
    tempfilt, coeffs, temp_sed, pz = readEazyBinary(MAIN_OUTPUT_FILE='photz', \
                                                OUTPUT_DIRECTORY='./OUTPUT', \
                                                CACHE_FILE = 'Same')

    Read Eazy BINARY_OUTPUTS files into structure data.

    If the BINARY_OUTPUTS files are not in './OUTPUT', provide either a relative or absolute path
    in the OUTPUT_DIRECTORY keyword.

    By default assumes that CACHE_FILE is MAIN_OUTPUT_FILE+'.tempfilt'.
    Specify the full filename if otherwise.
    """

    # root='COSMOS/OUTPUT/cat3.4_default_lines_zp33sspNoU'

    root = OUTPUT_DIRECTORY + '/' + MAIN_OUTPUT_FILE

    ###### .tempfilt
    if CACHE_FILE == 'Same':
        CACHE_FILE = root + '.tempfilt'

    if os.path.exists(CACHE_FILE) is False:
        print(('File, %s, not found.' % (CACHE_FILE)))
        return -1, -1, -1, -1

    f = open(CACHE_FILE, 'rb')

    s = np.fromfile(file=f, dtype=np.int32, count=4)
    NFILT = s[0]
    NTEMP = s[1]
    NZ = s[2]
    NOBJ = s[3]
    tempfilt = np.fromfile(file=f, dtype=np.double, count=NFILT * NTEMP * NZ).reshape((NZ, NTEMP, NFILT)).transpose()
    lc = np.fromfile(file=f, dtype=np.double, count=NFILT)
    zgrid = np.fromfile(file=f, dtype=np.double, count=NZ)
    fnu = np.fromfile(file=f, dtype=np.double, count=NFILT * NOBJ).reshape((NOBJ, NFILT)).transpose()
    efnu = np.fromfile(file=f, dtype=np.double, count=NFILT * NOBJ).reshape((NOBJ, NFILT)).transpose()

    f.close()

    tempfilt = {'NFILT': NFILT, 'NTEMP': NTEMP, 'NZ': NZ, 'NOBJ': NOBJ, \
                'tempfilt': tempfilt, 'lc': lc, 'zgrid': zgrid, 'fnu': fnu, 'efnu': efnu}

    ###### .coeff
    f = open(root + '.coeff', 'rb')

    s = np.fromfile(file=f, dtype=np.int32, count=4)
    NFILT = s[0]
    NTEMP = s[1]
    NZ = s[2]
    NOBJ = s[3]
    coeffs = np.fromfile(file=f, dtype=np.double, count=NTEMP * NOBJ).reshape((NOBJ, NTEMP)).transpose()
    izbest = np.fromfile(file=f, dtype=np.int32, count=NOBJ)
    tnorm = np.fromfile(file=f, dtype=np.double, count=NTEMP)

    f.close()

    coeffs = {'NFILT': NFILT, 'NTEMP': NTEMP, 'NZ': NZ, 'NOBJ': NOBJ, \
              'coeffs': coeffs, 'izbest': izbest, 'tnorm': tnorm}

    ###### .temp_sed
    f = open(root + '.temp_sed', 'rb')
    s = np.fromfile(file=f, dtype=np.int32, count=3)
    NTEMP = s[0]
    NTEMPL = s[1]
    NZ = s[2]
    templam = np.fromfile(file=f, dtype=np.double, count=NTEMPL)
    temp_seds = np.fromfile(file=f, dtype=np.double, count=NTEMPL * NTEMP).reshape((NTEMP, NTEMPL)).transpose()
    da = np.fromfile(file=f, dtype=np.double, count=NZ)
    db = np.fromfile(file=f, dtype=np.double, count=NZ)

    f.close()

    temp_sed = {'NTEMP': NTEMP, 'NTEMPL': NTEMPL, 'NZ': NZ, \
                'templam': templam, 'temp_seds': temp_seds, 'da': da, 'db': db}

    ###### .pz
    if os.path.exists(root + '.pz'):
        f = open(root + '.pz', 'rb')
        s = np.fromfile(file=f, dtype=np.int32, count=2)
        NZ = s[0]
        NOBJ = s[1]
        chi2fit = np.fromfile(file=f, dtype=np.double, count=NZ * NOBJ).reshape((NOBJ, NZ)).transpose()

        ### This will break if APPLY_PRIOR No
        s = np.fromfile(file=f, dtype=np.int32, count=1)

        if len(s) > 0:
            NK = s[0]
            kbins = np.fromfile(file=f, dtype=np.double, count=NK)
            priorzk = np.fromfile(file=f, dtype=np.double, count=NZ * NK).reshape((NK, NZ)).transpose()
            kidx = np.fromfile(file=f, dtype=np.int32, count=NOBJ)
            pz = {'NZ': NZ, 'NOBJ': NOBJ, 'NK': NK, 'chi2fit': chi2fit, 'kbins': kbins, 'priorzk': priorzk,
                  'kidx': kidx}
        else:
            pz = None

        f.close()

    else:
        pz = None

    if False:
        f = open(root + '.zbin', 'rb')
        s = np.fromfile(file=f, dtype=np.int32, count=1)
        NOBJ = s[0]
        z_a = np.fromfile(file=f, dtype=np.double, count=NOBJ)
        z_p = np.fromfile(file=f, dtype=np.double, count=NOBJ)
        z_m1 = np.fromfile(file=f, dtype=np.double, count=NOBJ)
        z_m2 = np.fromfile(file=f, dtype=np.double, count=NOBJ)
        z_peak = np.fromfile(file=f, dtype=np.double, count=NOBJ)
        f.close()

    ###### Done.
    return tempfilt, coeffs, temp_sed, pz


def getEazyPz(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same', binaries=None,
              get_prior=False, get_chi2=False):
    """
    zgrid, pz = getEazyPz(idx, \
                      MAIN_OUTPUT_FILE='photz', \
                      OUTPUT_DIRECTORY='./OUTPUT', \
                      CACHE_FILE='Same', binaries=None)

    Get Eazy p(z) for object #idx.

    To avoid re-reading the binary files, supply binaries = (tempfilt, pz)

    """
    if binaries is None:
        tempfilt, coeffs, temp_seds, pz = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                                         OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                                         CACHE_FILE=CACHE_FILE)
    else:
        tempfilt, pz = binaries

    if pz is None:
        return None, None

    ###### Get p(z|m) from prior grid
    kidx = pz['kidx'][idx]
    # print kidx, pz['priorzk'].shape
    if (kidx > 0) & (kidx < pz['priorzk'].shape[1]):
        prior = pz['priorzk'][:, kidx]
    else:
        prior = np.ones(pz['NZ'])

    if get_chi2:
        if get_prior:
            if get_prior:
                return tempfilt['zgrid'], pz['chi2fit'][:, idx], prior
            else:
                return tempfilt['zgrid'], pz['chi2fit'][:, idx]

    ###### Convert Chi2 to p(z)
    pzi = np.exp(-0.5 * (pz['chi2fit'][:, idx] - min(pz['chi2fit'][:, idx]))) * prior  # *(1+tempfilt['zgrid'])

    if np.sum(pzi) > 0:
        pzi /= np.trapz(pzi, tempfilt['zgrid'])

    ###### Done
    if get_prior:
        return tempfilt['zgrid'], pzi, prior
    else:
        return tempfilt['zgrid'], pzi