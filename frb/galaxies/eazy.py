""" Module to faciliate scripting of EAZY analysis"""

import os
import warnings
from pkg_resources import resource_filename
from distutils import spawn
import subprocess

import numpy as np

from astropy.table import Table

from frb.surveys import catalog_utils

from IPython import embed

_template_list = ['br07_default','br07_goods','cww+kin','eazy_v1.0','eazy_v1.1_lines','eazy_v1.2_dusty','eazy_v1.3','pegase','pegase13']
_default_prior = 'prior_R_zmax7'
_acceptable_combos = [1,2,99,-1,'a']

# This syncs to our custom FILTERS.RES.latest file
frb_to_eazy_filters = dict(GMOS_S_r=349,
                           LRISb_V=346,
                           LRISr_I=345,
                           NOT_z=348,
                           NIRI_J=257,
                           DES_g=294,
                           DES_r=295,
                           DES_i=296,
                           DES_z=297,
                           DES_Y=298,
                           )

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
    catfile = os.path.join(input_dir, '{}.cat'.format(name))
    param_file = os.path.join(input_dir, 'zphot.param.{}'.format(name))
    translate_file = os.path.join(input_dir, 'zphot.translate.{}'.format(name))
    #
    return catfile, param_file, translate_file


def eazy_setup(input_dir, template_dir):
    """
    Setup for EAZY

    Args:
        input_dir (str):
            Path to perosonal eazy inputs/ folder  (can be relative)
        template_dir:
            Path to templates/ folder in EAZY software package

    Returns:

    """
    if not os.path.isdir(input_dir):
        os.mkdir(input_dir)
    # Copy over templates
    os.system('cp -rp {:s} {:s}'.format(template_dir, os.path.join(input_dir, '..')))
    # And link
    os.system('ln -s {:s} {:s}'.format('../templates', os.path.join(input_dir, 'templates')))
    # And FILTER files
    filter_info = os.path.join(resource_filename('frb', 'data'), 'analysis', 'EAZY', 'FILTER.RES.latest.info')
    os.system('cp -rp {:s} {:s}'.format(filter_info, input_dir))
    filter_latest = os.path.join(resource_filename('frb', 'data'), 'analysis', 'EAZY', 'FILTER.RES.latest')
    os.system('cp -rp {:s} {:s}'.format(filter_latest, input_dir))


def eazy_input_files(photom, input_dir, name, out_dir, prior_filter=None):
    """
    Write to disk a series of files needed to run EAZY
      - catalog file
      - translation file
      - param file

    Args:
        photom (dict):
            Held by an FRBGalaxy object
        input_dir (str):
            Path to eazy inputs/ folder  (can be relative)
            This is where eazy will be run
        name (str):
            Name of the source being analzyed
        out_dir (str):
            Path to eazy OUTPUT folder *relative* to the input_dir
        prior_filter (str, optional):
            If provided, use the flux in this filter for EAZY's prior
    """

    # Output filenames
    catfile, param_file, translate_file = eazy_filenames(input_dir, name)

    # Check output dir
    full_out_dir = os.path.join(input_dir, out_dir)
    if not os.path.isdir(full_out_dir):
        warnings.warn("Output directory {} does not exist, creating it!".format(full_out_dir))
        os.mkdir(full_out_dir)

    # Prior
    if prior_filter is not None:
        if prior_filter[-1] not in ['r', 'R']:
            raise IOError("Not prepared for this type of prior filter")

    # Generate the translate file
    filters = []
    codes = []
    for filt in photom.keys():
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
    outfile = os.path.join(input_dir, translate_file)
    with open(outfile, 'w') as f:
        for code, filt in zip(codes, filters):
            f.write('{} {} \n'.format(filt, code))
    print("Wrote: {}".format(outfile))

    # Catalog file
    # Generate a simple table
    phot_tbl = Table()
    phot_tbl[filters[0]] = [photom[filters[0]]]
    for filt in filters[1:]:
        phot_tbl[filt] = photom[filt]
    # Convert --
    fluxtable = catalog_utils.convert_mags_to_flux(phot_tbl, fluxunits='uJy')
    # Write
    newfs, newv = [], []
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
    with open(default_file, 'r') as df:
        df_lines = df.readlines()

    with open(param_file, 'w') as f:
        for dfline in df_lines:
            if 'CATALOG_FILE' in dfline:
                line = dfline.replace('REPLACE.cat', base_cat)
            elif prior_filter is not None and 'APPLY_PRIOR' in dfline:
                line = dfline.replace('n', 'y', 1)
            elif prior_filter is not None and 'PRIOR_FILTER' in dfline:
                line = dfline.replace('999', str(frb_to_eazy_filters[prior_filter]), 1)
            elif prior_filter is not None and 'PRIOR_FILE' in dfline:
                line = dfline  # Deal with this if we do anything other than r
            elif prior_filter is not None and 'PRIOR_ABZP' in dfline:
                line = dfline  # Deal with this if we do anything other than r
            elif 'Directory to put output files in' in dfline:  # Relative to the Input directory
                line = dfline[0:10]+dfline[10:].replace('OUTPUT', out_dir, -1)
            else:
                line = dfline
            # Write
            f.write(line)
    print("Wrote param file: {}".format(param_file))


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
    command_line = [path_to_eazy, '-p', os.path.basename(param_file), '-t', translate_file]
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


def run(catalog, infile='zphot.param', templates='eazy_v1.3', combo="a", translate_file=None,
        prior=_default_prior, prior_band=_default_band, magnitudes=False, outfolder=None,
        zmin=0.050, zmax=7.000, zstep=0.0010, prior_ABZP=23.9, logfile=None, n_min_col=5):
    """
    Runs EAZY for a given set of input parameters and saves a logfile with EAZY's console dump.
    Parameters
    ----------
    catalog: str
        Path to a photometric catalog file. ASCII table.
    infile: str, optional
        Path to the EAZY input parameter file to be created.
        Defaults to 'zphot.param' in the current working directory
    templates: str, optional
        Template set name to be used. Should be one of
        'br07_deafult','br07_goods','cww+kin','eazy_v1.0',
        'eazy_v1.1_lines','eazy_v1.2_dusty','eazy_v1.3','pegase',
        'pegase13'.
    combo: int or str, optional
        Combinations of templates to be used for analysis. Can be
        one of 1,2,99,-1 and 'a'. Read EAZY's zphot.param.default
        file for details
    translate_file: str, optional
        Translate file that instructs EAZY to recognise catalog column
        names. Defaults to the translate file supplied in this
        repository for DES+WISE photometry.
    prior: str, optional
        Name of the prior file found in the EAZY templates folder.
        Default value is 'prior_R_zmax7'.
    prior_band: int, optional
        The band on whose photometry the prior will be applied.
        EAZY by default has photometry-redshift priors for R and
        K band. Here, the default value is 295, which corresponds to
        DECam R band. Read FILTER.RES.info in the EAZY folder for
        values corresponding to other bands (not advisable to use as
        with wrong prior.)
    magnitudes: bool, optional
        True if catalog contains magnitudes as opposed to F_nu values.
        True by default.
    outfolder: str, optional
        Path to the folder where EAZY outputs will be dumped to.
        Default behaviour is to create a folder with current date and time
        as name in the current working directory.
    zmin: float, optional
        Minimum search redshift for EAZY. Default value is 0.05.
    zmax: float, optional
        Maximum search redshift for EAZY. Be careful about the prior
        file not having information beyond a redshift less than zmax.
    zstep: float, optional
        Step size of the redshift grid. (z_{i+1} = z_i+zstep*(1+z_i)).
        Default value is 0.001.
    prior_ABZP: float, optional
        Zero point redshift for the band on which prior will be applied.
        Default value is for DECam r (https://cdcvs.fnal.gov/redmine/projects/des-sci-verification/wiki/Photometry)
    logfile: str, optional
        Path to logfile with EAZY's console dump and error messages
        if any. Defaults to 'logfile.log' within the output folder.
    """
    # Validity checks.
    assert _os.path.isfile(catalog), "Invalid catalog path. The input catalog file can't be found."

    # TODO Should allow users to create their own template file in the EAZY folder to read it.
    assert templates in _template_list, "Invalid template set. Should belong to {:s}".format(str(_template_list))
    assert combo in _acceptable_combos, "Invalid combo value. Should belong to {:s}".format(str(_acceptable_combos))

    # Create an output folder
    if outfolder is None:
        outfolder = str(_dt.now().date()) + '_' + str(_dt.now().time()) + '_out'
    # Create output directory
    if not _os.path.isdir(outfolder):
        _os.mkdir(outfolder)

    # Create a log file
    if logfile is None:
        logfile = outfolder + '/logfile.log'

    print('Creating an input file for EAZY to read out from')
    if not _os.path.isfile('zphot.param'):
        result = _sub.run(_eazyexec, stdout=_sub.PIPE, stderr=_sub.PIPE)  # Run EAZY and capture output
        with open(logfile, "w+") as fstream:  # Write to logfile
            fstream.write(result.stdout.decode('utf-8'))
        # Read table
        in_tab = _pd.read_csv('zphot.param.default', delim_whitespace=True, comment="#",
                              header=None, names=('eazy_par', 'par_val'), skiprows=[11, 12, 13, 14])
        # Delete temporary file
        _os.system('rm -f zphot.param.default')
    else:
        with open(logfile, "w+") as fstream:  # Write to logfile
            fstream.write("Creating Input file\n")
        in_tab = _pd.read_table('zphot.param', delim_whitespace=True, comment="#",
                                header=None, names=('eazy_par', 'par_val'))

    # Change default parameters to reflect current values
    in_tab.par_val[in_tab.eazy_par == 'FILTERS_RES'] = _eazy_path + '/inputs/FILTER.RES.latest'
    in_tab.par_val[in_tab.eazy_par == 'TEMPLATES_FILE'] = _eazy_path + '/templates/' + templates + ".spectra.param"
    in_tab.par_val[in_tab.eazy_par == 'TEMP_ERR_FILE'] = _eazy_path + '/templates/TEMPLATE_ERROR.eazy_v1.0'
    in_tab.par_val[in_tab.eazy_par == 'TEMPLATE_COMBOS'] = combo
    in_tab.par_val[in_tab.eazy_par == 'WAVELENGTH_FILE'] = _eazy_path + '/templates/EAZY_v1.1_lines/lambda_v1.1.def'
    in_tab.par_val[in_tab.eazy_par == 'LAF_FILE'] = _eazy_path + '/templates/LAFcoeff.txt'
    in_tab.par_val[in_tab.eazy_par == 'DLA_FILE'] = _eazy_path + '/templates/DLAcoeff.txt'
    in_tab.par_val[in_tab.eazy_par == 'CATALOG_FILE'] = catalog
    if magnitudes:
        in_tab.par_val[in_tab.eazy_par == 'MAGNITUDES'] = 'y'
    else:
        in_tab.par_val[in_tab.eazy_par == 'MAGNITUDES'] = 'n'
    in_tab.par_val[in_tab.eazy_par == 'N_MIN_COLORS'] = n_min_col
    in_tab.par_val[in_tab.eazy_par == 'OUTPUT_DIRECTORY'] = outfolder
    in_tab.par_val[in_tab.eazy_par == 'APPLY_PRIOR'] = 1
    in_tab.par_val[in_tab.eazy_par == 'PRIOR_FILE'] = _eazy_path + '/templates/' + prior + '.dat'
    in_tab.par_val[in_tab.eazy_par == 'PRIOR_FILTER'] = prior_band
    in_tab.par_val[in_tab.eazy_par == 'PRIOR_ABZP'] = prior_ABZP
    in_tab.par_val[in_tab.eazy_par == 'Z_MIN'] = zmin
    in_tab.par_val[in_tab.eazy_par == 'Z_MAX'] = zmax
    in_tab.par_val[in_tab.eazy_par == 'Z_STEP'] = zstep
    in_tab.par_val[in_tab.eazy_par == 'H0'] = _p15.H0.value
    in_tab.par_val[in_tab.eazy_par == 'OMEGA_M'] = _p15.Om0
    in_tab.par_val[in_tab.eazy_par == 'OMEGA_L'] = _p15.Ode0

    # Create infile
    in_tab.to_csv(infile, header=False, index=False, sep="\t")

    # Just use the default translate file if not given
    if translate_file is None:
        translate_file = _os.path.relpath('../zphot.translate')
        _warnings.warn_explicit(
            "Using default translate file from {:s}. May not work for your columns.".format(translate_file),
            category=UserWarning, filename='run_eazy.py', lineno=150)

    # Run EAZY
    print("Running EAZY ...")
    # TODO: add some progressbar here
    eazy_out = _sub.run([_eazyexec, '-p', infile, '-t', translate_file], stdout=_sub.PIPE, stderr=_sub.PIPE)
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
    return eazy_out.returncode


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

