""" Module to faciliate scripting of EAZY analysis"""

import os
import warnings
from pkg_resources import resource_filename
from distutils import spawn
import subprocess

import numpy as np
import pandas

from astropy.table import Table
from astropy.cosmology import Planck15

from frb.surveys import catalog_utils

from IPython import embed

# Necessary because people might not execute eazy from src
# but might have a copy in bin or some other location.
try:
    _eazy_root = os.environ['EAZYDIR']
except KeyError:
    import pdb; pdb.set_trace()
    raise AssertionError('Please define the variable EAZYDIR in your environment pointing to the EAZY folder.')

_template_list = ['br07_default','br07_goods','cww+kin','eazy_v1.0','eazy_v1.1_lines','eazy_v1.2_dusty','eazy_v1.3','pegase','pegase13']
_acceptable_priors = ['prior_R_zmax7', 'prior_K_zmax7', 'prior_R_extend', 'prior_K_extend'] # F160W_TAO not included yet.
_default_prior = 'prior_R_zmax7'
_acceptable_combos = [1,2,99,-1,'a']

# This syncs to our custom FILTERS.RES.latest file
frb_to_eazy_filters = dict(GMOS_S_r=349,
                           LRISb_V=346,
                           LRISr_I=345,
                           NOT_z=348,
                           NIRI_J=257,
                           DECaL_g=294, # Turns out these are legacy transmission curves
                           DECaL_r=295,
                           DECaL_z=296,
                           DES_u=351, # Added DR1 filter curves
                           DES_g=352,
                           DES_r=353,
                           DES_i=354,
                           DES_z=355,
                           DES_y=356,
                           SDSS_u=156,
                           SDSS_g=157,
                           SDSS_r=158,
                           SDSS_i=159,
                           SDSS_z=160,
                           WISE_W1=244,
                           WISE_W2=245,
                           WISE_W3=246,
                           WISE_W4=247
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
                     templates='eazy_v1.3', combo="a", cosmo=None,
                     magnitudes=False, prior=_default_prior,
                     zmin=0.050, zmax=7.000, zstep=0.0010, prior_ABZP=23.9,
                     n_min_col=5):
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
            Defaults to Planck15
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
    """
    #
    if cosmo is None:
        cosmo = Planck15

    # Output filenames
    catfile, param_file, translate_file = eazy_filenames(input_dir, name)

    # Check output dir
    full_out_dir = os.path.join(input_dir, out_dir)
    if not os.path.isdir(full_out_dir):
        warnings.warn("Output directory {} does not exist, creating it!".format(full_out_dir))
        os.mkdir(full_out_dir)

    # Prior
    if prior_filter is not None:
        assert prior in _acceptable_priors, "Allowed priors are {}".format(_acceptable_priors)
        if prior_filter[-1] not in ['r', 'R', 'k', 'K']:
            raise IOError("Not prepared for this type of prior filter")
    
    # Test combo
    assert combo in _acceptable_combos, "Allowed values of 'combo' are {}".format(_acceptable_combos)

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
    with open(translate_file, 'w') as f:
        for code, filt in zip(codes, filters):
            f.write('{} {} \n'.format(filt, code))
    print("Wrote: {}".format(translate_file))

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

def getEazySED(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same', scale_flambda=1.e-17, verbose=False, individual_templates=False):
    """
lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
     getEazySED(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same')
    
    Get best-fit Eazy template for object number 'idx' from the specified Eazy output files. 

    Output variables are as follows:
        
        lambdaz: full best-fit template (observed) wavelength, interpolated at WAVELENGTH_GRID
        temp_sed:          "        "              flux (F_lambda)
        lci: filter pivot wavelengths
        fobs: observed fluxes, including zeropoint offsets if used, F_lambda
        efobs: observed flux errors,    "            "        "        "
    """
    tempfilt, coeffs, temp_seds, pz = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE = CACHE_FILE)
    
    ##### Apply zeropoint factors
    param = EazyParam(PARAM_FILE=OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.param')
    fnumbers = np.zeros(len(param.filters), dtype=np.int)
    for i in range(len(fnumbers)):
        fnumbers[i] = int(param.filters[i].fnumber)
    
    zpfile = OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zeropoint'
    if os.path.exists(zpfile):
        zpfilts, zpf_file = np.loadtxt(zpfile, unpack=True, dtype=np.str)                                    
        zpf = np.ones(tempfilt['NFILT'])
        for i in range(len(zpfilts)):
            match = fnumbers == int(zpfilts[i][1:])
            zpf[match] = np.float(zpf_file[i])
    else:
        zpf = np.ones(tempfilt['NFILT'])

    zpfactors = np.dot(zpf.reshape(tempfilt['NFILT'],1),\
                       np.ones(tempfilt['NOBJ']).reshape(1,tempfilt['NOBJ']))

    if verbose:
        print(zpf)
        
    tempfilt['fnu'] *= zpfactors
    tempfilt['efnu'] *= zpfactors
    
    lci = tempfilt['lc'].copy()
    
    params = EazyParam(PARAM_FILE=OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.param')
    abzp = np.float(params['PRIOR_ABZP'])
        
    # fobs = tempfilt['fnu'][:,idx]/(lci/5500.)**2*flam_factor
    # efobs = tempfilt['efnu'][:,idx]/(lci/5500.)**2*flam_factor
    ### Physical f_lambda fluxes, 10**-17 ergs / s / cm2 / A
    if scale_flambda:
        flam_factor = 10**(-0.4*(params['PRIOR_ABZP']+48.6))*3.e18/scale_flambda
    else:
        flam_factor = 5500.**2
    
    missing = (tempfilt['fnu'][:,idx] < -99) | (tempfilt['efnu'][:,idx] < 0)
    fobs = tempfilt['fnu'][:,idx]/lci**2*flam_factor
    efobs = tempfilt['efnu'][:,idx]/lci**2*flam_factor
    fobs[missing] = -99
    efobs[missing] = -99
    #print lci, tempfilt['fnu'][:,idx], tempfilt['efnu'][:,idx]
    
    ##### Broad-band SED
    obs_sed = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][idx]],\
                     coeffs['coeffs'][:,idx])/(lci)**2*flam_factor
    
    zi = tempfilt['zgrid'][coeffs['izbest'][idx]]
    
    ###### Full template SED, observed frame
    lambdaz = temp_seds['templam']*(1+zi)
    temp_sed = np.dot(temp_seds['temp_seds'],coeffs['coeffs'][:,idx])
    if individual_templates:
        temp_sed = temp_seds['temp_seds']*coeffs['coeffs'][:,idx]
    
    temp_sed /= (1+zi)**2
    
    temp_sed *= (1/5500.)**2*flam_factor
    
    ###### IGM absorption
    lim1 = np.where(temp_seds['templam'] < 912)
    lim2 = np.where((temp_seds['templam'] >= 912) & (temp_seds['templam'] < 1026))
    lim3 = np.where((temp_seds['templam'] >= 1026) & (temp_seds['templam'] < 1216))
    
    if lim1[0].size > 0: temp_sed[lim1] *= 0.
    if lim2[0].size > 0: temp_sed[lim2] *= 1.-temp_seds['db'][coeffs['izbest'][idx]]
    if lim3[0].size > 0: temp_sed[lim3] *= 1.-temp_seds['da'][coeffs['izbest'][idx]]
        
    ###### Done
    return lambdaz, temp_sed, lci, obs_sed, fobs, efobs

def plotExampleSED(idx=20, writePNG=True, MAIN_OUTPUT_FILE = 'photz', OUTPUT_DIRECTORY = 'OUTPUT', CACHE_FILE = 'Same', lrange=[3000,8.e4], axes=None, individual_templates=False, fnu=False, show_pz=True, snlim=2, scale_flambda=1.e-17, setrc=True, show_rest=False):
    """
PlotSEDExample(idx=20)

    Plot an example Eazy best-fit SED.
    """

    #zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    #zout = catIO.Readfile(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    zout = catIO.Table(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')    
    #qz = np.where(zout.z_spec > 0)[0]
    print((zout.filename))
    qz = np.arange(len(zout['id']))
    
    if show_rest:
        z_peak = zout['z_peak'][idx]
        xrest = 1+z_peak
        rest_label = r'$\mathrm{rest}$'
    else:
        xrest = 1.
        rest_label = ''
        
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        getEazySED(qz[idx], MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                          OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                          CACHE_FILE = CACHE_FILE, individual_templates=individual_templates, scale_flambda=scale_flambda)
    
    zgrid, pz = getEazyPz(qz[idx], MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                   OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                   CACHE_FILE = CACHE_FILE)
    ##### plot defaults
    #rc('font',**{'family':'serif','serif':['Times']})
    if setrc:
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.serif'] = ['Helvetica']
        #plt.rcParams['ps.useafm'] = True
        plt.rcParams['patch.linewidth'] = 0.
        plt.rcParams['patch.edgecolor'] = 'black'
        #plt.rcParams['text.usetex'] = True
        plt.rcParams['text.usetex'] = False
        #plt.rcParams['text.latex.preamble'] = ''

    ##### start plot
    if axes is None:
        fig = plt.figure(figsize=[8,4],dpi=100)
        fig.subplots_adjust(wspace=0.18, hspace=0.0,left=0.09,bottom=0.15,right=0.98,top=0.98)
    
    #### Plot parameters
    plotsize=35
    alph=0.9
    
    if fnu:
        temp_sed *= (lambdaz / 5500.)**2
        fobs *= (lci/5500.)**2
        efobs *= (lci/5500.)**2
        obs_sed *= (lci/5500.)**2

    #### Full best-fit template
    if axes is None:
        ax = fig.add_subplot(121)
        axp = fig.add_subplot(122)
    else:
        ax = axes[0]
        axp = axes[1]
        
    if individual_templates:
        ax.plot(lambdaz/xrest, temp_sed, linewidth=1.0, color='blue',alpha=0.4)
        ax.plot(lambdaz/xrest, temp_sed.sum(axis=1), linewidth=1.0, color='blue',alpha=alph)
    else:
        ax.plot(lambdaz/xrest, temp_sed, linewidth=1.5, color='blue',alpha=alph*0.8, zorder=-3)
    
    #### template fluxes integrated through the filters
    ax.scatter(lci/xrest, obs_sed,
               c='red',marker='o',s=plotsize,alpha=alph, zorder=-1)

    #### Observed fluxes w/ errors
    #ax.errorbar(lci,fobs,yerr=efobs,ecolor=None,
    #           color='black',fmt='o',alpha=alph)
    #
    # ax.errorbar(lci, fobs, yerr=efobs, ecolor='black',
    #            color='black',fmt='o',alpha=alph, markeredgecolor='black', markerfacecolor='None', markeredgewidth=1.5, ms=8, zorder=1)

    highsn = fobs/efobs > snlim
    ax.errorbar(lci[highsn]/xrest, fobs[highsn], yerr=efobs[highsn], ecolor='black',
               color='black',fmt='o',alpha=alph, markeredgecolor='black', markerfacecolor='None', markeredgewidth=1.5, ms=8, zorder=2)
    #
    ax.errorbar(lci[~highsn]/xrest, fobs[~highsn], yerr=efobs[~highsn], ecolor='0.7',
               color='black',fmt='o',alpha=alph, markeredgecolor='0.7', markerfacecolor='None', markeredgewidth=1.5, ms=8, zorder=1)
    
    for i in range(len(lci)):
        print(('%f %e %e %e' %(lci[i], obs_sed[i], fobs[i], efobs[i])))
        
    #### Set axis range and titles
    ax.semilogx()
    ax.set_xlim(lrange[0],lrange[1])
    ax.set_ylim(-0.05*max(obs_sed),1.1*max(fobs))
    ax.set_xlabel(r'$\lambda%s$ [$\AA$]' %(rest_label))
    ax.set_ylabel(r'$f_\lambda$')
    
    ##### P(z)
    if (pz is not None) & (show_pz):            
        axp.plot(zgrid, pz, linewidth=1.0, color='orange',alpha=alph)
        axp.fill_between(zgrid,pz,np.zeros(zgrid.size),color='yellow')

        if zout['z_spec'][qz[idx]] > 0:
            axp.plot(zout['z_spec'][qz[idx]]*np.ones(2), np.array([0,1e6]),color='red',alpha=0.4)

        #### Set axis range and titles
        axp.set_xlim(0,np.ceil(np.max(zgrid)))
        axp.set_ylim(0,1.1*max(pz))
        axp.set_xlabel(r'$z$')
        axp.set_ylabel(r'$p(z)$')
        
    if (writePNG is not False) & (axes is None):
        if isinstance(writePNG, str):
            out=writePNG
        else:
            out='/tmp/test.pdf'
            
        fig.savefig(out,dpi=100)

    if axes is None:
        return fig #[ax, axp]