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
    catfile = os.path.join(input_dir, '{}.cat'.format(name))
    param_file = os.path.join(input_dir, 'zphot.param.{}'.format(name))
    #
    return catfile, param_file

def eazy_input_files(photom, input_dir, name, out_dir, prior_filter=None):
    """
    Write to disk a series of files needed to run EAZY
      - catalog file
      - translation file
      - param file

    Args:
        photom (dict):
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
    catfile, param_file = eazy_filenames(input_dir, name)

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
    outfile = os.path.join(input_dir, 'zphot.translate')
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
    _, param_file = eazy_filenames(input_dir, name)

    # Find the eazy executable
    path_to_eazy = spawn.find_executable('eazy')
    if path_to_eazy is None:
        raise ValueError("You must have eazy in your Unix path..")
    # Run it!
    command_line = [path_to_eazy, '-p', os.path.basename(param_file)]
    with open(logfile, 'w') as f:
        retval = subprocess.call(command_line, stdout=f, stderr=f, cwd=input_dir)


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

