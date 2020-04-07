""" Module to faciliate scripting of EAZY analysis"""

import os
import warnings
from pkg_resources import resource_filename

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


def eazy_input_files(photom, input_dir, name, prior_filter=None):

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
    catfile = os.path.join(input_dir, '{}.cat'.format(name))
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

    new_default_file = os.path.join(input_dir, 'zphot.param.{}'.format(name))
    with open(new_default_file, 'w') as f:
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
            else:
                line = dfline
            # Write
            f.write(line)
    print("Wrote param file: {}".format(new_default_file))




def run_eazy_on_photom(photom):
    pass

