""" Methods related to nebular line analysis, e.g. dust extinction, SFR"""
import pdb

from pkg_resources import resource_filename

import numpy as np

from scipy.interpolate import interp1d

from astropy.table import Table

def calc_dust_extinct(neb_lines, method, curve='MW'):

    # Load the extinction curve
    if curve == 'MW':
        dust_file = resource_filename('frb', 'data/Dust/MW_dust.dat')
        MW_dust = Table.read(dust_file, format='ascii')
    else:
        raise IOError("Not ready for this extinction curve!")

    # Generate function for interpolating
    alAV = interp1d(MW_dust['wave'], MW_dust['Al_AV'])

    # Which ratio?
    if method == 'Ha/Hb':
        Ha_Hb_intrin = 2.8 # Osterbrock
        wave1 = 6564.6  # redder
        wave2 = 4862.7
        #
        F1_obs = neb_lines['Ha']
        F2_obs = neb_lines['Hb']
        #
        pair = True
    else:
        print("Not prepared for this method of analysis: {}".format(method))
        raise IOError("See docs for available options")

    if not pair:
        raise IOError("Not ready for this mode")

    #
    a1AV = alAV(wave1)
    a2AV = alAV(wave2)

    fratio_obs = F1_obs/F2_obs

    AV = 2.5 * np.log10(Ha_Hb_intrin/fratio_obs) / (a1AV - a2AV)

    EBV = AV / 3.1

    return EBV
