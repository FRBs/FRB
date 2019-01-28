""" Define allowed quantities for FRB galaxies"""

##############################################################
# Photometry

# Filters

valid_filters = []  # Uncertainty is valid with '_err' add-on, eg. W1_err

# DES
DES_bands = ['g', 'r', 'i', 'z', 'Y']
for band in DES_bands:
    valid_filters.append('DES_{:s}'.format(band))

# WISE
WISE_bands = ['W1', 'W2', 'W3', 'W4']
for band in WISE_bands:
    valid_filters.append('{:s}'.format(band))

##############################################################
# Derived quantities

valid_derived_photom = [
    'Mstar',           # Stellar mass; linear in Msun CIGALE
    'f_AGN',           # Fraction of AGN contribution to light; CIGALE
    'u-r',             # Rest-frame; CIGALE
    'M_r',             # Absolute magnitude, r-band rest-frame; CIGALE+
    'SFR_photom',      # SFR in Msun/yr from photometry; CIGALE
    'EBV_photom',      # E(B-V) from photometry; CIGALE
    'Z_photom',        # Metallicity from photometry; CIGALE
    ]
