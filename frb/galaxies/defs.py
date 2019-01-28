""" Define allowed quantities for FRB galaxies

  Uncertainty is valid for any quantity with '_err' add-on, eg. W1_err
  Am also likely to add _flg for each as well
"""

##############################################################
# Redshift

valid_z = [
    'z',       # Preferred redshift, may derived from one of several ways
    'z_phot',  # Photometric redshift
    'z_spec',  # Spectroscopic redshift
]

##############################################################
# Photometry

# Filters
valid_filters = []

# DES
DES_bands = ['g', 'r', 'i', 'z', 'Y']
for band in DES_bands:
    valid_filters.append('DES_{:s}'.format(band))

# WISE
WISE_bands = ['W1', 'W2', 'W3', 'W4']
for band in WISE_bands:
    valid_filters.append('{:s}'.format(band))

##############################################################
# Line measurements

valid_neb_lines = [
    'Ha',  # Halpha flux erg/s/cm^2; pPXF
    'Hb',  # Hbeta flux erg/s/cm^2; pPXF
    'Hg',  # Hgamma flux erg/s/cm^2; pPXF
    '[NII] 6583',  # [NII] 6583 flux erg/s/cm^2; pPXF
    '[OII] 3726',  # [OII] flux erg/s/cm^2; pPXF
    '[OII] 3729',  # [OII] flux erg/s/cm^2; pPXF
    '[OIII] 5007',  # [OII] 5007 flux erg/s/cm^2; pPXF
]


##############################################################
# Derived quantities

valid_derived_photom = [
    'Mstar',           # Stellar mass; linear in Msun CIGALE
    'f_AGN',           # Fraction of AGN contribution to light; CIGALE
    'u-r',             # Rest-frame; CIGALE
    'Lnu_r',           # Specific luminosity (J/s/Hz); CIGALE; cosmology dependent
    'M_r',             # Absolute magnitude, r-band rest-frame; CIGALE+
    'SFR_photom',      # SFR in Msun/yr from photometry; CIGALE
    'EBV_photom',      # E(B-V) from photometry; CIGALE
    'Z_photom',        # Metallicity from photometry; CIGALE
    ]

valid_derived_nebular = [
    'EBV_nebular',     # E(B-V) from nebular line analysis (e.g. Ha/Hb)
    'SFR_nebular',     # SFR in Msun/yr from nebular emission (e.g. Halpha); pPXF+
    ]
