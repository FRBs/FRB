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
    'z_FRB',   # FRB redshift
]

##############################################################
# Error Ellipse
valid_e = [
    'a',       # Major axis
    'b',       # Minor axis
    'theta',   # Rotation of the major axis E from N (deg)
    'cl',      # Confidence level of the ellipse
    ]

##############################################################
# Photometry

# Filters
valid_filters = []

# DES
DES_bands = ['g', 'r', 'i', 'z', 'Y']
for band in DES_bands:
    valid_filters.append('DES_{:s}'.format(band))

# VLT
VLT_bands = ['g', 'I']
for band in VLT_bands:
    valid_filters.append('VLT_{:s}'.format(band))

# WISE
WISE_bands = ['W1', 'W2', 'W3', 'W4']
for band in WISE_bands:
    valid_filters.append('{:s}'.format(band))
    
valid_photom = valid_filters

##############################################################
# Line measurements -- Use linetools naming only!!!

valid_neb_lines = [
    'Halpha',  # Halpha flux erg/s/cm^2; pPXF
    'Hbeta',  # Hbeta flux erg/s/cm^2; pPXF
    'Hgamma',  # Hgamma flux erg/s/cm^2; pPXF
    '[NII] 6548',  # [NII] 6584 flux erg/s/cm^2; 
    '[NII] 6584',  # [NII] 6584 flux erg/s/cm^2; pPXF
    '[OII] 3726',  # [OII] flux erg/s/cm^2; pPXF
    '[OII] 3729',  # [OII] flux erg/s/cm^2; pPXF
    '[OIII] 4959',  # [OII] 4959 flux erg/s/cm^2; 
    '[OIII] 5007',  # [OII] 5007 flux erg/s/cm^2; pPXF
]

##############################################################
# Morphology

valid_morphology = [
    'reff_ang',   # Effective radius in arcsec; Galfit
    'reff_kpc',   # Effective radius in kpc; Galfit
    'n',          # Sersic index; Galfit
    'PA',         # Position angle (deg); Galfit
    'b/a',        # Ellipticity; Galfit
]

##############################################################
# Derived quantities

valid_derived_photom = [
    'Mstar',           # Stellar mass; linear in Msun CIGALE
    'Mstar_spec',      # Stellar mass from pPXF; linear in Msun
    'f_AGN',           # Fraction of AGN contribution to light; CIGALE
    'u-r',             # Rest-frame; CIGALE
    'Lnu_r',           # Specific luminosity (J/s/Hz); CIGALE; cosmology dependent
    'M_r',             # Absolute magnitude, r-band rest-frame; CIGALE+
    'SFR_photom',      # SFR in Msun/yr from photometry; CIGALE
    'EBV_photom',      # E(B-V) from photometry; CIGALE
    'EBV_spec',        # E(B-V) from spectral SED; pPXF
    'Z_photom',        # Metallicity from photometry; CIGALE
    'Z_spec',          # Metallicity from spectra; pPXF
    ]

valid_derived_nebular = [
    'AV_nebular',      # AV from nebular line analysis (e.g. Ha/Hb)
    'SFR_nebular',     # SFR in Msun/yr from nebular emission (e.g. Halpha); pPXF+
    ]

valid_derived = valid_derived_photom + valid_derived_nebular
