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

# SDSS
SDSS_bands = ['u', 'g', 'r', 'i', 'z']
for band in SDSS_bands:
    valid_filters.append('SDSS_{:s}'.format(band))
    
# DES
DES_bands = ['g', 'r', 'i', 'z', 'Y']
for band in DES_bands:
    valid_filters.append('DES_{:s}'.format(band))

# DECaLS
DECaL_bands = ['g', 'r', 'z']
for band in DECaL_bands:
    valid_filters.append('DECaL_{:s}'.format(band))

#PanSTARRS
PanSTARRS_bands = ['g','r','i','z','y']
for band in PanSTARRS_bands:
    valid_filters.append('Pan-STARRS_{:s}'.format(band))

# VLT
VLT_bands = ['u', 'g', 'I', 'z']
for band in VLT_bands:
    valid_filters.append('VLT_FORS2_{:s}'.format(band))

# GMOS
#south
GMOS_bands = ['u', 'g', 'r', 'i', 'z']
for band in GMOS_bands:
    valid_filters.append('GMOS_S_{:s}'.format(band))
#north
for band in GMOS_bands:
    valid_filters.append('GMOS_N_{:s}'.format(band))

#NOT
NOT_bands = ['g','r','i','z']
for band in NOT_bands:
    valid_filters.append('NOT_{:s}'.format(band))

#NIRI
NIRI_bands = ['J']
for band in NIRI_bands:
    valid_filters.append('NIRI_{:s}'.format(band))

#LRIS
LRISb_bands = ['U', 'G', 'V', 'B']
for band in LRISb_bands:
    valid_filters.append('LRISb_{:s}'.format(band))

LRISr_bands = ['V', 'R', 'I']
for band in LRISr_bands:
    valid_filters.append('LRISr_{:s}'.format(band))

# VISTA (VIRCAM)
VISTA_bands = ['Y','J','H','Ks']
for band in VISTA_bands:
    valid_filters.append('VISTA_{:s}'.format(band))

# HST instruments
# WFC3
WFC3_bands = ['F300X', 'F110W', 'F160W', 'F763M']
for band in WFC3_bands:
    valid_filters.append('WFC3_{:s}'.format(band))

# WISE
WISE_bands = ['W1', 'W2', 'W3', 'W4']
for band in WISE_bands:
    valid_filters.append('WISE_{:s}'.format(band))

# Spitzer
Spitzer_bands = ['3.6', '4.5']
for band in Spitzer_bands:
    valid_filters.append('Spitzer_{:s}'.format(band))

valid_flux = [entry+'_flux' for entry in valid_filters]

valid_photom = valid_filters + ['EBV']  # Galactic

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
    'ra',         # RA centroid inferred from Galfit
    'dec',        # DEC centroid inferred from Galfit
    'n',          # Sersic index from Galfit
]

##############################################################
# Offsets
valid_offsets = [
    'ang_best',   # Angular offset in arcsec from localization centroid to galaxy
    'ang_avg',    # Angular offset in arcsec averaging over localization
    'physical',   # Physical offset in kpc;  Uses ang_best
    ]


##############################################################
# Positional  (Astrometric and Source) Errors
valid_positional_error = [
    'ra_astrometric',   # error for astrometric tie in RA; arcsec
    'dec_astrometric',  # error for astrometric tie in Dec; arcsec
    'ra_source',        # RA error for source position (e.g. from source extractor); arcsec
    'dec_source',       # Dec error for source position; arcsec
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
    'age_mass',        # Age weighted mass from CIGALE
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
