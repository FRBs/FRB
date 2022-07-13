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
NOT_bands = ['u', 'g','r','i','z']
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

#MMT
MMIRS_bands = ['J','H','K']
for band in MMIRS_bands:
    valid_filters.append('MMIRS_{:s}'.format(band))

#2MASS
MASS_bands = ['J','H','K']
for band in MASS_bands:
    valid_filters.append('2MASS_{:s}'.format(band))

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

NSC_bands = ['u','g', 'r', 'i', 'z', 'Y', 'VR']
for band in NSC_bands:
    valid_filters.append('NSC_{:s}'.format(band))

DECam_bands = ['u','g', 'r', 'i', 'z', 'Y', 'VR']
for band in DECam_bands:
    valid_filters.append("DECam_{:s}".format(band))


# For upper limits, the flux is 3sigma and the error is set to -99.0
valid_flux = [entry+'_flux' for entry in valid_filters]
valid_ref = [entry+'_ref' for entry in valid_filters]

valid_photom = valid_filters + ['EBV']  # Galactic

##############################################################
# Line measurements -- Use ppxf naming only!
#   Not corrected for internal dust extinction

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
    '[SII] 6716',  # [SII] 6716 flux erg/s/cm^2; pPXF
    '[SII] 6731',  # [SII] 6731 flux erg/s/cm^2; pPXF
]
valid_neb_ref = [entry+'_ref' for entry in valid_neb_lines]

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
    'SFR_radio',       # SFR in Msun/yr from radio photometry
    'EBV_photom',      # E(B-V) from photometry; CIGALE
    'EBV_spec',        # E(B-V) from spectral SED; pPXF
    'Z_photom',        # Metallicity from photometry; CIGALE
    'Z_spec',          # Metallicity from spectra; pPXF
    ]

valid_derived_nebular = [
    'AV_nebular',      # AV from nebular line analysis (e.g. Ha/Hb)
    'SFR_nebular',     # SFR in Msun/yr from nebular emission (e.g. Halpha); pPXF+
    ]

# Local measurements (usually from high-spatial resolution imaging)
valid_derived_local = [
    'UV_filter', # name of UV filter
    'IR_filter', # name of UV filter
    'UV_SB', # comment='UV SB in UV_filter at FRB position')
    'IR_SB', # comment='IR SB in IR_filter at FRB position')
    'UV_global', # comment=UV photometry of entire galaxy
    'IR_global', # comment=IR photometry of entire galaxy
    'SSFR', # comment='Surface density of Star formation at FRB position (Msun / yr / kpc^2)')
    'SMStar', # comment='Surface density of Stars at FRB position')
    'reff_iso', # comment='Effective Isophotal radius (angular)')
    'UVff', # comment='Fractional flux at FRB location in UV_filter')
    'IRff', # comment='Fractional flux at FRB location in IR_filter')
    'IRfe', # comment='Enclosed flux at FRB location in IR_filter')
    'IRlim', # comment='Limiting IR flux (mag) at FRB location after subtracting off galfit model.')
]


valid_derived = valid_derived_photom + valid_derived_nebular + valid_derived_local
valid_derived_ref = [entry+'_ref' for entry in valid_derived]
