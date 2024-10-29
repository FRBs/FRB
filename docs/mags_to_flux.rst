*****************************
Magnitudes to flux conversion
*****************************

Most public photometry tables report brightness of objects
in terms of magnitudes as opposed to flux units (Jy or mJy).
However, using software such as CIGALE and EAZY are made
less ambiguous if fluxes are used. To enable this, we provide
the funciton `frb.surveys.catalog_utils.convert_mags_to_flux`