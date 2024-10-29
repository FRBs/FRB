*****************************
Magnitudes to flux conversion
*****************************

Most public photometry tables report brightness of objects
in terms of magnitudes as opposed to flux units (Jy or mJy).
However, using software such as CIGALE and EAZY is made
less ambiguous if fluxes are supplied. To enable this, we provide
the function `frb.surveys.catalog_utils.convert_mags_to_flux`

The definition for magnitude errors arises from the following
equation:

.. math::
  
  m\pm \delta m = ZP-2.5 \log_{10}(f\pm \delta f) 
  \Rightarrow m = ZP-2.5 \log_{10}(f)
  \Rightarrow \delta m = 2.5 \log_{10}(1+\delta f/f)

While the math:`\delta m` formula above is exact, most photometry
software use a first order Taylor expansion in math:`\delta f/f`
assuming the error is typically small. i.e. 

.. math::
  
  \delta m \approx \frac{2.5}{\log_{10}}(\delta f/f)

This will obviously break down when looking at really faint objects,
barely at the threshold of detection above the noise floor.
`convert_mags_to_flux` explicitly assumes the catalog mag errors are
produced using this approximation. However, if you have reason to believe
the exact formula was used, you can set `exact_mag_err=True` to get the
correct flux error.