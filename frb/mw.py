""" Module for calculations related to the Galaxy
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

from frb.halos.models import ModifiedNFW

from ne2001 import density

from IPython import embed

def ismDM(coord):
    """
    Calculate the dispersion measure (DM) contribution from the interstellar medium (ISM)
    along a given line of sight in galactic coordinates.


    Parameters:
        coord (astropy.coordinates.SkyCoord): The input sky coordinate to transform 
            into galactic coordinates.


    Returns:
        float: The dispersion measure (DM) value calculated for the ISM along the 
        specified line of sight up to a distance of 100 parsecs.
    """
    gcoord = coord.transform_to('galactic')
    l, b = gcoord.l.value, gcoord.b.value
    
    ne = density.ElectronDensity()#**PARAMS)

    ismDM = ne.DM(l, b, 100.)
    
    # Return
    return ismDM

def haloDM(coord, f_diffuse=0.75, zero=True):
    """
    Uses a Modified NFW profile with parameters from Boylan-Kolchin et al. 2013.



    Args:
        coord (astropy.coordinates.SkyCoord): Sky coordinate for the line of sight.
        f_diffuse (float, optional): Fraction of the halo mass in diffuse gas. Default 0.75.
        zero (bool, optional): If True, zero out the inner 10 kpc. Default True.



    Returns:
        float: DM contribution from the MW halo in pc/cm^3.
    """
    gcoord = coord.transform_to('galactic')
    l, b = gcoord.l.value, gcoord.b.value

    # MW
    Mhalo = np.log10(1.5e12) # Boylan-Kolchin et al. 2013
    c = 7.7
    mnfw_2 = ModifiedNFW(log_Mhalo=Mhalo, f_hot=f_diffuse, y0=2, alpha=2, c=c)
    # Zero out inner 10kpc?
    if zero:
        mnfw_2.zero_inner_ne = 10.  # kpc
    params = dict(F=1., e_density=1.)
    model_ne = density.NEobject(mnfw_2.ne, **params)
    haloDM = model_ne.DM(l, b, mnfw_2.r200.value)
    #
    return haloDM
