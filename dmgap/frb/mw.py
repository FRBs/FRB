""" Module for calculations related to the Galaxy
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import numpy as np

from astropy import units
import warnings

from frb.halos import ModifiedNFW

from ne2001 import density

def ismdm(coord):
    gcoord = coord.transform_to('galactic')
    l, b = gcoord.l.value, gcoord.b.value
    
    ne = density.ElectronDensity()#**PARAMS)
    ismdm = ne.dm(l, b, 100.)
    
    # Return
    return ismdm

def halodm(coord, f_diffuse=0.75, zero=True):

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
    halodm = model_ne.dm(l, b, mnfw_2.r200.value)
    #
    return halodm
