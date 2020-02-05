# PACKAGE IMPORT
import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import matplotlib
import matplotlib.pyplot as plt
import ne2001.density
from frb import igm
from pdf_defs import make_pdf, rv_amount
import suftware as sw
from astropy.cosmology import WMAP9 as cosmo
import random

# z to Mpc
z_stepsize = 0.00001
z_max = 0.1
z = np.arange(0,z_max,z_stepsize)
Mpc = cosmo.comoving_distance(z).value 
Mpc_val = 20.

def Mpc_from_z(z,Mpc,Mpc_val):
    ordering = np.argsort(z)
    dist_interpol = sp.interpolate.interp1d(Mpc[ordering], z[ordering], bounds_error=False, fill_value='extrapolate')
    Mpc_val_ = dist_interpol(Mpc_val)
    return Mpc_val_

z_DMcos = Mpc_from_z(z,Mpc,Mpc_val)

dm_ave_fn = igm.average_dm(z_DMcos, cosmo=None, cumul=False, neval=1000, mu=1.3)


print(Mpc_from_z(z,Mpc,Mpc_val))
print(dm_ave_fn)
