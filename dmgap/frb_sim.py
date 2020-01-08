# PACKAGE IMPORT
import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import matplotlib
import matplotlib.pyplot as plt
import ne2001.density
from frb import igm
from scipy.stats import lognorm
from pdf_defs import make_pdf, rv_amount
import suftware as sw
from astropy.cosmology import WMAP9 as cosmo
import random

# VARIABLES
z_obs = np.asarray([0.0337, 0.1927, 0.3212, 0.4755, 0.29, 0.66, 0.1178, 0.378])
z_stepsize =  0.001
z_max = 1.5
z_grid = np.arange(0,z_max,z_stepsize)
alpha_param = 2

# Make PDF using DEFT
z_density = sw.DensityEstimator(z_obs, alpha=alpha_param, bounding_box=[0,3])
z_func = z_density.evaluate(z_grid)
z_func = z_func/np.sum(z_func)
z = make_pdf(distribution=z_func,num_of_draws=rv_amount,grid=z_grid,stepsize=z_stepsize)
z_draws = z[0]
z_distribution = z[1]

# Calculate average dm
dm_ave_fn = igm.average_dm(z_max, cosmo=None, cumul=True, neval=1000, mu=1.3)

dm_ave_ = dm_ave_fn[0].value
z_ave = dm_ave_fn[1]
dm_int_ave = sp.interpolate.CubicSpline(z_ave,dm_ave_) #spline
dm_ave = dm_int_ave(z_draws)

# Find cosmic dm
F = 0.2
sigma_dm = np.minimum(F*z_draws**(-0.5),0.5)
rand_draws_gauss = np.random.normal(loc=0, scale=1, size=rv_amount)
rand_draws_gauss = np.maximum(rand_draws_gauss,-1)
rand_draws_gauss = np.minimum(rand_draws_gauss,1)
dm_cos = np.maximum(dm_ave + rand_draws_gauss*dm_ave*sigma_dm,0)

# Estimate halo
dm_halo  = 50.

# DM_host
#use loglikelihood function for this
mu_host = 40.
sigma_host = 1.
dm_host = lognorm(sigma_host, loc=mu_host, scale=15).rvs(size=rv_amount)

dm_frb_sim = dm_cos+dm_halo+dm_host
dm_stepsize =  0.01
dm_max = 3000
dm_grid = np.arange(0,dm_max,dm_stepsize)

np.save('sim_outputs/dm_frb_sim.npy', dm_frb_sim)
np.save('sim_outputs/dm_grid.npy', dm_grid)
print('FRB PDF simulation saved')


