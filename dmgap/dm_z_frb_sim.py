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
from astropy import units
import random

# Data upload
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')

z_stepsize =  0.001
z_max = 3
z_grid = np.arange(0,z_max,z_stepsize)

dm_frb = np.array(frbcat_df['dmdiff'])* units.pc/units.cm**3

z = []
for i in range(len(dm_frb)):
    z_ = igm.z_from_DM(dm_frb[i], corr_nuisance=False)
    z = np.append(z,z_)

# Make PDF using DEFT
z_density = sw.DensityEstimator(z, alpha=3, bounding_box=[0,3])
z_func_ = z_density.evaluate(z_grid)
z_func = z_func_/np.sum(z_func_)
z = make_pdf(distribution=z_func,num_of_draws=rv_amount,grid=z_grid,stepsize=z_stepsize)
z_draws = z[0]
z_distribution = z[1]

# plt.hist(z,bins=100)
# plt.plot(z_grid,z_func,color='#f46036',label=r'$z$(DM)')
# plt.show()

# Calculate average dm
dm_ave_fn = igm.average_DM(z_max, cosmo=None, cumul=True, neval=1000, mu=1.3)

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

matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
plt.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=14)
plt.ylabel('PDF',fontsize=14)
plt.xlim(0,2000)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('sim_outputs/DM_FRB_z.png', dpi=300)
plt.show()

np.save('sim_outputs/dm_frb_sim.npy', dm_frb_sim)
np.save('sim_outputs/dm_grid.npy', dm_grid)
print('FRB PDF simulation saved')


