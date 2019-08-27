# PACKAGE IMPORT
import numpy as np
import scipy as sp
from scipy.stats import norm, lognorm
from scipy import signal
import pandas as pd
from pandas import DataFrame as df
import sklearn
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('agg')
import ne2001.density
from frb import igm
from frb import mw
from astropy.cosmology import WMAP9 as cosmo
import random
from DM_definitions import rv_amount, z_known_draws, FRB_distribution_from_SFR, make_pdf, make_kde_funtion

import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 

"VARIABLES"
z_stepsize =  0.001
z_max = 1.5
z_grid = np.arange(0,z_max,z_stepsize)

"Find distribution of FRBs along z (going like the SFR)"
z_FRB_likelihood_SFR = FRB_distribution_from_SFR(z_grid=z_grid)
z_SFR_function = make_pdf(distribution=z_FRB_likelihood_SFR,num_of_draws=rv_amount,grid=z_grid[1:],stepsize=z_stepsize) #pdf from z_FRB_likelihood_SFR
z_SFR_draws = z_SFR_function[0]
z_SFR_distribution  = z_SFR_function[1]

"Make KDE distribution with gaussian kernel and known z values, then use this distribution to make a PDF"
z_kde_func_ = make_kde_funtion(grid=z_grid[1:], draws=z_known_draws, min_bandwidth=0.05, max_bandwidth=0.5, bandwidth_stepsize=0.05, cv=3, kernel='gaussian')
z_kde_func = z_kde_func_/np.sum(z_kde_func_) #normalised
z_kde = make_pdf(distribution=z_kde_func,num_of_draws=rv_amount,grid=z_grid[1:],stepsize=z_stepsize)
z_draws = z_kde[0]
z_distribution = z_kde[1] #z_kde_func scaled

"Calculate average DM using X's frb package"
DM_ave_fn = igm.average_DM(z_max, cosmo=None, cumul=True, neval=1000, mu=1.3)
DM_ave_ = DM_ave_fn[0].value
z_ave = DM_ave_fn[1]
DM_int_ave = sp.interpolate.CubicSpline(z_ave,DM_ave_) #spline
DM_ave = DM_int_ave(z_draws)

"Find DM_cos"
F = 0.2
sigma_DM = np.minimum(F*z_draws**(-0.5),0.5)
rand_draws_gauss = np.random.normal(loc=0, scale=1, size=rv_amount)
rand_draws_gauss = np.maximum(rand_draws_gauss,-1)
rand_draws_gauss = np.minimum(rand_draws_gauss,1)
DM_cos = np.maximum(DM_ave + rand_draws_gauss*DM_ave*sigma_DM,0)

"DM_host"
#use loglikelihood function for this
mu_hostDM = 40.
sigma_hostDM = 1.
DM_host = lognorm(sigma_hostDM, loc=mu_hostDM, scale=15).rvs(size=rv_amount)
params_host = lognorm.fit(DM_host)
x_host_grid = np.arange(0,max(DM_host)) #redshift range of interest
pdf_host_fitted = lognorm.pdf(x_host_grid, params_host[0], loc=params_host[1], scale=params_host[2])
x_host_fitted = np.arange(len(pdf_host_fitted))

"DM_MWhalo"
#use delta function for this
DM_MWhalo = 60.

"DM_FRB"
DM_frb = DM_MWhalo+DM_host+DM_cos

"SAVE DM_FRB SIULATION TO FILE"
np.save('DM_outputs/DM_frb.npy', DM_frb)

"PLOTS"

"FRB_distribution_SFR"
# plt.hist(z_SFR_draws, density=True, bins=100, histtype='stepfilled', alpha=0.5)
# plt.plot(z_grid[1:], z_SFR_distribution, alpha=0.5)
# plt.xlabel('FRB distribution (z)', fontsize=14)
# plt.ylabel('PDF', fontsize=14)
# plt.title(r'$z_{FRB}$ PDF')
# plt.show()

"z_KDE"
# fig = plt.figure()
# plt.hist(z_SFR_draws, density=True, bins=500, histtype='stepfilled', alpha=0.5, label=r'$z_{SFR}$')
# plt.hist(z_known_draws, density=True, bins=len(z_known_draws)+1, histtype='stepfilled', alpha=0.5, label=r'$z_{FRB, observed}$',color='r')
# plt.plot(z_grid[1:], z_kde[1], linewidth=1, alpha=1, label=r'KDE $z_{FRB, observed}$', color='r')
# plt.plot(z_grid[1:], z_SFR_distribution, linewidth=1, alpha=1, label=r'KDE $z_{FRB, SFR}$', color='blue')
# for xc in z_known_draws:
#     plt.axvline(x=xc,linewidth=1,color='k',linestyle='--')
# plt.xlabel('$z$', fontsize=14)
# plt.ylabel('PDF', fontsize=14)
# plt.legend(fontsize=14)
# plt.show()
# plt.savefig('DM_outputs/z_PDF')

"Plot DM_ave"
# plt.plot(z_ave,DM_int_ave(z_ave),lw=2,label='Average DM Spline')
# plt.plot(z_ave, DM_ave_,lw=2,linestyle='dashed',label='Average DM')
# plt.xlabel('$z$', fontsize=14)
# plt.ylabel('$DM_{ave,cosmic}$', fontsize=14)
# plt.legend(loc=2, prop={'size': 10})
# plt.show()

"Plot DM_cos PDF"
# plt.hist(DM_cos, density=True, bins=1000, histtype='stepfilled', alpha=0.5, label='Gaussian',color='b')
# plt.xlabel(r'$DM_{cosmic}$',fontsize=14)
# plt.ylabel('PDF', fontsize=14)
# plt.legend(fontsize=14)
# plt.show()

"Plot DM_host PDF"
# plt.hist(DM_host, density=True, bins=1000, histtype='stepfilled', alpha=0.5)
# plt.plot(x_host_fitted,pdf_host_fitted, alpha=0.5)
# plt.xlabel('$DM_{host}$', fontsize=14)
# plt.ylabel('PDF', fontsize=14)
# plt.xlim(70,150)
# plt.show()