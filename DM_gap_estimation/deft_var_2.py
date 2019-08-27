# PACKAGE IMPORT
import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import random
import matplotlib
from matplotlib import rcParams
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from DM_definitions import DM_known_draws, make_pdf, make_kde_funtion
from DM_kde_gen import DM_grid, DM_frb, DM_rand_draws
import suftware as sw
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

matplotlib.rc('xtick', labelsize=30) 
matplotlib.rc('ytick', labelsize=30) 
# rcParams['figure.figsize'] = 15,8

DM_kde_rand_100 = np.load('variance_plots/100kde.npy')
DM_kde_rand_1000 = np.load('variance_plots/1000kde.npy')
DM_kde_obs = np.load('variance_plots/obskde.npy')
var_by_dm_100 = np.load('variance_plots/var100kde.npy')
var_by_dm_1000 = np.load('variance_plots/var1000kde.npy')
var_by_dm_obs = np.load('variance_plots/varobskde.npy')

DM_rand_draws_sample_100 = np.asarray(random.sample(list(DM_frb),100))
DM_rand_draws_sample_1000 = np.asarray(random.sample(list(DM_frb),1000))
density_100 = sw.DensityEstimator(DM_rand_draws_sample_100, alpha=4)
density_1000 = sw.DensityEstimator(DM_rand_draws_sample_1000, alpha=4)
density_obs = sw.DensityEstimator(DM_known_draws, alpha=4)

# Evaluate optimal density
new_values_100 = density_100.evaluate(DM_grid)
new_values_1000 = density_1000.evaluate(DM_grid)
new_values_obs = density_obs.evaluate(DM_grid)

# Evaluate sampled densities
new_sampled_values_100 = density_100.evaluate_samples(DM_grid)
new_sampled_values_1000 = density_1000.evaluate_samples(DM_grid)
new_sampled_values_obs = density_obs.evaluate_samples(DM_grid)


var100 = []
for i in range(len(DM_grid)):
    var100_ = np.std(new_sampled_values_100[i])
    var100 = np.append(var100,var100_)

var1000 = []
for i in range(len(DM_grid)):
    var1000_ = np.std(new_sampled_values_1000[i])
    var1000 = np.append(var1000,var1000_)

varobs = []
for i in range(len(DM_grid)):
    varobs_ = np.std(new_sampled_values_obs[i])
    varobs = np.append(varobs,varobs_)
   

"100 draws plot"

plt.figure(figsize=[20,15])
plt.plot(DM_grid, new_values_100, color='blue',linewidth=1,label='DEFT 100 samples')
plt.plot(DM_grid, DM_kde_rand_100, color='purple',linewidth=1,label='KDE 100 samples')
plt.fill_between(DM_grid, new_values_100-var100, new_values_100+var100, alpha=0.2,color='blue')
plt.fill_between(DM_grid, DM_kde_rand_100-var_by_dm_100, DM_kde_rand_100+var_by_dm_100, alpha=0.2,color='purple')
plt.xlim(0,1500)
plt.legend(fontsize=40)
plt.xlabel('$DM$')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('bootstrap_outputs/DEFT/variances_100.png')
plt.show()

"1000 draws plot"

plt.figure(figsize=[20,15])
plt.plot(DM_grid, new_values_1000, color='blue',linewidth=1,label='DEFT 1000 samples')
plt.plot(DM_grid, DM_kde_rand_1000, color='purple',linewidth=1,label='KDE 1000 samples')
plt.fill_between(DM_grid, new_values_1000-var1000, new_values_1000+var1000, alpha=0.2,color='blue')
plt.fill_between(DM_grid, DM_kde_rand_1000-var_by_dm_1000, DM_kde_rand_1000+var_by_dm_1000, alpha=0.2,color='purple')
plt.xlim(0,1500)
plt.legend(fontsize=40)
plt.xlabel('$DM$')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('bootstrap_outputs/DEFT/variances_1000.png')
plt.show()

"Observed plot"

plt.figure(figsize=[20,15])
plt.plot(DM_grid, new_values_obs, color='blue',linewidth=1,label='DEFT Observed')
plt.plot(DM_grid, DM_kde_obs, color='purple',linewidth=1,label='KDE Observed')
plt.fill_between(DM_grid, new_values_obs-varobs, new_values_obs+varobs, alpha=0.2,color='blue')
plt.fill_between(DM_grid, DM_kde_obs-var_by_dm_obs, DM_kde_obs+var_by_dm_obs, alpha=0.2,color='purple')
plt.xlim(0,1500)
plt.legend(fontsize=40)
plt.xlabel('$DM$')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('bootstrap_outputs/DEFT/variances_obs.png')
plt.show()

# plt.figure(figsize=[20,15])
# plt.plot(DM_grid, new_values_obs, color='blue',linewidth=1,label='DEFT observed samples')
# plt.plot(DM_grid, new_values_100, color='blue',linewidth=1,label='DEFT 100 samples')
# plt.plot(DM_grid, new_values_1000, color='purple',linewidth=1,label='DEFT 1000 samples')
# plt.plot(DM_grid,varobs, linewidth=1.5, label='Observed',color='k',linestyle='dashed')
# plt.plot(DM_grid,var100, linewidth=1.5, label='100 samples', color='purple')
# plt.plot(DM_grid,var1000, linewidth=1.5, label='1000 samples',color='green')
# plt.xlim(0,2000)
# plt.legend(fontsize=40)
# plt.xlabel('$DM$')
# plt.ylabel('PDF')
# plt.tight_layout()
# plt.savefig('bootstrap_outputs/DEFT/variances_plot.png')
# plt.show()



# plt.figure(figsize=[20,15])
# # Plot optimal and posterior-sampled densities
# plt.plot(DM_grid, new_sampled_values_1000, color='dodgerblue', alpha=.1)
# plt.plot(DM_grid, new_values_1000, color='blue')
# plt.xlabel('$DM$', fontsize=30)
# plt.ylabel('PDF', fontsize=30)
# plt.xlim(0,200)
# # plt.ylim(0,0.001)
# plt.tight_layout()
# plt.savefig('bootstrap_outputs/sim_DEFT_cropped_1000.png')
# plt.show()



"Plot DM_FRB rand draws"
# plt.plot(DM_grid, DM_kde_rand_distribution/np.sum(DM_kde_rand_distribution), linewidth=1, alpha=1, label=r'DM KDE', color='purple')
# plt.hist(DM_kde_draws_rand, density=True, bins=1000, histtype='stepfilled', alpha=0.5, label=r'DM draws from KDE',color='purple')
# plt.hist(DM_frb, density=True, bins=500, histtype='stepfilled', alpha=0.5,label=r'DM simulated')
# # for xc in DM_rand_draws:
# #     plt.scatter(xc,0,marker='o',s=2,c='k')
# plt.xlabel('$DM$', fontsize=30)
# plt.ylabel('PDF', fontsize=30)
# plt.legend(fontsize=14)
# plt.xlim(0,2500)
# plt.savefig('DM_outputs/DM_kde_choice.png')
# plt.show()



"Seaborn plot"
# sns.distplot( DM_known_draws, hist = True, kde = True, rug = True, bins=100,
#              color = 'darkblue', 
#              kde_kws={'linewidth': 2, 'bw':99},
#              rug_kws={'color': 'black'})
# plt.xlabel('$DM_{Observed}$', fontsize=28)
# plt.ylabel('PDF', fontsize=28)
# plt.legend(fontsize=28)
# plt.xlim(0,2500)
# plt.tight_layout()
# plt.savefig('DM_outputs/DM_kde_known_sns.png')
# plt.show()

# "Seaborn plot"
# sns.distplot(DM_rand_draws, hist = True, kde = True, rug = True, bins=100,
#              color = 'darkblue', 
#              kde_kws={'linewidth': 2, 'bw':25},
#              rug_kws={'color': 'black'})
# plt.xlabel('$DM_{Simulated}$', fontsize=18)
# plt.ylabel('PDF', fontsize=18)
# plt.xlim(0,1500)
# # plt.ylim(0,0.01)
# plt.tight_layout()
# plt.savefig('DM_outputs/DM_kde_choice_sns_1000.png')
# plt.show()

# print(kde_original)