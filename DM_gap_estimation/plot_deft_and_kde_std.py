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

random.seed(1919)

"KDE: load datasets"
kde_100 = np.load('kde_std_data/kde_100.npy')
kde_1000 = np.load('kde_std_data/kde_1000.npy')
kde_observed = np.load('kde_std_data/kde_observed.npy')
kde_std_100 = np.load('kde_std_data/std_100.npy')
kde_std_1000 = np.load('kde_std_data/std_1000.npy')
kde_std_observed = np.load('kde_std_data/std_observed.npy')

"DEFT: create small datasets"
deft_sample_100 = np.asarray(random.sample(list(DM_frb),100))
deft_sample_1000  = np.asarray(random.sample(list(DM_frb),1000))
deft_density_100 = sw.DensityEstimator(deft_sample_100, alpha=4)
deft_density_1000 = sw.DensityEstimator(deft_sample_1000, alpha=4)
deft_density_observed = sw.DensityEstimator(DM_known_draws, alpha=4)

# Evaluate optimal density
deft_optimal_100 = deft_density_100.evaluate(DM_grid)
deft_optimal_1000 = deft_density_1000.evaluate(DM_grid)
deft_optimal_observed = deft_density_observed.evaluate(DM_grid)

# Evaluate sampled densities
deft_sampled_100 = deft_density_100.evaluate_samples(DM_grid)
deft_optimal_1000 = deft_density_1000.evaluate_samples(DM_grid)
deft_sampled_observed = deft_density_observed.evaluate_samples(DM_grid)

deft_std_100 = []
for i in range(len(DM_grid)):
    deft_std_100_ = np.std(deft_sampled_100[i])
    deft_std_100 = np.append(deft_std_100,deft_std_100_)

deft_std_1000 = []
for i in range(len(DM_grid)):
    deft_std_1000_ = np.std(deft_optimal_1000[i])
    deft_std_1000 = np.append(deft_std_1000,deft_std_1000_)

deft_std_observed = []
for i in range(len(DM_grid)):
    deft_std_observed_ = np.std(deft_sampled_observed[i])
    deft_std_observed = np.append(deft_std_observed,deft_std_observed_)

"100 draws plot"

plt.figure(figsize=[20,15])
plt.plot(DM_grid, deft_optimal_100, color='blue',linewidth=1,label='DEFT 100 samples')
plt.plot(DM_grid, kde_100, color='purple',linewidth=1,label='KDE 100 samples')
plt.fill_between(DM_grid, deft_optimal_100-deft_std_100, deft_optimal_100+deft_std_100, alpha=0.2,color='blue')
plt.fill_between(DM_grid, kde_100-kde_std_100, kde_100+kde_std_100, alpha=0.2,color='purple')
plt.xlim(0,1500)
plt.legend(fontsize=40)
plt.xlabel('$DM$')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('bootstrap_outputs/DEFT_vs_KDE/std_100.png')
plt.show()

"1000 draws plot"

plt.figure(figsize=[20,15])
plt.plot(DM_grid, deft_optimal_1000, color='blue',linewidth=1,label='DEFT 1000 samples')
plt.plot(DM_grid, kde_1000, color='purple',linewidth=1,label='KDE 1000 samples')
plt.fill_between(DM_grid, deft_optimal_1000-deft_std_1000, deft_optimal_1000+deft_std_1000, alpha=0.2,color='blue')
plt.fill_between(DM_grid, kde_1000-kde_std_100, kde_1000+kde_std_100, alpha=0.2,color='purple')
plt.xlim(0,1500)
plt.legend(fontsize=40)
plt.xlabel('$DM$')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('bootstrap_outputs/DEFT_vs_KDE/std_1000.png')
plt.show()

"Observed plot"

plt.figure(figsize=[20,15])
plt.plot(DM_grid, deft_optimal_observed, color='blue',linewidth=1,label='DEFT Observed')
plt.plot(DM_grid, kde_observed, color='purple',linewidth=1,label='KDE Observed')
plt.fill_between(DM_grid, deft_optimal_observed-deft_std_observed, deft_optimal_observed+deft_std_observed, alpha=0.2,color='blue')
plt.fill_between(DM_grid, kde_observed-kde_std_observed, kde_observed+kde_std_observed, alpha=0.2,color='purple')
plt.xlim(0,1500)
plt.legend(fontsize=40)
plt.xlabel('$DM$')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('bootstrap_outputs/DEFT_vs_KDE/std_obs.png')
plt.show()
