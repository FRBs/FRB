# PACKAGE IMPORT
import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import random
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
from percentile_defs import find_quantile

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 
plt.figure(figsize=(12,8))

# Data upload 
dm_frb_sim = np.load('sim_outputs/dm_frb_sim.npy')
dm_grid = np.load('sim_outputs/dm_grid.npy')

deft_optimal_frb_len = np.load('kde_and_deft_data/deft_frb_len_optimal.npy')
deft_optimal_1000 = np.load('kde_and_deft_data/deft_1000_optimal.npy')
deft_optimal_2000 = np.load('kde_and_deft_data/deft_2000_optimal.npy')
deft_optimal_4000 = np.load('kde_and_deft_data/deft_4000_optimal.npy')
deft_optimal_10000 = np.load('kde_and_deft_data/deft_10000_optimal.npy')

deft_sampled_frb_len = np.load('kde_and_deft_data/deft_frb_len_sampled.npy')
deft_sampled_1000 = np.load('kde_and_deft_data/deft_1000_sampled.npy')
deft_sampled_2000 = np.load('kde_and_deft_data/deft_2000_sampled.npy')
deft_sampled_4000 = np.load('kde_and_deft_data/deft_4000_sampled.npy')
deft_sampled_10000 = np.load('kde_and_deft_data/deft_10000_sampled.npy')

# DEFT with 10000 draws
plt.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
plt.plot(dm_grid,deft_optimal_10000,color='crimson',label=r'DEFT DM$_\mathrm{FRB}$ (n=10000)')
plt.plot(dm_grid,deft_sampled_10000,color='crimson',linewidth=0.5,alpha=.1)
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
plt.ylabel('PDF',fontsize=22)
plt.xlim(0,2000)
plt.legend(fontsize=22)
plt.tight_layout()
plt.savefig('sim_outputs/DM_10000_samples.png', dpi=300)
plt.show()