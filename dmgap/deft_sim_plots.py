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
# import seaborn as sns
# sns.set()
# sns.set_style('darkgrid')

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.arange(0,6000,0.01)

dmdiff_frb_optimal = np.load('obs_outputs/dmdiff_frb_optimal.npy')
dmdiff_frb_ensemble = np.load('obs_outputs/dmdiff_frb_ensemble.npy')

dm_frb_sim = np.load('sim_outputs/dm_frb_sim.npy')
dm_grid = np.load('sim_outputs/dm_grid.npy')

deft_optimal_frb_len = np.load('kde_and_deft_data/deft_frb_len_optimal.npy')
deft_optimal_1000 = np.load('kde_and_deft_data/deft_1000_optimal.npy')
deft_optimal_5000 = np.load('kde_and_deft_data/deft_5000_optimal.npy')
deft_optimal_10000 = np.load('kde_and_deft_data/deft_10000_optimal.npy')

deft_sampled_frb_len = np.load('kde_and_deft_data/deft_frb_len_sampled.npy')
deft_sampled_1000 = np.load('kde_and_deft_data/deft_1000_sampled.npy')
deft_sampled_5000 = np.load('kde_and_deft_data/deft_5000_sampled.npy')
deft_sampled_10000 = np.load('kde_and_deft_data/deft_10000_sampled.npy')

deft_std_frb_len = np.load('kde_and_deft_data/deft_std_frb_len.npy')
deft_std_1000 = np.load('kde_and_deft_data/deft_std_1000.npy')
deft_std_5000 = np.load('kde_and_deft_data/deft_std_5000.npy')
deft_std_10000 = np.load('kde_and_deft_data/deft_std_10000.npy')

# Plots
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex='col', sharey='row', figsize=(8,12))
ax1.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated $\Delta DM_{FRB}$')
ax1.plot(dm_grid, deft_optimal_frb_len, color='crimson',linewidth=1,label=r'$\Delta DM_{FRB}$ sim ($n=76$)')
ax1.fill_between(dm_grid, deft_optimal_frb_len-deft_std_frb_len, deft_optimal_frb_len+deft_std_frb_len, alpha=.4,color='crimson')
ax1.set(ylabel='PDF')
ax1.yaxis.label.set_size(14)
ax1.legend(fontsize=14)
ax1.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)

ax2.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated $\Delta DM_{FRB}$')
ax2.plot(dm_grid, deft_optimal_1000, color='crimson',linewidth=1,label=r'$\Delta DM_{FRB}$ sim ($n=1000$)')
ax2.fill_between(dm_grid, deft_optimal_1000-deft_std_1000, deft_optimal_1000+deft_std_1000, alpha=.4,color='crimson')
ax2.set(ylabel='PDF')
ax2.legend(fontsize=14)
ax2.yaxis.label.set_size(14)
ax2.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)

ax3.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated $\Delta DM_{FRB}$')
ax3.plot(dm_grid, deft_optimal_5000, color='crimson',linewidth=1,label=r'$\Delta DM_{FRB}$ sim ($n=5000$)')
ax3.fill_between(dm_grid, deft_optimal_5000-deft_std_5000, deft_optimal_5000+deft_std_5000, alpha=.4,color='crimson')
ax3.set(ylabel='PDF')
ax3.legend(fontsize=14)
ax3.yaxis.label.set_size(14)
ax3.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)

ax4.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated $\Delta DM_{FRB}$')
ax4.plot(dm_grid, deft_optimal_10000, color='crimson',linewidth=1,label=r'$\Delta DM_{FRB}$ sim ($n=10000$)')
ax4.fill_between(dm_grid, deft_optimal_10000-deft_std_10000, deft_optimal_10000+deft_std_10000, alpha=.4,color='crimson')
ax4.legend(fontsize=14)
ax4.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax4.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
ax4.yaxis.label.set_size(14)
ax4.xaxis.label.set_size(14)
plt.tight_layout()
plt.savefig('sim_outputs/dm_sim_vary_n.png', dpi=300)
plt.show()

# matplotlib.rc('xtick', labelsize=10) 
# matplotlib.rc('ytick', labelsize=10) 
# plt.figure(figsize=(10,5))
# plt.hist(np.array(frbcat_df['dmdiff']),density=True,histtype='stepfilled',bins=40,alpha=0.5,label=r'Observed $\Delta$ DM$_{FRB}$',color='teal')
# plt.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_{FRB}$')
# plt.plot(dm_grid, deft_optimal_frb_len, color='crimson',linewidth=1,label=r'DM$_{FRB}$ sim ($n=76$)')
# plt.plot(grid_frb, dmdiff_frb_optimal*100, color='teal',linewidth=1,label=r'$\Delta$ DM$_{FRB}$ obs ($n=76$)')
# plt.fill_between(dm_grid, deft_optimal_frb_len-deft_std_frb_len, deft_optimal_frb_len+deft_std_frb_len, alpha=.4,color='crimson')
# plt.ylabel('PDF',fontsize=14)
# plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=14)
# plt.xlim([0,2000])
# plt.legend(fontsize=14)
# plt.tight_layout()
# plt.savefig('sim_outputs/dm_sim_with_obs.png', dpi=300)
# plt.show()
