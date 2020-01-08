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

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.arange(0,6000,0.01)

dm_frb_sim = np.load('sim_outputs/dm_frb_sim.npy')
dm_grid = np.load('sim_outputs/dm_grid.npy')

deft_optimal_100 = np.load('kde_and_deft_data/deft_100_optimal.npy')
deft_optimal_1000 = np.load('kde_and_deft_data/deft_1000_optimal.npy')
deft_optimal_2000 = np.load('kde_and_deft_data/deft_2000_optimal.npy')
deft_optimal_4000 = np.load('kde_and_deft_data/deft_4000_optimal.npy')

deft_sampled_100 = np.load('kde_and_deft_data/deft_100_sampled.npy')
deft_sampled_1000 = np.load('kde_and_deft_data/deft_1000_sampled.npy')
deft_sampled_2000 = np.load('kde_and_deft_data/deft_2000_sampled.npy')
deft_sampled_4000 = np.load('kde_and_deft_data/deft_4000_sampled.npy')

deft_std_100 = np.load('kde_and_deft_data/deft_std_100.npy')
deft_std_1000 = np.load('kde_and_deft_data/deft_std_1000.npy')
deft_std_2000 = np.load('kde_and_deft_data/deft_std_2000.npy')
deft_std_4000 = np.load('kde_and_deft_data/deft_std_4000.npy')

kde_optimal_100 = np.load('kde_and_deft_data/kde_100_optimal.npy')
kde_optimal_1000 = np.load('kde_and_deft_data/kde_1000_optimal.npy')
kde_optimal_2000 = np.load('kde_and_deft_data/kde_2000_optimal.npy')
kde_optimal_4000 = np.load('kde_and_deft_data/kde_4000_optimal.npy')

kde_bootstrapped_100 = np.load('kde_and_deft_data/kde_100_bootstrapped.npy')
kde_bootstrapped_1000 = np.load('kde_and_deft_data/kde_1000_bootstrapped.npy')
kde_bootstrapped_2000 = np.load('kde_and_deft_data/kde_2000_bootstrapped.npy')
kde_bootstrapped_4000 = np.load('kde_and_deft_data/kde_4000_bootstrapped.npy')

kde_std_100 = np.load('kde_and_deft_data/kde_std_100.npy')
kde_std_1000 = np.load('kde_and_deft_data/kde_std_1000.npy')
kde_std_2000 = np.load('kde_and_deft_data/kde_std_2000.npy')
kde_std_4000 = np.load('kde_and_deft_data/kde_std_4000.npy')

# Plots
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex='col', sharey='row', figsize=(10,12))
ax1.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax1.plot(dm_grid, kde_optimal_100, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=100$)')
ax1.fill_between(dm_grid, kde_optimal_100-kde_std_100, kde_optimal_100+kde_std_100, alpha=.3,color='teal')
ax1.plot(dm_grid, deft_optimal_100, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=100$)')
ax1.fill_between(dm_grid, deft_optimal_100-deft_std_100, deft_optimal_100+deft_std_100, alpha=.4,color='crimson')
# ax1.set_yticklabels([])
ax1.legend(fontsize=16)
ax1.set(ylabel='PDF')
ax1.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax1.xaxis.label.set_size(16)
ax1.yaxis.label.set_size(16)

ax2.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax2.plot(dm_grid, kde_optimal_1000, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=1000$)')
ax2.fill_between(dm_grid, kde_optimal_1000-kde_std_1000, kde_optimal_1000+kde_std_1000, alpha=.3,color='teal')
ax2.plot(dm_grid, deft_optimal_1000, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=1000$)')
ax2.fill_between(dm_grid, deft_optimal_1000-deft_std_1000, deft_optimal_1000+deft_std_1000, alpha=.4,color='crimson')
# ax2.set_yticklabels([])
ax2.set(ylabel='PDF')
ax2.legend(fontsize=16)
ax2.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax2.xaxis.label.set_size(16)
ax2.yaxis.label.set_size(16)

ax3.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax3.plot(dm_grid, kde_optimal_2000, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=2000$)')
ax3.fill_between(dm_grid, kde_optimal_2000-kde_std_2000, kde_optimal_2000+kde_std_2000, alpha=.3,color='teal')
ax3.plot(dm_grid, deft_optimal_2000, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=2000$)')
ax3.fill_between(dm_grid, deft_optimal_2000-deft_std_2000, deft_optimal_2000+deft_std_2000, alpha=.4,color='crimson')
ax3.set(ylabel='PDF')
ax3.legend(fontsize=14)
ax3.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax3.xaxis.label.set_size(16)
ax3.yaxis.label.set_size(16)

ax4.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax4.plot(dm_grid, kde_optimal_4000, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=4000$)')
ax4.fill_between(dm_grid, kde_optimal_4000-kde_std_4000, kde_optimal_4000+kde_std_4000, alpha=.3,color='teal')
ax4.plot(dm_grid, deft_optimal_4000, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=4000$)')
ax4.fill_between(dm_grid, deft_optimal_4000-deft_std_4000, deft_optimal_4000+deft_std_4000, alpha=.4,color='crimson')
# ax4.set_yticklabels([])
ax4.legend(fontsize=16)
ax4.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax4.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
ax4.xaxis.label.set_size(16)
ax4.yaxis.label.set_size(16)

plt.tight_layout()
plt.savefig('sim_outputs/dm_sim_vary_n.png', dpi=300)
plt.show()
