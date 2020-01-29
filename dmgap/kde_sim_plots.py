import numpy as np
import pandas as pd
from pandas import DataFrame as df
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.arange(0,6000,0.01)

dm_frb_sim = np.load('sim_output_data/dm_frb_sim.npy')
dm_grid = np.load('sim_output_data/dm_grid.npy')

optimal_100 = np.load('sim_output_data/100_optimal.npy')
optimal_1000 = np.load('sim_output_data/1000_optimal.npy')
optimal_2000 = np.load('sim_output_data/2000_optimal.npy')
optimal_5000 = np.load('sim_output_data/5000_optimal.npy')
optimal_10000 = np.load('sim_output_data/10000_optimal.npy')

sampled_100 = np.load('sim_output_data/100_sampled.npy')
sampled_1000 = np.load('sim_output_data/1000_sampled.npy')
sampled_2000 = np.load('sim_output_data/2000_sampled.npy')
sampled_5000 = np.load('sim_output_data/5000_sampled.npy')
sampled_10000 = np.load('sim_output_data/10000_sampled.npy')

std_100 = np.load('sim_output_data/std_100.npy')
std_1000 = np.load('sim_output_data/std_1000.npy')
std_2000 = np.load('sim_output_data/std_2000.npy')
std_5000 = np.load('sim_output_data/std_5000.npy')
std_10000 = np.load('sim_output_data/std_10000.npy')

kde_optimal_100 = np.load('kde_outputs/kde_100_optimal.npy')
kde_optimal_1000 = np.load('kde_outputs/kde_1000_optimal.npy')
kde_optimal_2000 = np.load('kde_outputs/kde_2000_optimal.npy')
kde_optimal_5000 = np.load('kde_outputs/kde_5000_optimal.npy')
kde_optimal_10000 = np.load('kde_outputs/kde_10000_optimal.npy')

kde_bootstrapped_100 = np.load('kde_outputs/kde_100_bootstrapped.npy')
kde_bootstrapped_1000 = np.load('kde_outputs/kde_1000_bootstrapped.npy')
kde_bootstrapped_2000 = np.load('kde_outputs/kde_2000_bootstrapped.npy')
kde_bootstrapped_5000 = np.load('kde_outputs/kde_5000_bootstrapped.npy')
kde_bootstrapped_10000 = np.load('kde_outputs/kde_10000_bootstrapped.npy')

kde_std_100 = np.load('kde_outputs/kde_std_100.npy')
kde_std_1000 = np.load('kde_outputs/kde_std_1000.npy')
kde_std_2000 = np.load('kde_outputs/kde_std_2000.npy')
kde_std_5000 = np.load('kde_outputs/kde_std_5000.npy')
kde_std_10000 = np.load('kde_outputs/kde_std_10000.npy')

# Plots
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex='col', sharey='row', figsize=(8,12))
ax1.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax1.plot(dm_grid, kde_optimal_100, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=100$)')
ax1.fill_between(dm_grid, kde_optimal_100-kde_std_100, kde_optimal_100+kde_std_100, alpha=.3,color='teal')
ax1.plot(dm_grid, optimal_100, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=100$)')
ax1.fill_between(dm_grid, optimal_100-std_100, optimal_100+std_100, alpha=.4,color='crimson')
ax1.legend(fontsize=14)
ax1.set(ylabel='PDF')
ax1.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax1.xaxis.label.set_size(14)
ax1.yaxis.label.set_size(14)

ax2.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax2.plot(dm_grid, kde_optimal_1000, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=1000$)')
ax2.fill_between(dm_grid, kde_optimal_1000-kde_std_1000, kde_optimal_1000+kde_std_1000, alpha=.3,color='teal')
ax2.plot(dm_grid, optimal_1000, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=1000$)')
ax2.fill_between(dm_grid, optimal_1000-std_1000, optimal_1000+std_1000, alpha=.4,color='crimson')
ax2.set(ylabel='PDF')
ax2.legend(fontsize=14)
ax2.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax2.xaxis.label.set_size(14)
ax2.yaxis.label.set_size(14)

ax3.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax3.plot(dm_grid, kde_optimal_5000, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=5000$)')
ax3.fill_between(dm_grid, kde_optimal_5000-kde_std_5000, kde_optimal_5000+kde_std_5000, alpha=.3,color='teal')
ax3.plot(dm_grid, optimal_5000, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=5000$)')
ax3.fill_between(dm_grid, optimal_5000-std_5000, optimal_5000+std_5000, alpha=.4,color='crimson')
ax3.set(ylabel='PDF')
ax3.legend(fontsize=14)
ax3.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax3.xaxis.label.set_size(14)
ax3.yaxis.label.set_size(14)

ax4.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax4.plot(dm_grid, kde_optimal_10000, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=10000$)')
ax4.fill_between(dm_grid, kde_optimal_10000-kde_std_10000, kde_optimal_10000+kde_std_10000, alpha=.3,color='teal')
ax4.plot(dm_grid, optimal_10000, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=10000$)')
ax4.fill_between(dm_grid, optimal_10000-std_10000, optimal_10000+std_10000, alpha=.4,color='crimson')
ax4.legend(fontsize=14)
ax4.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax4.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
ax4.xaxis.label.set_size(14)
ax4.yaxis.label.set_size(14)

plt.tight_layout()
plt.savefig('kde_outputs/kde_sim_vary_n.png', dpi=300)
plt.show()
