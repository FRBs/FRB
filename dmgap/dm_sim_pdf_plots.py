import numpy as np
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt

# Data upload 
dm_frb_sim = np.load('sim_output_data/dm_frb_sim.npy')
dm_grid = np.load('sim_output_data/dm_grid.npy')

optimal_100 = np.load('sim_output_data/optimal_100.npy')
sampled_100 = np.load('sim_output_data/sampled_100.npy')
optimal_1000 = np.load('sim_output_data/optimal_1000.npy')
sampled_1000 = np.load('sim_output_data/sampled_1000.npy')
optimal_10000 = np.load('sim_output_data/optimal_10000.npy')
sampled_10000 = np.load('sim_output_data/sampled_10000.npy')
std_100 = np.load('sim_output_data/std_100.npy')
std_1000 = np.load('sim_output_data/std_1000.npy')
std_10000 = np.load('sim_output_data/std_10000.npy')

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12)
plt.figure(figsize=(12,8))
plt.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='k', label=r'$\Delta$DM$_\mathrm{FRB,sim}$')
plt.plot(dm_grid,optimal_100,color='slateblue',label=r'$\Delta$DM$_\mathrm{FRB}$ ($n=100$)')
plt.fill_between(dm_grid, optimal_100-std_100, optimal_100+std_100, alpha=.4,color='slateblue')
plt.plot(dm_grid,optimal_1000,color='teal',label=r'$\Delta$DM$_\mathrm{FRB}$ ($n=1000$)')
plt.fill_between(dm_grid, optimal_1000-std_1000, optimal_1000+std_1000, alpha=.4,color='teal')
plt.plot(dm_grid,optimal_10000,color='crimson',label=r'$\Delta$DM$_\mathrm{FRB}$ ($n=10000$)')
plt.fill_between(dm_grid, optimal_10000-std_10000, optimal_10000+std_10000, alpha=.4,color='crimson')
plt.axvline(x=min(dm_frb_sim),linewidth=1,color='k',linestyle='--', label=r'DM$_\mathrm{gap,sim}=90$pc cm$^{-3}$')
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
plt.ylabel('PDF',fontsize=22)
plt.xlim(-200,2000)
plt.ylim(0,0.0017)
plt.legend(fontsize=22,loc=1)
plt.tight_layout()
plt.savefig('sim_output_figs/DM_n_samples.png', dpi=300)
plt.show()