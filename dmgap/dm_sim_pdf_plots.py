import numpy as np
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 
plt.figure(figsize=(12,8))

# Data upload 
dm_frb_sim = np.load('sim_output_data/dm_frb_sim.npy')
dm_grid = np.load('sim_output_data/dm_grid.npy')

optimal_10000 = np.load('sim_output_data/10000_optimal.npy')
sampled_10000 = np.load('sim_output_data/10000_sampled.npy')

# DEFT with 10000 draws
plt.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
plt.plot(dm_grid,optimal_10000,color='crimson',label=r'DEFT DM$_\mathrm{FRB}$ (n=10000)')
plt.plot(dm_grid,sampled_10000,color='crimson',linewidth=0.5,alpha=.1)
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
plt.ylabel('PDF',fontsize=22)
plt.xlim(0,2000)
plt.ylim(0,0.0015)
plt.legend(fontsize=22)
plt.tight_layout()
plt.savefig('sim_output_figs/DM_10000_samples.png', dpi=300)
plt.show()