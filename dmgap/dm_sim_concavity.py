import numpy as np
import pandas as pd
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw

# Data upload 
dm_frb_sim = np.load('sim_output_data/dm_frb_sim.npy')
dm_grid = np.load('sim_output_data/dm_grid.npy')

optimal_100 = np.load('sim_output_data/frb_100_optimal.npy')
optimal_1000 = np.load('sim_output_data/frb_1000_optimal.npy')
optimal_10000 = np.load('sim_output_data/frb_10000_optimal.npy')

sampled_100 = np.load('sim_output_data/frb_100_sampled.npy')
sampled_1000 = np.load('sim_output_data/frb_1000_sampled.npy')
sampled_10000 = np.load('sim_output_data/frb_10000_sampled.npy')

# # CDFs
# cdf_frb = np.cumsum(deltaDM_frb_optimal)
# cdf_frb_pp = np.diff(np.diff(cdf_frb)) #second derivative
# grid_cdf_frb = grid_frb[2:]
# y_cdf_gap = max(cdf_frb_pp)
# print(y_cdf_gap)

# max_cdf_y = y_cdf_gap  # Find the maximum y value
# max_cdf_x = grid_cdf_frb[cdf_frb_pp.argmax()]
# print(max_cdf_x)

dm_grid_pp = dm_grid[35000:90000]

# 100
frb_sim_pp_100 = np.diff(np.diff(optimal_100))[35000:80000]
frb_sim_gap_100 = dm_grid_pp[frb_sim_pp_100.argmax()]

frb_sim_gaps_100 = []
for i in range(len(sampled_100[1])):
    frb_sim_pp_ensemble = np.diff(np.diff(sampled_100[34998:80000,i]))
    frb_sim_gaps_ = dm_grid_pp[frb_sim_pp_ensemble.argmax()]
    frb_sim_gaps_100 = np.append(frb_sim_gaps_100,frb_sim_gaps_)
frb_sim_gaps_100 = np.append(frb_sim_gap_100,frb_sim_gaps_100)

alpha_param = 3
num_ensembles = 100
b_box_100 = [min(frb_sim_gaps_100)-100,max(frb_sim_gaps_100)+200]
density_gap_100 = sw.DensityEstimator(frb_sim_gaps_100, alpha=alpha_param, bounding_box=b_box_100, num_posterior_samples=num_ensembles)
gap_100_pdf = density_gap_100.evaluate(dm_grid_pp)

# 1000
frb_sim_pp_1000 = np.diff(np.diff(optimal_1000))[35000:80000]
frb_sim_gap_1000 = dm_grid_pp[frb_sim_pp_1000.argmax()]

frb_sim_gaps_1000 = []
for i in range(len(sampled_1000[1])):
    frb_sim_pp_ensemble = np.diff(np.diff(sampled_1000[34998:80000,i]))
    frb_sim_gaps_ = dm_grid_pp[frb_sim_pp_ensemble.argmax()]
    frb_sim_gaps_1000 = np.append(frb_sim_gaps_1000,frb_sim_gaps_)
frb_sim_gaps_1000 = np.append(frb_sim_gap_1000,frb_sim_gaps_1000)

b_box_1000 = [min(frb_sim_gaps_1000)-50,max(frb_sim_gaps_1000)+50]
density_gap_1000 = sw.DensityEstimator(frb_sim_gaps_1000, alpha=alpha_param, bounding_box=b_box_1000, num_posterior_samples=num_ensembles)
gap_1000_pdf = density_gap_1000.evaluate(dm_grid_pp)

# 10000
frb_sim_pp_10000 = np.diff(np.diff(optimal_10000))[35000:80000]
frb_sim_gap_10000 = dm_grid_pp[frb_sim_pp_10000.argmax()]

frb_sim_gaps_10000 = []
for i in range(len(sampled_10000[1])):
    frb_sim_pp_ensemble = np.diff(np.diff(sampled_10000[34998:80000,i]))
    frb_sim_gaps_ = dm_grid_pp[frb_sim_pp_ensemble.argmax()]
    frb_sim_gaps_10000 = np.append(frb_sim_gaps_10000,frb_sim_gaps_)
frb_sim_gaps_10000 = np.append(frb_sim_gap_10000,frb_sim_gaps_10000)

b_box_10000 = [min(frb_sim_gaps_10000)-50,max(frb_sim_gaps_10000)+50]
density_gap_10000 = sw.DensityEstimator(frb_sim_gaps_10000, alpha=alpha_param, bounding_box=b_box_10000, num_posterior_samples=num_ensembles)
gap_10000_pdf = density_gap_10000.evaluate(dm_grid_pp)


matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
plt.hist(frb_sim_gaps_100, density=True, bins=60, histtype='stepfilled', alpha=.4, color='crimson')
plt.plot(dm_grid_pp,gap_100_pdf,color='crimson',label=r'DM$_\mathrm{gap}$ (n=100)')
plt.hist(frb_sim_gaps_1000, density=True, bins=60, histtype='stepfilled', alpha=.4, color='green')
plt.plot(dm_grid_pp,gap_1000_pdf,color='green',label=r'DM$_\mathrm{gap}$ (n=1000)')
plt.hist(frb_sim_gaps_10000, density=True, bins=40, histtype='stepfilled', alpha=.2, color='blue')
plt.plot(dm_grid_pp,gap_10000_pdf,color='blue',label=r'DM$_\mathrm{gap}$ (n=10000)')
plt.axvline(x=90,linewidth=1,color='k',linestyle='--', label='DM$_\mathrm{gap,sim}=90$pc cm$^{-3}$')
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=14)
plt.ylabel('PDF',fontsize=12)
plt.xlim(0,200)
plt.legend(fontsize=12,loc=1)
plt.tight_layout()
plt.savefig('sim_output_figs/DM_gap.png', dpi=300)
plt.show()