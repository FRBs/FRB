import numpy as np
import pandas as pd
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw
from collections import OrderedDict
from percentile_defs import find_percentile

# Data upload 
dm_frb_sim = np.load('sim_output_data/dm_frb_sim.npy')
dm_grid = np.load('sim_output_data/dm_grid.npy')

# Constraints to region of interest for max concavity
dm_min = -50
dm_max = 300
dm_min_idx = (np.abs(dm_grid - dm_min)).argmin()
dm_max_idx = (np.abs(dm_grid - dm_max)).argmin()

optimal_100 = np.load('sim_output_data/optimal_100.npy')
optimal_100 = optimal_100[dm_min_idx:dm_max_idx]
optimal_1000 = np.load('sim_output_data/optimal_1000.npy')
optimal_1000 = optimal_1000[dm_min_idx:dm_max_idx]
optimal_5000 = np.load('sim_output_data/optimal_5000.npy')
optimal_5000 = optimal_5000[dm_min_idx:dm_max_idx]
optimal_10000 = np.load('sim_output_data/optimal_10000.npy')
optimal_10000 = optimal_10000[dm_min_idx:dm_max_idx]

sampled_100 = np.load('sim_output_data/sampled_100.npy')
sampled_100 = sampled_100[dm_min_idx:dm_max_idx]
sampled_1000 = np.load('sim_output_data/sampled_1000.npy')
sampled_1000 = sampled_1000[dm_min_idx:dm_max_idx]
sampled_5000 = np.load('sim_output_data/sampled_5000.npy')
sampled_5000 = sampled_5000[dm_min_idx:dm_max_idx]
sampled_10000 = np.load('sim_output_data/sampled_10000.npy')
sampled_10000 = sampled_10000[dm_min_idx:dm_max_idx]

dm_grid_pp = dm_grid[dm_min_idx:dm_max_idx]

# 100
frb_sim_pp_100 = np.diff(np.diff(optimal_100))
frb_sim_gap_100 = dm_grid_pp[frb_sim_pp_100.argmax()]

frb_sim_gaps_100 = []
for i in range(len(sampled_100[1])):
    frb_sim_pp_ensemble = np.diff(np.diff(sampled_100[:,i]))
    frb_sim_gaps_ = dm_grid_pp[frb_sim_pp_ensemble.argmax()]
    frb_sim_gaps_100 = np.append(frb_sim_gaps_100,frb_sim_gaps_)
frb_sim_gaps_100 = np.append(frb_sim_gap_100,frb_sim_gaps_100)
frb_sim_gaps_100 = list(set(frb_sim_gaps_100))

print(len(frb_sim_gaps_100))

alpha_param = 3
num_ensembles = 100
b_box_100 = [min(frb_sim_gaps_100)-100,max(frb_sim_gaps_100)+200]
density_gap_100 = sw.DensityEstimator(frb_sim_gaps_100, alpha=alpha_param, bounding_box=b_box_100, num_posterior_samples=num_ensembles)
gap_100_pdf = density_gap_100.evaluate(dm_grid_pp)
gap_std_100 = np.std(frb_sim_gaps_100)

plt.plot(dm_grid_pp,gap_100_pdf)
plt.hist(frb_sim_gaps_100,bins=20,density=True,alpha=.2)
plt.show()

# # 1000
# frb_sim_pp_1000 = np.diff(np.diff(optimal_1000))
# frb_sim_gap_1000 = dm_grid_pp[frb_sim_pp_1000.argmax()]

# frb_sim_gaps_1000 = []
# for i in range(len(sampled_1000[1])):
#     frb_sim_pp_ensemble = np.diff(np.diff(sampled_1000[:,i]))
#     frb_sim_gaps_ = dm_grid_pp[frb_sim_pp_ensemble.argmax()]
#     frb_sim_gaps_1000 = np.append(frb_sim_gaps_1000,frb_sim_gaps_)
# frb_sim_gaps_1000 = np.append(frb_sim_gap_1000,frb_sim_gaps_1000)

# b_box_1000 = [min(frb_sim_gaps_1000)-50,max(frb_sim_gaps_1000)+50]
# density_gap_1000 = sw.DensityEstimator(frb_sim_gaps_1000, alpha=alpha_param, bounding_box=b_box_1000, num_posterior_samples=num_ensembles)
# gap_1000_pdf = density_gap_1000.evaluate(dm_grid_pp)
# gap_std_1000 = np.std(frb_sim_gaps_1000)

# # 5000
# frb_sim_pp_5000 = np.diff(np.diff(optimal_5000))
# frb_sim_gap_5000 = dm_grid_pp[frb_sim_pp_5000.argmax()]

# frb_sim_gaps_5000 = []
# for i in range(len(sampled_5000[1])):
#     frb_sim_pp_ensemble = np.diff(np.diff(sampled_5000[:,i]))
#     frb_sim_gaps_ = dm_grid_pp[frb_sim_pp_ensemble.argmax()]
#     frb_sim_gaps_5000 = np.append(frb_sim_gaps_5000,frb_sim_gaps_)
# frb_sim_gaps_5000 = np.append(frb_sim_gap_5000,frb_sim_gaps_5000)

# b_box_5000 = [min(frb_sim_gaps_5000)-50,max(frb_sim_gaps_5000)+50]
# density_gap_5000 = sw.DensityEstimator(frb_sim_gaps_5000, alpha=alpha_param, bounding_box=b_box_5000, num_posterior_samples=num_ensembles)
# gap_5000_pdf = density_gap_5000.evaluate(dm_grid_pp)
# gap_std_5000 = np.std(frb_sim_gaps_5000)

# # 10000
# frb_sim_pp_10000 = np.diff(np.diff(optimal_10000))
# frb_sim_gap_10000 = dm_grid_pp[frb_sim_pp_10000.argmax()]

# frb_sim_gaps_10000 = []
# for i in range(len(sampled_10000[1])):
#     frb_sim_pp_ensemble = np.diff(np.diff(sampled_10000[:,i]))
#     frb_sim_gaps_ = dm_grid_pp[frb_sim_pp_ensemble.argmax()]
#     frb_sim_gaps_10000 = np.append(frb_sim_gaps_10000,frb_sim_gaps_)
# frb_sim_gaps_10000 = np.append(frb_sim_gap_10000,frb_sim_gaps_10000)

# b_box_10000 = [min(frb_sim_gaps_10000)-50,max(frb_sim_gaps_10000)+50]
# density_gap_10000 = sw.DensityEstimator(frb_sim_gaps_10000, alpha=alpha_param, bounding_box=b_box_10000, num_posterior_samples=num_ensembles)
# gap_10000_pdf = density_gap_10000.evaluate(dm_grid_pp)
# gap_std_10000 = np.std(frb_sim_gaps_10000)

# # Upper and lower gap constraints
# p = 95. #percentile
# p_lower = (100-p)/2
# p_upper = p+(100-p)/2

# gap_100_lower = find_percentile(dm_grid_pp, gap_100_pdf, p_lower) 
# gap_100_upper = find_percentile(dm_grid_pp, gap_100_pdf, p_upper)
# gap_1000_lower = find_percentile(dm_grid_pp, gap_1000_pdf, p_lower) 
# gap_1000_upper = find_percentile(dm_grid_pp, gap_1000_pdf, p_upper)
# gap_5000_lower = find_percentile(dm_grid_pp, gap_5000_pdf, p_lower) 
# gap_5000_upper = find_percentile(dm_grid_pp, gap_5000_pdf, p_upper)
# gap_10000_lower = find_percentile(dm_grid_pp, gap_10000_pdf, p_lower)
# gap_10000_upper = find_percentile(dm_grid_pp, gap_10000_pdf, p_upper)

# print('100 Optimal:',frb_sim_gap_100)
# print('100 lower bound:', gap_100_lower) 
# print('100 upper bound:', gap_100_upper) 
# print('Std is:', gap_std_100)
# print('_________')
# print('1000 Optimal:',frb_sim_gap_1000)
# print('1000 lower bound:', gap_1000_lower) 
# print('1000 upper bound:', gap_1000_upper) 
# print('Std is:', gap_std_1000)
# print('_________')
# print('5000 Optimal:',frb_sim_gap_5000)
# print('5000 lower bound:', gap_5000_lower) 
# print('5000 upper bound:', gap_5000_upper) 
# print('Std is:', gap_std_5000)
# print('_________')
# print('10000 Optimal:',frb_sim_gap_10000)
# print('10000 lower bound:', gap_10000_lower) 
# print('10000 upper bound:', gap_10000_upper) 
# print('Std is:', gap_std_10000)
# print('_________')

# matplotlib.rc('xtick', labelsize=12) 
# matplotlib.rc('ytick', labelsize=12) 
# plt.figure(figsize=(12,8))
# plt.hist(frb_sim_gaps_100, density=True, bins=100, histtype='stepfilled', alpha=.2, color='blue')
# plt.plot(dm_grid_pp,gap_100_pdf,color='blue',label=r'$n=100$')
# plt.hist(frb_sim_gaps_1000, density=True, bins=80, histtype='stepfilled', alpha=.2, color='teal')
# plt.plot(dm_grid_pp,gap_1000_pdf,color='teal',label=r'$n=1000$')
# plt.hist(frb_sim_gaps_5000, density=True, bins=80, histtype='stepfilled', alpha=.2, color='green')
# plt.plot(dm_grid_pp,gap_5000_pdf,color='green',label=r'$n=5000$')
# plt.hist(frb_sim_gaps_10000, density=True, bins=40, histtype='stepfilled', alpha=.3, color='crimson')
# plt.plot(dm_grid_pp,gap_10000_pdf,color='crimson',label=r'$n=10000$')
# plt.axvline(x=90,color='k',linestyle='--',linewidth=1,label=r'$\Delta$DM$_\mathrm{FRB}$=90')
# plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
# plt.ylabel('PDF',fontsize=22)
# plt.xlim(-50,300)
# plt.ylim(0,0.075)
# plt.legend(fontsize=22,loc=1)
# plt.tight_layout()
# plt.savefig('sim_output_figs/DM_FRB_sim.png', dpi=300)
# plt.show()