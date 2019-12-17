import sys
sys.path.insert(1, '/eplatts_UCSC_Server/dm_gap/ne2001-master/src')
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw
from percentile_defs import find_quantile

matplotlib.rc('xtick', labelsize=22) 
matplotlib.rc('ytick', labelsize=22) 

# Data upload 
dm_frb_sim = np.load('sim_outputs/dm_frb_sim.npy')
dm_grid = np.load('sim_outputs/dm_grid.npy')

deft_optimal_frb_len = np.load('kde_and_deft_data/deft_frb_len_optimal.npy')
# deft_optimal_frb_len = deft_optimal_frb_len/np.sum(deft_optimal_frb_len)
deft_optimal_100 = np.load('kde_and_deft_data/deft_100_optimal.npy')
# print(sum(deft_optimal_100))
# deft_optimal_100 = deft_optimal_100/np.sum(deft_optimal_100)
deft_optimal_1000 = np.load('kde_and_deft_data/deft_1000_optimal.npy')
# deft_optimal_1000 = deft_optimal_1000/np.sum(deft_optimal_1000)
deft_optimal_2000 = np.load('kde_and_deft_data/deft_2000_optimal.npy')
# deft_optimal_2000 = deft_optimal_2000/np.sum(deft_optimal_2000)
deft_optimal_4000 = np.load('kde_and_deft_data/deft_4000_optimal.npy')
# deft_optimal_4000 = deft_optimal_4000/np.sum(deft_optimal_4000)
deft_optimal_10000 = np.load('kde_and_deft_data/deft_10000_optimal.npy')

deft_sampled_frb_len = np.load('kde_and_deft_data/deft_frb_len_sampled.npy')
# deft_sampled_frb_len = deft_sampled_frb_len/np.sum(deft_sampled_frb_len)
deft_sampled_100 = np.load('kde_and_deft_data/deft_100_sampled.npy')
# deft_sampled_100 = deft_sampled_100/np.sum(deft_sampled_100)
deft_sampled_1000 = np.load('kde_and_deft_data/deft_1000_sampled.npy')
# deft_sampled_1000 = deft_sampled_1000/np.sum(deft_sampled_1000)
deft_sampled_2000 = np.load('kde_and_deft_data/deft_2000_sampled.npy')
# deft_sampled_2000 = deft_sampled_2000/np.sum(deft_sampled_2000)
deft_sampled_4000 = np.load('kde_and_deft_data/deft_4000_sampled.npy')
# deft_sampled_4000 = deft_sampled_4000/np.sum(deft_sampled_4000)
deft_sampled_10000 = np.load('kde_and_deft_data/deft_10000_sampled.npy')

ci_interval = .95 #confidence
ci_cut_off = np.sqrt(ci_interval) #for joint probability: probability that <ci_interval% of events will occur simultaneously

frb_sim_1000 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_1000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_1000[0:,i],quantile=1-ci_cut_off)
    frb_sim_1000 = np.append(frb_sim_1000,frb_sim_1000_)
    
frb_sim_2000 = []
for i in range(len(deft_sampled_2000[1])):
    frb_sim_2000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_2000[0:,i],quantile=1-ci_cut_off)
    frb_sim_2000 = np.append(frb_sim_2000,frb_sim_2000_)

frb_sim_4000 = []
for i in range(len(deft_sampled_4000[1])):
    frb_sim_4000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_4000[0:,i],quantile=1-ci_cut_off)
    frb_sim_4000 = np.append(frb_sim_4000,frb_sim_4000_)

frb_sim_10000 = []
for i in range(len(deft_sampled_10000[1])):
    frb_sim_10000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_10000[0:,i],quantile=1-ci_cut_off)
    frb_sim_10000 = np.append(frb_sim_10000,frb_sim_10000_)

print(find_quantile(grid=dm_grid,distribution=deft_optimal_1000,quantile=1-ci_cut_off))
print(find_quantile(grid=dm_grid,distribution=deft_optimal_2000,quantile=1-ci_cut_off))
print(find_quantile(grid=dm_grid,distribution=deft_optimal_4000,quantile=1-ci_cut_off))
print(find_quantile(grid=dm_grid,distribution=deft_optimal_10000,quantile=1-ci_cut_off))

# plt.hist(frb_sim_1000,bins=50, alpha=0.5)
# plt.hist(frb_sim_2000,bins=50, alpha=0.5)
# plt.hist(frb_sim_4000,bins=50, alpha=0.5)
# plt.hist(frb_sim_10000,bins=50,alpha=0.5)
# plt.axvline(x=50,color='k')

# plt.show()

# psr_optimal_dm = find_quantile(grid=grid_psr,distribution=dmdiff_psr_optimal,quantile=ci_cut_off)

# # Save results in dataframe
# df_dm_ci = pd.DataFrame(columns=['CI', 'dm_halo_psr', 'dm_halo_frb_num_obs', 'dm_halo_frb_1000', 'dm_halo_frb_2000', 'dm_halo_frb_4000'])
# print(df_dm_ci)