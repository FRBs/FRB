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
import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 

# Data upload 
dmdiff_psr_optimal = np.load('obs_outputs/dmdiff_psr_optimal.npy')
grid_psr = np.arange(-500,500,0.01)
dm_grid = np.load('sim_outputs/dm_grid.npy')

deft_optimal_frb_len = np.load('kde_and_deft_data/deft_frb_len_optimal.npy')
deft_optimal_frb_len = deft_optimal_frb_len/np.sum(deft_optimal_frb_len)
deft_optimal_100 = np.load('kde_and_deft_data/deft_100_optimal.npy')
deft_optimal_100 = deft_optimal_100/np.sum(deft_optimal_100)
deft_optimal_1000 = np.load('kde_and_deft_data/deft_1000_optimal.npy')
deft_optimal_1000 = deft_optimal_1000/np.sum(deft_optimal_1000)
deft_optimal_2000 = np.load('kde_and_deft_data/deft_2000_optimal.npy')
deft_optimal_2000 = deft_optimal_2000/np.sum(deft_optimal_2000)
deft_optimal_4000 = np.load('kde_and_deft_data/deft_4000_optimal.npy')
deft_optimal_4000 = deft_optimal_4000/np.sum(deft_optimal_4000)

deft_sampled_frb_len = np.load('kde_and_deft_data/deft_frb_len_sampled.npy')
deft_sampled_frb_len = deft_sampled_frb_len/np.sum(deft_sampled_frb_len)
deft_sampled_100 = np.load('kde_and_deft_data/deft_100_sampled.npy')
deft_sampled_100 = deft_sampled_100/np.sum(deft_sampled_100)
deft_sampled_1000 = np.load('kde_and_deft_data/deft_1000_sampled.npy')
deft_sampled_1000 = deft_sampled_1000/np.sum(deft_sampled_1000)
deft_sampled_2000 = np.load('kde_and_deft_data/deft_2000_sampled.npy')
deft_sampled_2000 = deft_sampled_2000/np.sum(deft_sampled_2000)
deft_sampled_4000 = np.load('kde_and_deft_data/deft_4000_sampled.npy')
deft_sampled_4000 = deft_sampled_4000/np.sum(deft_sampled_4000)

deft_std_frb_len = np.load('kde_and_deft_data/deft_std_frb_len.npy')
deft_std_frb_len = deft_std_frb_len/np.sum(deft_std_frb_len)
deft_std_100 = np.load('kde_and_deft_data/deft_std_100.npy')
deft_std_100 = deft_std_100/np.sum(deft_std_100)
deft_std_1000 = np.load('kde_and_deft_data/deft_std_1000.npy')
deft_std_1000 = deft_std_1000/np.sum(deft_std_1000)
deft_std_2000 = np.load('kde_and_deft_data/deft_std_2000.npy')
deft_std_2000 = deft_std_2000/np.sum(deft_std_2000)
deft_std_4000 = np.load('kde_and_deft_data/deft_std_4000.npy')
deft_std_4000 = deft_std_4000/np.sum(deft_std_4000)

ci_interval = [.95, .96, .97, .98, .99] #confidence
ci_cut_off = np.sqrt(ci_interval) #for joint probability: probability that <ci_interval% of events will occur simultaneously

frb_sim_frb_len = []
for i in range(len(ci_cut_off)):
    frb_sim_frb_len_ = find_quantile(grid=dm_grid,distribution=deft_optimal_frb_len,quantile=1-ci_cut_off[i])
    frb_sim_frb_len = np.append(frb_sim_frb_len,frb_sim_frb_len_)

frb_sim_1000 = []
for i in range(len(ci_cut_off)):
    frb_sim_1000_ = find_quantile(grid=dm_grid,distribution=deft_optimal_1000,quantile=1-ci_cut_off[i])
    frb_sim_1000 = np.append(frb_sim_1000,frb_sim_1000_)

frb_sim_2000 = []
for i in range(len(ci_cut_off)):
    frb_sim_2000_ = find_quantile(grid=dm_grid,distribution=deft_optimal_2000,quantile=1-ci_cut_off[i])
    frb_sim_2000 = np.append(frb_sim_2000,frb_sim_2000_)

frb_sim_4000 = []
for i in range(len(ci_cut_off)):
    frb_sim_4000_ = find_quantile(grid=dm_grid,distribution=deft_optimal_4000,quantile=1-ci_cut_off[i])
    frb_sim_4000 = np.append(frb_sim_4000,frb_sim_4000_)

print(frb_sim_frb_len)
print(frb_sim_1000)
print(frb_sim_2000)
print(frb_sim_4000)

psr_optimal_dm = find_quantile(grid=grid_psr,distribution=dmdiff_psr_optimal,quantile=ci_cut_off)

# Save results in dataframe
df_dm_ci = pd.DataFrame(columns=['CI', 'dm_halo_psr', 'dm_halo_frb_num_obs', 'dm_halo_frb_1000', 'dm_halo_frb_2000', 'dm_halo_frb_4000'])
print(df_dm_ci)