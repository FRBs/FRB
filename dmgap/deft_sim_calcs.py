# PACKAGE IMPORT
import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import random
import suftware as sw
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
psrcat_df = pd.read_csv('transient_data/psrcat_df.csv')

dm_frb_sim = np.load('sim_outputs/dm_frb_sim.npy')
dm_grid = np.load('sim_outputs/dm_grid.npy')

alpha_param = 3

deft_sample_frb_len = np.asarray(random.sample(list(dm_frb_sim),len(frbcat_df['dmdiff'])))
deft_density_frb_len = sw.DensityEstimator(deft_sample_frb_len, alpha=alpha_param, bounding_box=[0,2500])
deft_sample_100 = np.asarray(random.sample(list(dm_frb_sim),100))
deft_density_100 = sw.DensityEstimator(deft_sample_100, alpha=alpha_param, bounding_box=[0,2500])
deft_sample_1000 = np.asarray(random.sample(list(dm_frb_sim),1000))
deft_density_1000 = sw.DensityEstimator(deft_sample_1000, alpha=alpha_param, bounding_box=[0,2500])
deft_sample_2000 = np.asarray(random.sample(list(dm_frb_sim),2000))
deft_density_2000 = sw.DensityEstimator(deft_sample_2000, alpha=alpha_param, bounding_box=[0,2500])
deft_sample_4000 = np.asarray(random.sample(list(dm_frb_sim),4000))
deft_density_4000 = sw.DensityEstimator(deft_sample_4000, alpha=alpha_param, bounding_box=[0,2500])
deft_sample_10000 = np.asarray(random.sample(list(dm_frb_sim),10000))
deft_density_10000 = sw.DensityEstimator(deft_sample_10000, alpha=alpha_param, bounding_box=[0,2500])

# Evaluate optimal density
deft_optimal_frb_len = deft_density_frb_len.evaluate(dm_grid)
deft_optimal_100 = deft_density_100.evaluate(dm_grid)
deft_optimal_1000 = deft_density_1000.evaluate(dm_grid)
deft_optimal_2000 = deft_density_2000.evaluate(dm_grid)
deft_optimal_4000 = deft_density_4000.evaluate(dm_grid)
deft_optimal_10000 = deft_density_10000.evaluate(dm_grid)
print('Density calculations done.')

# Evaluate sampled densities
deft_sampled_frb_len = deft_density_frb_len.evaluate_samples(dm_grid)
deft_sampled_100 = deft_density_100.evaluate_samples(dm_grid)
deft_sampled_1000 = deft_density_1000.evaluate_samples(dm_grid)
deft_sampled_2000 = deft_density_2000.evaluate_samples(dm_grid)
deft_sampled_4000 = deft_density_4000.evaluate_samples(dm_grid)
deft_sampled_10000 = deft_density_10000.evaluate_samples(dm_grid)
print('Sampled density calculations done.')

# Find std 
deft_std_frb_len = []
for i in range(len(dm_grid)):
    deft_std_frb_len_ = np.std(deft_sampled_frb_len[i])
    deft_std_frb_len = np.append(deft_std_frb_len,deft_std_frb_len_)

deft_std_100 = []
for i in range(len(dm_grid)):
    deft_std_100_ = np.std(deft_sampled_100[i])
    deft_std_100 = np.append(deft_std_100,deft_std_100_)

deft_std_1000 = []
for i in range(len(dm_grid)):
    deft_std_1000_ = np.std(deft_sampled_1000[i])
    deft_std_1000 = np.append(deft_std_1000,deft_std_1000_)

deft_std_2000 = []
for i in range(len(dm_grid)):
    deft_std_2000_ = np.std(deft_sampled_2000[i])
    deft_std_2000 = np.append(deft_std_2000,deft_std_2000_)

deft_std_4000 = []
for i in range(len(dm_grid)):
    deft_std_4000_ = np.std(deft_sampled_4000[i])
    deft_std_4000 = np.append(deft_std_4000,deft_std_4000_)
print('Standard deviation calculations done.')

# Save DEFT datasets
np.save('kde_and_deft_data/deft_frb_len_optimal.npy',deft_optimal_frb_len)
np.save('kde_and_deft_data/deft_100_optimal.npy',deft_optimal_100)
np.save('kde_and_deft_data/deft_1000_optimal.npy',deft_optimal_1000)
np.save('kde_and_deft_data/deft_2000_optimal.npy',deft_optimal_2000)
np.save('kde_and_deft_data/deft_4000_optimal.npy',deft_optimal_4000)
np.save('kde_and_deft_data/deft_10000_optimal.npy',deft_optimal_10000)

np.save('kde_and_deft_data/deft_frb_len_sampled.npy',deft_sampled_frb_len)
np.save('kde_and_deft_data/deft_100_sampled.npy',deft_sampled_100)
np.save('kde_and_deft_data/deft_1000_sampled.npy',deft_sampled_1000)
np.save('kde_and_deft_data/deft_2000_sampled.npy',deft_sampled_2000)
np.save('kde_and_deft_data/deft_4000_sampled.npy',deft_sampled_4000)
np.save('kde_and_deft_data/deft_10000_sampled.npy',deft_sampled_10000)


np.save('kde_and_deft_data/deft_std_frb_len.npy',deft_std_frb_len)
np.save('kde_and_deft_data/deft_std_100.npy',deft_std_100)
np.save('kde_and_deft_data/deft_std_1000.npy',deft_std_1000)
np.save('kde_and_deft_data/deft_std_2000.npy',deft_std_2000)
np.save('kde_and_deft_data/deft_std_4000.npy',deft_std_4000)
