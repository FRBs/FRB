import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import random
import suftware as sw

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
psrcat_df = pd.read_csv('transient_data/psrcat_df.csv')

dm_frb_sim = np.load('sim_output_data/dm_frb_sim.npy')
dm_grid = np.load('sim_output_data/dm_grid.npy')

alpha_param = 3
b_box = [-600,6000]
num_ensembles = 1000

sample_100 = np.asarray(random.sample(list(dm_frb_sim),100))
density_100 = sw.DensityEstimator(sample_100, alpha=alpha_param, bounding_box=b_box, num_posterior_samples=num_ensembles)
sample_1000 = np.asarray(random.sample(list(dm_frb_sim),1000))
density_1000 = sw.DensityEstimator(sample_1000, alpha=alpha_param, bounding_box=b_box, num_posterior_samples=num_ensembles)
sample_10000 = np.asarray(random.sample(list(dm_frb_sim),10000))
density_10000 = sw.DensityEstimator(sample_10000, alpha=alpha_param, bounding_box=b_box, num_posterior_samples=num_ensembles)

# Evaluate optimal density
optimal_100 = density_100.evaluate(dm_grid)
optimal_1000 = density_1000.evaluate(dm_grid)
optimal_10000 = density_10000.evaluate(dm_grid)
print('Density calculations done.')

# Evaluate sampled densities
sampled_100 = density_100.evaluate_samples(dm_grid)
sampled_1000 = density_1000.evaluate_samples(dm_grid)
sampled_10000 = density_10000.evaluate_samples(dm_grid)
print('Sampled density calculations done.')

# Evaluate standard deviations
std_100 = []
for i in range(len(dm_grid)):
    std_100_ = np.std(sampled_100[i])
    std_100 = np.append(std_100,std_100_)

std_1000 = []
for i in range(len(dm_grid)):
    std_1000_ = np.std(sampled_1000[i])
    std_1000 = np.append(std_1000,std_1000_)

std_10000 = []
for i in range(len(dm_grid)):
    std_10000_ = np.std(sampled_10000[i])
    std_10000 = np.append(std_10000,std_10000_)

print('Standard deviation calculations done.')

# Save DEFT datasets
np.save('sim_output_data/optimal_100.npy',optimal_100)
np.save('sim_output_data/optimal_1000.npy',optimal_1000)
np.save('sim_output_data/optimal_10000.npy',optimal_10000)

np.save('sim_output_data/sampled_100.npy',sampled_100)
np.save('sim_output_data/sampled_1000.npy',sampled_1000)
np.save('sim_output_data/sampled_10000.npy',sampled_10000)

np.save('sim_output_data/std_100.npy',std_100)
np.save('sim_output_data/std_1000.npy',std_1000)
np.save('sim_output_data/std_10000.npy',std_10000)
