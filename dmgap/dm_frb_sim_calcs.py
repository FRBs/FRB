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

# Save DEFT datasets
np.save('sim_output_data/frb_100_optimal.npy',optimal_100)
np.save('sim_output_data/frb_1000_optimal.npy',optimal_1000)
np.save('sim_output_data/frb_10000_optimal.npy',optimal_10000)

np.save('sim_output_data/frb_100_sampled.npy',sampled_100)
np.save('sim_output_data/frb_1000_sampled.npy',sampled_1000)
np.save('sim_output_data/frb_10000_sampled.npy',sampled_10000)
