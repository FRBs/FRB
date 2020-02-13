import numpy as np
import random
import pandas as pd
from pandas import DataFrame as df
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from pdf_defs import make_kde_funtion
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 

random.seed(1919)

# Data upload 
dm_frb_sim = np.load('sim_output_data/dm_frb_sim.npy')
dm_grid = np.load('sim_output_data/dm_grid.npy')

dm_sample_100 = np.sort(np.asarray(random.sample(list(dm_frb_sim),100)))
dm_sample_1000 = np.sort(np.asarray(random.sample(list(dm_frb_sim),1000)))
dm_sample_5000 = np.sort(np.asarray(random.sample(list(dm_frb_sim),5000)))
dm_sample_10000 = np.sort(np.asarray(random.sample(list(dm_frb_sim),10000)))

num_resamples = 50
with NumpyRNGContext(1):
    boot_100 = bootstrap(dm_sample_100, num_resamples)
    boot_1000 = bootstrap(dm_sample_1000, num_resamples)
    boot_5000 = bootstrap(dm_sample_5000, num_resamples)
    boot_10000 = bootstrap(dm_sample_10000, num_resamples)

# Sample size = 100
kde_100 = []
for i in range(num_resamples):
    kde_func = make_kde_funtion(grid=dm_grid, draws = boot_100[i], min_bandwidth=40, max_bandwidth=100, bandwidth_stepsize=5, cv=5, kernel='gaussian')
    kde_100.append(kde_func)
kde_100 = np.array(kde_100)
kde_100_opt =  make_kde_funtion(grid=dm_grid, draws = dm_sample_100, min_bandwidth=40, max_bandwidth=100, bandwidth_stepsize=5, cv=5, kernel='gaussian') #np.mean(kde_100, axis=0).reshape(-1)
kde_100_std = np.std(kde_100, axis=0).reshape(-1)

np.save('kde_output_data/kde_100.npy',kde_100)
np.save('kde_output_data/kde_100_opt.npy',kde_100_opt)
np.save('kde_output_data/kde_100_std.npy',kde_100_std)
np.save('kde_output_data/dm_sample_100.npy',dm_sample_100)
print('KDE done for 100 draws')

# # Sample size = 1000
# kde_1000 = []
# for i in range(num_resamples):
#     kde_func = make_kde_funtion(grid=dm_grid, draws = boot_1000[i], min_bandwidth=10, max_bandwidth=80, bandwidth_stepsize=5, cv=5, kernel='gaussian')
#     kde_1000.append(kde_func)
# kde_1000 = np.array(kde_1000)
# kde_1000_opt = np.mean(kde_1000, axis=0).reshape(-1)
# kde_1000_std = np.std(kde_1000, axis=0).reshape(-1)

# np.save('kde_output_data/kde_1000.npy',kde_1000)
# np.save('kde_output_data/kde_1000_opt.npy',kde_1000_opt)
# np.save('kde_output_data/kde_1000_std.npy',kde_1000_std)
# np.save('kde_output_data/dm_sample_1000.npy',dm_sample_1000)
# print('KDE done for 1000 draws')

# # Sample size = 5000
# kde_5000 = []
# for i in range(num_resamples):
#     kde_func = make_kde_funtion(grid=dm_grid, draws = boot_5000[i], min_bandwidth=10, max_bandwidth=60, bandwidth_stepsize=5, cv=5, kernel='gaussian')
#     kde_5000.append(kde_func)
# kde_5000 = np.array(kde_5000)
# kde_5000_opt = np.mean(kde_5000, axis=0).reshape(-1)
# kde_5000_std = np.std(kde_5000, axis=0).reshape(-1)

# np.save('kde_output_data/kde_5000.npy',kde_5000)
# np.save('kde_output_data/kde_5000_opt.npy',kde_5000_opt)
# np.save('kde_output_data/kde_5000_std.npy',kde_5000_std)
# np.save('kde_output_data/dm_sample_5000.npy',dm_sample_5000)
# print('KDE done for 5000 draws')

# # Sample size = 10000
# kde_10000 = []
# for i in range(num_resamples):
#     kde_func = make_kde_funtion(grid=dm_grid, draws = boot_10000[i], min_bandwidth=5, max_bandwidth=10, bandwidth_stepsize=1, cv=5, kernel='gaussian')
#     kde_10000.append(kde_func)
# kde_10000 = np.array(kde_10000)
# kde_10000_opt = np.mean(kde_10000, axis=0).reshape(-1)
# kde_10000_std = np.std(kde_10000, axis=0).reshape(-1)

# np.save('kde_output_data/kde_10000.npy',kde_10000)
# np.save('kde_output_data/kde_10000_opt.npy',kde_10000_opt)
# np.save('kde_output_data/kde_10000_std.npy',kde_10000_std)
# np.save('kde_output_data/dm_sample_10000.npy',dm_sample_10000)
# print('KDE done for 10000 draws')

