# PACKAGE IMPORT
import numpy as np
import scipy as sp
import random
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import statistics as stats
import pandas as pd
from pandas import DataFrame as df
from pdf_defs import make_kde_funtion, find_kde_ensemble

import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 
 
random.seed(1919)

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')

dm_frb_sim = np.load('sim_outputs/dm_frb_sim.npy')
dm_grid = np.load('sim_outputs/dm_grid.npy')

# Sample size = 100
dm_sample_100 = np.sort(np.asarray(random.sample(list(dm_frb_sim),100)))
n_100 = 100
reps_100 = 100
resample_100 = np.sort(np.random.choice(dm_frb_sim, (n_100, reps_100),replace=True))
kde_optimal_100 = make_kde_funtion(grid=dm_grid, draws = dm_sample_100, min_bandwidth=50, max_bandwidth=100, bandwidth_stepsize=1, cv=10, kernel='gaussian')

kde_100_bootstrapped = []
for i in range(len(frbcat_df['dmdiff'])):
    kde_predictions = make_kde_funtion(grid=dm_grid, draws = resample_100[i], min_bandwidth=50, max_bandwidth=100, bandwidth_stepsize=1, cv=5, kernel='gaussian')
    kde_100_bootstrapped.append([kde_predictions])

kde_100_bootstrapped = np.array(kde_100_bootstrapped)
kde_std_100 = np.std(kde_100_bootstrapped, axis=0).reshape(-1)

np.save('kde_and_deft_data/kde_100_optimal.npy',kde_optimal_100)
np.save('kde_and_deft_data/kde_100_bootstrapped.npy',kde_100_bootstrapped)
np.save('kde_and_deft_data/kde_std_100.npy',kde_std_100)
print('KDE done for 100 draws')

# Sample size = 1000
dm_sample_1000 = np.sort(np.asarray(random.sample(list(dm_frb_sim),1000)))
n_1000 = 100
reps_1000 = 100
resample_1000 = np.sort(np.random.choice(dm_frb_sim, (n_1000, reps_1000),replace=True))
kde_optimal_1000 = make_kde_funtion(grid=dm_grid, draws = dm_sample_1000, min_bandwidth=30, max_bandwidth=50, bandwidth_stepsize=1, cv=5, kernel='gaussian')

kde_1000_bootstrapped = []
for i in range(len(frbcat_df['dmdiff'])):
    kde_predictions = make_kde_funtion(grid=dm_grid, draws = resample_1000[i], min_bandwidth=30, max_bandwidth=50, bandwidth_stepsize=1, cv=5, kernel='gaussian')
    kde_1000_bootstrapped.append([kde_predictions])

kde_1000_bootstrapped = np.array(kde_1000_bootstrapped)
kde_std_1000 = np.std(kde_1000_bootstrapped, axis=0).reshape(-1)

np.save('kde_and_deft_data/kde_1000_optimal.npy',kde_optimal_1000)
np.save('kde_and_deft_data/kde_1000_bootstrapped.npy',kde_1000_bootstrapped)
np.save('kde_and_deft_data/kde_std_1000.npy',kde_std_1000)
print('KDE done for 1000 draws')

# Sample size = 2000
dm_sample_2000 = np.sort(np.asarray(random.sample(list(dm_frb_sim),2000)))
n_2000 = 100
reps_2000 = 100
resample_2000 = np.sort(np.random.choice(dm_frb_sim, (n_2000, reps_2000),replace=True))
kde_optimal_2000 = make_kde_funtion(grid=dm_grid, draws = dm_sample_2000, min_bandwidth=20, max_bandwidth=30, bandwidth_stepsize=1, cv=5, kernel='gaussian')

kde_2000_bootstrapped = []
for i in range(len(frbcat_df['dmdiff'])):
    kde_predictions = make_kde_funtion(grid=dm_grid, draws = resample_2000[i], min_bandwidth=20, max_bandwidth=30, bandwidth_stepsize=1, cv=5, kernel='gaussian')
    kde_2000_bootstrapped.append([kde_predictions])

kde_2000_bootstrapped = np.array(kde_2000_bootstrapped)
kde_std_2000 = np.std(kde_2000_bootstrapped, axis=0).reshape(-1)

np.save('kde_and_deft_data/kde_2000_optimal.npy',kde_optimal_2000)
np.save('kde_and_deft_data/kde_2000_bootstrapped.npy',kde_2000_bootstrapped)
np.save('kde_and_deft_data/kde_std_2000.npy',kde_std_2000)
print('KDE done for 2000 draws')

# Sample size = 4000
dm_sample_4000 = np.sort(np.asarray(random.sample(list(dm_frb_sim),4000)))
n_4000 = 100
reps_4000 = 100
resample_4000 = np.sort(np.random.choice(dm_frb_sim, (n_4000, reps_4000),replace=True))
kde_optimal_4000 = make_kde_funtion(grid=dm_grid, draws = dm_sample_4000, min_bandwidth=20, max_bandwidth=30, bandwidth_stepsize=1, cv=5, kernel='gaussian')

kde_4000_bootstrapped = []
for i in range(len(frbcat_df['dmdiff'])):
    kde_predictions = make_kde_funtion(grid=dm_grid, draws = resample_4000[i], min_bandwidth=20, max_bandwidth=30, bandwidth_stepsize=1, cv=5, kernel='gaussian')
    kde_4000_bootstrapped.append([kde_predictions])

kde_4000_bootstrapped = np.array(kde_4000_bootstrapped)
kde_std_4000 = np.std(kde_4000_bootstrapped, axis=0).reshape(-1)

np.save('kde_and_deft_data/kde_4000_optimal.npy',kde_optimal_4000)
np.save('kde_and_deft_data/kde_4000_bootstrapped.npy',kde_4000_bootstrapped)
np.save('kde_and_deft_data/kde_std_4000.npy',kde_std_4000)
print('KDE done for 4000 draws')




