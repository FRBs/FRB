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

# Sample size = 5000
dm_sample_5000 = np.sort(np.asarray(random.sample(list(dm_frb_sim),5000)))
n_5000 = 100
reps_5000 = 100
resample_5000 = np.sort(np.random.choice(dm_frb_sim, (n_5000, reps_5000),replace=True))
kde_optimal_5000 = make_kde_funtion(grid=dm_grid, draws = dm_sample_5000, min_bandwidth=20, max_bandwidth=30, bandwidth_stepsize=1, cv=5, kernel='gaussian')

kde_5000_bootstrapped = []
for i in range(len(frbcat_df['dmdiff'])):
    kde_predictions = make_kde_funtion(grid=dm_grid, draws = resample_5000[i], min_bandwidth=20, max_bandwidth=30, bandwidth_stepsize=1, cv=5, kernel='gaussian')
    kde_5000_bootstrapped.append([kde_predictions])

kde_5000_bootstrapped = np.array(kde_5000_bootstrapped)
kde_std_5000 = np.std(kde_5000_bootstrapped, axis=0).reshape(-1)

np.save('kde_and_deft_data/kde_5000_optimal.npy',kde_optimal_5000)
np.save('kde_and_deft_data/kde_5000_bootstrapped.npy',kde_5000_bootstrapped)
np.save('kde_and_deft_data/kde_std_5000.npy',kde_std_5000)
print('KDE done for 5000 draws')

# Sample size = 10000
dm_sample_10000 = np.sort(np.asarray(random.sample(list(dm_frb_sim),10000)))
n_10000 = 100
reps_10000 = 100
resample_10000 = np.sort(np.random.choice(dm_frb_sim, (n_10000, reps_10000),replace=True))
kde_optimal_10000 = make_kde_funtion(grid=dm_grid, draws = dm_sample_10000, min_bandwidth=20, max_bandwidth=30, bandwidth_stepsize=1, cv=5, kernel='gaussian')

kde_10000_bootstrapped = []
for i in range(len(frbcat_df['dmdiff'])):
    kde_predictions = make_kde_funtion(grid=dm_grid, draws = resample_10000[i], min_bandwidth=20, max_bandwidth=30, bandwidth_stepsize=1, cv=5, kernel='gaussian')
    kde_10000_bootstrapped.append([kde_predictions])

kde_10000_bootstrapped = np.array(kde_10000_bootstrapped)
kde_std_10000 = np.std(kde_10000_bootstrapped, axis=0).reshape(-1)

np.save('kde_and_deft_data/kde_10000_optimal.npy',kde_optimal_10000)
np.save('kde_and_deft_data/kde_10000_bootstrapped.npy',kde_10000_bootstrapped)
np.save('kde_and_deft_data/kde_std_10000.npy',kde_std_10000)
print('KDE done for 10000 draws')




