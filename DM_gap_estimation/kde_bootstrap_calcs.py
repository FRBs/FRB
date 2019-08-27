# PACKAGE IMPORT
import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import random
import matplotlib
from matplotlib import rcParams
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from DM_definitions import DM_known_draws, make_kde_funtion
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')
# rcParams['figure.figsize'] = 22,11
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 

random.seed(1919)

"DM_FRB AND SAMPLES"
num_kde_samples = 10*6
DM_frb = np.load('DM_outputs/DM_frb.npy')
DM_grid = np.arange(len(DM_frb))
DM_stepsize = 0.1
num_frb_draws = 1000
DM_rand_draws = np.asarray(random.sample(list(DM_frb),num_frb_draws))

"100 RANDOM draws"
DM_rand_100 = np.asarray(random.sample(list(DM_frb),100))
n_100 = len(DM_rand_100)
reps = 100
resample_rand_100 = np.random.choice(DM_rand_100, (n_100, reps))

#original draws
DM_kde_rand_100 = make_kde_funtion(grid=DM_grid, draws = DM_rand_100, min_bandwidth=25, max_bandwidth=100, bandwidth_stepsize=1, cv=30, kernel='gaussian')

#resample draws
plt.figure(figsize=[20,15])
DM_100_bootstrapped = []
for i in range(len(resample_rand_100)):
    kde_predictions = make_kde_funtion(grid=DM_grid, draws = resample_rand_100[i], min_bandwidth=25, max_bandwidth=100, bandwidth_stepsize=1, cv=30, kernel='gaussian')
    DM_100_bootstrapped.append([kde_predictions])
    plt.plot(DM_grid, kde_predictions, color='dodgerblue', alpha=.2, linewidth=1.5)
"Plot bootstrapped kde"
plt.plot(DM_grid, DM_kde_rand_100, color='blue',label='KDE 100 draws', linewidth=1)
plt.xlabel('$DM$', fontsize=30)
plt.ylabel('PDF', fontsize=30)
plt.xlim(0,1800)
plt.tight_layout()
plt.savefig('bootstrap_outputs/KDE_bootstrapped/KDE_100.png')
plt.show()

all_kde_predictions_100 = np.array(DM_100_bootstrapped)
kde_std_100 = np.std(all_kde_predictions_100, axis=0).reshape(-1)

np.save('kde_std_data/std_100.npy',kde_std_100)
np.save('kde_std_data/kde_100.npy',DM_kde_rand_100)

print('100 DONE')

"1000 RANDOM DRAWS"
DM_rand_1000 = np.asarray(random.sample(list(DM_frb),1000))
n_1000 = len(DM_rand_1000)
reps = 100
resample_rand_1000 = np.random.choice(DM_rand_1000, (n_1000, reps))

#original draws
kde_1000 = make_kde_funtion(grid=DM_grid, draws = DM_rand_1000, min_bandwidth=25, max_bandwidth=100, bandwidth_stepsize=1, cv=30, kernel='gaussian')

#resample draws
plt.figure(figsize=[20,15])
DM_1000_bootstrapped = []
for i in range(len(resample_rand_1000)):
    kde_predictions = make_kde_funtion(grid=DM_grid, draws = resample_rand_1000[i], min_bandwidth=25, max_bandwidth=100, bandwidth_stepsize=1, cv=30, kernel='gaussian')
    DM_1000_bootstrapped.append([kde_predictions])
    plt.plot(DM_grid, kde_predictions,color='dodgerblue', alpha=.2, linewidth=1.5)
"Plot bootstrapped kde"
plt.plot(DM_grid, kde_1000, color='blue',label='KDE 1000 draws', linewidth=1)
plt.xlabel('$DM$', fontsize=30)
plt.ylabel('PDF', fontsize=30)
plt.xlim(0,1800)
plt.tight_layout()
plt.savefig('bootstrap_outputs/KDE_bootstrapped/KDE_1000.png')
plt.show()

all_kde_predictions_1000 = np.array(DM_1000_bootstrapped)
kde_std_1000 = np.std(all_kde_predictions_1000, axis=0).reshape(-1)

np.save('kde_std_data/std_1000.npy',kde_std_1000)
np.save('kde_std_data/kde_1000.npy',kde_1000)

print('1000 DONE')

"OBSERVED"
n_observed = len(DM_known_draws)
reps_observed = 100
resample_observed = np.random.choice(DM_known_draws, (n_observed, reps_observed))

#original draws
DM_kde_known = make_kde_funtion(grid=DM_grid, draws = DM_known_draws, min_bandwidth=25, max_bandwidth=100, bandwidth_stepsize=1, cv=30, kernel='gaussian')

#resample draws
plt.figure(figsize=[20,15])
DM_observed_bootstrapped = []
for i in range(len(resample_observed)):
    kde_predictions = make_kde_funtion(grid=DM_grid, draws = resample_observed[i], min_bandwidth=25, max_bandwidth=100, bandwidth_stepsize=1, cv=30, kernel='gaussian')
    DM_observed_bootstrapped.append([kde_predictions])
    plt.plot(DM_grid, kde_predictions,color='dodgerblue', alpha=.2, linewidth=1.5)
"Plot bootstrapped kde"
plt.plot(DM_grid, DM_kde_known, color='blue')
plt.xlabel('$DM$', fontsize=30)
plt.ylabel('PDF', fontsize=30)
plt.xlim(0,1800)
plt.tight_layout()
plt.savefig('bootstrap_outputs/KDE_bootstrapped/KDE_observed.png')
plt.show()

all_kde_predictions_obs = np.array(DM_observed_bootstrapped)
kde_std_observed = np.std(all_kde_predictions_obs, axis=0).reshape(-1)

np.save('kde_std_data/std_observed.npy', kde_std_observed)
np.save('kde_std_data/kde_observed.npy', DM_kde_known)

print('OBSERVED DONE')


