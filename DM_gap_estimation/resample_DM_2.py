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
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from DM_definitions_boot import DM_known_draws, FRB_distribution_from_SFR, make_pdf, make_kde_funtion
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

"RAND 100"
DM_rand_100 = np.asarray(random.sample(list(DM_frb),100))
n_100 = len(DM_rand_100)
reps = 100
resample_rand_100 = np.random.choice(DM_rand_100, (n_100, reps))

plt.figure(figsize=[20,15])
var100 = []
for i in range(len(resample_rand_100)):
    DM_kde_rand_func_ = make_kde_funtion(grid=DM_grid, draws = resample_rand_100[i], bandwidth=60, kernel='gaussian')
    var100_ = np.var(DM_kde_rand_func_[i])
    var100 = np.append(var100,var100_)    
    # plt.plot(DM_grid, DM_kde_rand_func_,color='dodgerblue', alpha=.2,linewidth=1.5)
DM_kde_rand_100 = make_kde_funtion(grid=DM_grid, draws = DM_rand_100, bandwidth=60, kernel='gaussian')
# plt.plot(DM_grid, DM_kde_rand_100, color='blue',label='KDE 100 draws',linewidth=1)
# plt.xlabel('$DM$', fontsize=30)
# plt.ylabel('PDF', fontsize=30)
# plt.xlim(0,200)
# plt.tight_layout()
# plt.savefig('bootstrap_outputs/KDE_bootstrapped/sim_KDE_100_cropped.png')
# plt.show()

DM_var_func_100 = []
for i in range(len(resample_rand_100)):
    kde_predictions = make_kde_funtion(grid=DM_grid, draws = resample_rand_100[i], bandwidth=60, kernel='gaussian')
    DM_var_func_100.append([kde_predictions])

all_kde_predictions_100 = np.array(DM_var_func_100)
var_by_dm_100 = np.std(all_kde_predictions_100, axis=0).reshape(-1)

np.save('variance_plots/var100kde.npy',var_by_dm_100)
np.save('variance_plots/100kde.npy',DM_kde_rand_100)

print('100 DONE')


"RAND 1000"
DM_rand_1000 = np.asarray(random.sample(list(DM_frb),1000))
n_1000 = len(DM_rand_1000)
reps = 100
resample_rand_1000 = np.random.choice(DM_rand_1000, (n_1000, reps))
plt.figure(figsize=[20,15])
var1000 = []
for i in range(len(resample_rand_1000)):
    DM_kde_rand_func_ = make_kde_funtion(grid=DM_grid, draws = resample_rand_1000[i], bandwidth=60, kernel='gaussian')
    var1000_ = np.var(DM_kde_rand_func_[i])
    var1000 = np.append(var1000,var1000_)    
    # plt.plot(DM_grid, DM_kde_rand_func_,color='dodgerblue', alpha=.1,linewidth=1.5)
DM_kde_rand_1000 = make_kde_funtion(grid=DM_grid, draws = DM_rand_1000, bandwidth=60, kernel='gaussian')
# plt.plot(DM_grid, DM_kde_rand_1000, color='purple',label='KDE 1000 draws', linewidth=1)
# plt.xlabel('$DM$', fontsize=30)
# plt.ylabel('PDF', fontsize=30)
# plt.xlim(0,200)
# plt.tight_layout()
# plt.savefig('bootstrap_outputs/KDE_bootstrapped/sim_KDE_1000_cropped.png')
# plt.show()

DM_var_func_1000 = []
for i in range(len(resample_rand_1000)):
    kde_predictions = make_kde_funtion(grid=DM_grid, draws = resample_rand_1000[i], bandwidth=60, kernel='gaussian')
    DM_var_func_1000.append([kde_predictions])

all_kde_predictions_1000 = np.array(DM_var_func_1000)
var_by_dm_1000 = np.std(all_kde_predictions_1000, axis=0).reshape(-1)

np.save('variance_plots/var1000kde.npy',var_by_dm_1000)
np.save('variance_plots/1000kde.npy',DM_kde_rand_1000)

print('1000 DONE')

"KNOWN"
n_known = len(DM_known_draws)
reps_known = 100
resample_known = np.random.choice(DM_known_draws, (n_known, reps_known))
plt.figure(figsize=[20,15])
varobs = []
for i in range(len(resample_known)):
    DM_kde_known_func_ = make_kde_funtion(grid=DM_grid, draws = resample_known[i], bandwidth=80, kernel='gaussian')
    varobs_ = np.var(DM_kde_known_func_[i])
    varobs = np.append(varobs,varobs_)   
    plt.plot(DM_grid, DM_kde_known_func_,color='dodgerblue', alpha=.2,linewidth=1.5)
DM_kde_known = make_kde_funtion(grid=DM_grid, draws = DM_known_draws, bandwidth=80, kernel='gaussian')
# plt.plot(DM_grid, DM_kde_known, color='blue')
# plt.xlabel('$DM$', fontsize=30)
# plt.ylabel('PDF', fontsize=30)
# plt.xlim(0,200)
# plt.tight_layout()
# plt.savefig('bootstrap_outputs/KDE_bootstrapped/sim_KDE_obs_cropped.png')
# plt.show()

DM_var_func_obs = []
for i in range(len(resample_known)):
    kde_predictions = make_kde_funtion(grid=DM_grid, draws = resample_known[i], bandwidth=60, kernel='gaussian')
    DM_var_func_obs.append([kde_predictions])

all_kde_predictions_obs = np.array(DM_var_func_obs)
var_by_dm_obs = np.std(all_kde_predictions_obs, axis=0).reshape(-1)

np.save('variance_plots/varobskde.npy',var_by_dm_obs)
np.save('variance_plots/obskde.npy', DM_kde_known)

print('KNOWN DONE')


