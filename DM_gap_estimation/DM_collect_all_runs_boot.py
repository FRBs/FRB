import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import sklearn
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
from DM_definitions_boot import DM_known_draws, FRB_distribution_from_SFR, make_pdf, make_kde_funtion
from DM_kde_gen_boot import DM_grid
import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

matplotlib.rc('xtick', labelsize=30) 
matplotlib.rc('ytick', labelsize=30) 

plt.figure(figsize=[20,15])
DM_var_func = []
numpy_vars_dist = {}
for np_name in glob.glob('DM_gap_outputs/DM_rand_boot/100/*.npy'):
    numpy_vars_dist[np_name] = np.load(np_name)
    DM_var_func.append([numpy_vars_dist[np_name]])
    # plt.plot(DM_grid,numpy_vars_dist[np_name],linewidth=1,color='dodgerblue',alpha=.2)


# DM_var_func_100 = []
# for i in range(len(resample_rand_100)):
#     kde_predictions = make_kde_funtion(grid=DM_grid, draws = resample_rand_100[i], bandwidth=60, kernel='gaussian')
#     DM_var_func_100.append([kde_predictions])

all_kde_predictions = np.array(DM_var_func)
var_by_dm = np.var(all_kde_predictions, axis=0).reshape(-1)
# plt.plot(var_by_dm)
# plt.xlim(0,2000)
# plt.show()

###########
DM_var_func_1000 = []
numpy_vars_dist_1000 = {}
for np_name in glob.glob('DM_gap_outputs/DM_rand_boot/*.npy'):
    numpy_vars_dist_1000[np_name] = np.load(np_name)
    DM_var_func_1000.append([numpy_vars_dist_1000[np_name]])
    # plt.plot(DM_grid,numpy_vars_dist[np_name],linewidth=1,color='dodgerblue',alpha=.2)

all_kde_predictions_1000 = np.array(DM_var_func_1000)
var_by_dm_1000 = np.var(all_kde_predictions_1000, axis=0).reshape(-1)
# plt.plot(var_by_dm_1000)
# plt.xlim(0,2000)
# plt.show()
###########


# plt.plot(DM_grid,numpy_vars_dist[np_name],linewidth=2,color='blue')
# plt.xlabel('$DM$', fontsize=30)
# plt.ylabel('PDF', fontsize=30)
# plt.xlim(0,200)
# plt.tight_layout()
# plt.savefig('bootstrap_outputs/KDE_bootstrapped/bootstrap1000_cropped.png')
# plt.show()

############
DM_var_func_1000 = []
numpy_vars_dist_1000 = {}
for np_name in glob.glob('DM_gap_outputs/DM_rand_boot/*.npy'):
    numpy_vars_dist_1000[np_name] = np.load(np_name)
    DM_var_func_1000.append([numpy_vars_dist_1000[np_name]])
    # plt.plot(DM_grid,numpy_vars_dist[np_name],linewidth=1,color='dodgerblue',alpha=.2)

all_kde_predictions_1000 = np.array(DM_var_func_1000)
var_by_dm_1000 = np.var(all_kde_predictions_1000, axis=0).reshape(-1)



#############

"KNOWN"
n_known = len(DM_known_draws)
reps_known = 1000
resample_known = np.random.choice(DM_known_draws, (n_known, reps_known))

varobs = []
for i in range(len(resample_known)):
    DM_kde_known_func_ = make_kde_funtion(grid=DM_grid, draws = resample_known[i], bandwidth=80, kernel='gaussian')
    varobs_ = np.var(DM_kde_known_func_[i])
    varobs = np.append(varobs,varobs_)   
#     plt.plot(DM_grid, DM_kde_known_func_,color='dodgerblue', alpha=.1)
# DM_kde_known = make_kde_funtion(grid=DM_grid, draws = DM_known_draws, bandwidth=80, kernel='gaussian')
# plt.plot(DM_grid, DM_kde_known, color='blue')
# plt.xlabel('$DM$', fontsize=30)
# plt.ylabel('PDF', fontsize=30)
# plt.xlim(0,200)
# plt.savefig('bootstrap_outputs/sim_KDE_cropped_obs.png')
# plt.show()

DM_var_func_obs = []
for i in range(len(resample_known)):
    kde_predictions = make_kde_funtion(grid=DM_grid, draws = resample_known[i], bandwidth=80, kernel='gaussian')
    DM_var_func_obs.append([kde_predictions])

all_kde_predictions_obs = np.array(DM_var_func_obs)
var_by_dm_obs = np.var(all_kde_predictions_obs, axis=0).reshape(-1)

#Variance
plt.figure(figsize=[20,15])
# plt.plot(DM_grid,var_by_dm_obs, linewidth=1.5, label='Observed',color='k',linestyle='dashed')
plt.plot(DM_grid,var_by_dm, linewidth=1.5, label='100 samples', color='purple')
plt.plot(DM_grid,var_by_dm_1000, linewidth=1.5, label='1000 samples',color='green')
plt.xlim(0,2000)
plt.xlabel('$DM$')
plt.ylabel('Variance')
plt.legend(fontsize=40)
plt.tight_layout()
plt.savefig('bootstrap_outputs/KDE_bootstrapped/variances_bootstrapped_KDE.png')
plt.show()