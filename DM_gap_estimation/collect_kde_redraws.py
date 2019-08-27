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
from DM_redraw_definitions import DM_known_draws, FRB_distribution_from_SFR, make_pdf, make_kde_funtion
from DM_redraw_kde_gen import DM_grid
import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

matplotlib.rc('xtick', labelsize=30) 
matplotlib.rc('ytick', labelsize=30) 

plt.figure(figsize=[20,15])

"100 DRAWS"
DM_100_bootstrapped = []
numpy_vars_dist = {}
for np_name in glob.glob('DM_gap_outputs/DM_redraws/100/*.npy'):
    numpy_vars_dist[np_name] = np.load(np_name)
    DM_100_bootstrapped.append([numpy_vars_dist[np_name]])
    # plt.plot(DM_grid,numpy_vars_dist[np_name],linewidth=1,color='dodgerblue',alpha=.2)

all_kde_predictions = np.array(DM_100_bootstrapped)
kde_std_100 = np.std(all_kde_predictions, axis=0).reshape(-1)

"1000 DRAWS"
DM_1000_bootstrapped = []
numpy_vars_dist_1000 = {}
for np_name in glob.glob('DM_gap_outputs/DM_redraws/1000/*.npy'):
    numpy_vars_dist_1000[np_name] = np.load(np_name)
    DM_1000_bootstrapped.append([numpy_vars_dist_1000[np_name]])
    # plt.plot(DM_grid,numpy_vars_dist[np_name],linewidth=1,color='dodgerblue',alpha=.2)

all_kde_predictions_1000 = np.array(DM_1000_bootstrapped)
kde_std_1000 = np.std(all_kde_predictions_1000, axis=0).reshape(-1)