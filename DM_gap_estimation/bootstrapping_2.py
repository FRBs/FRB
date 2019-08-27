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
from DM_definitions import DM_known_draws, make_pdf, make_kde_funtion
# from DM_kde_gen import DM_frb, DM_rand_draws
import suftware as sw
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

rcParams['figure.figsize'] = 22,11
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 

kde_original = np.load('DM_outputs/Final_kernel/Final_kernel_choice_2000.npy')
DM_grid = np.load('DM_outputs/Final_kernel/DM_grid.npy')
DM_kde_draws_rand = np.load('DM_outputs/Final_kernel/DM_kde_rand_draws.npy')


# plt.plot(DM_grid,kde_original, linewidth=1, alpha=1, label=r'DM KDE', color='purple')

n = 2 #how many plots
reps = 10000 #how many samples
x = np.random.choice(DM_kde_draws_rand,(n, reps), replace=True)
print(np.shape(x))

for i in x:
    plt.hist(i, bins=100,alpha=0.5)
    # print(i)
    

plt.show()
#get percentile
# mb = np.percentile(x, 1, axis=0)
# mb.sort()
# lower, upper = np.percentile(mb, [2.5, 97.5])

# for v in (lower, upper):
#     plt.axvline(v, color='red',linewidth=1)

# y = 1/np.arange(1, n+1)[:, None]*np.cumsum(x, axis=0) #running average
# upper, lower = np.percentile(y, [2.5, 97.5], axis=1)
# plt.plot(np.arange(1, n+1)[:, None], y, c='grey', alpha=0.5,linewidth=0.5)
# plt.plot(np.arange(1, n+1), y[:, 0], c='red', linewidth=0.8)
# plt.plot(np.arange(1, n+1), upper, 'b', np.arange(1, n+1), lower, 'b',linewidth=0.8)

# Create figure
# plt.xlabel('$DM$', fontsize=30)
# plt.ylabel('PDF', fontsize=30)
# plt.xlim(0,1500)
# plt.ylim(0,0.0002)
# plt.tight_layout()
# plt.savefig('bootstrap_outputs/monte_int_multi_plots.png')
# plt.show()