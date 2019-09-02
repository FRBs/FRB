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
from DM_kde_gen import DM_grid
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')
# rcParams['figure.figsize'] = 22,11
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 

# Upload optimal DEFT
deft_optimal_100 = np.load('kde_and_deft_data/deft_100_optimal.npy')
deft_optimal_1000 = np.load('kde_and_deft_data/deft_1000_optimal.npy')
deft_optimal_observed = np.load('kde_and_deft_data/deft_observed_optimal.npy')

# Upload other viable DEFTs
deft_sampled_100 = np.load('kde_and_deft_data/deft_100_sampled.npy')
deft_sampled_1000 = np.load('kde_and_deft_data/deft_1000_sampled.npy')
deft_sampled_observed = np.load('kde_and_deft_data/deft_observed_sampled.npy')

plt.figure(figsize=[20,15])
plt.plot(DM_grid, deft_optimal_100, color='blue',linewidth=1)
plt.plot(DM_grid, deft_sampled_100, color='dodgerblue', alpha=.2, linewidth=1.5)
plt.xlim(0,150)
plt.xlabel('$DM$')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('bootstrap_outputs/DEFT/DEFT_100_cropped.png')
plt.show()

plt.figure(figsize=[20,15])
plt.plot(DM_grid, deft_optimal_1000, color='blue',linewidth=1)
plt.plot(DM_grid, deft_sampled_1000, color='dodgerblue', alpha=.2, linewidth=1.5)
plt.xlim(0,150)
plt.xlabel('$DM$')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('bootstrap_outputs/DEFT/DEFT_1000_cropped.png')
plt.show()

plt.figure(figsize=[20,15])
plt.plot(DM_grid, deft_optimal_observed, color='blue',linewidth=1)
plt.plot(DM_grid, deft_sampled_observed, color='dodgerblue', alpha=.2, linewidth=1.5)
plt.xlim(0,150)
plt.xlabel('$DM$')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('bootstrap_outputs/DEFT/DEFT_observed_cropped.png')
plt.show()



