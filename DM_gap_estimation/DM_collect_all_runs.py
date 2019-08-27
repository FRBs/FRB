import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import sklearn
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
from DM_definitions import DM_known_draws, FRB_distribution_from_SFR, make_pdf, make_kde_funtion
from DM_kde_gen import DM_grid
import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 

#numpy_vars_hist = {}
#for np_name in glob.glob('DM_gap_outputs/DM_rand_npy/rand_hist/*.npy'):
#    numpy_vars_hist[np_name] = np.load(np_name)
#    print(min(numpy_vars_hist[np_name]))

numpy_vars_dist = {}
for np_name in glob.glob('DM_gap_outputs/DM_rand_npy/rand_distri/*.npy'):
    numpy_vars_dist[np_name] = np.load(np_name)
    plt.plot(DM_grid,numpy_vars_dist[np_name],linewidth=0.5)

plt.xlabel('$DM$', fontsize=14)
plt.ylabel('PDF', fontsize=14)
plt.tight_layout()
plt.xlim(0,1500)
plt.savefig('DM_gap_outputs/DM_rand_npy/rand_distri/opt_plot.png')
plt.show()
