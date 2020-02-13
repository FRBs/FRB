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

random.seed(1919)

# Data upload 
psrcat_df = pd.read_csv('transient_data/psrcat_df.csv')
grid_psr = np.arange(-500,500,0.001)
np.save('kde_output_data/grid_psr.npy', grid_psr)

frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.arange(-1000,6000,0.01)
np.save('kde_output_data/grid_frb.npy', grid_frb)

delta_dm_psr = psrcat_df['deltaDM']
delta_dm_frb = frbcat_df['deltaDM']

num_resamples = 100
with NumpyRNGContext(1):
    boot_obs_psr = bootstrap(delta_dm_psr, num_resamples)
    boot_obs_frb = bootstrap(delta_dm_frb, num_resamples)

kde_psr = []
for i in range(num_resamples):
    kde_func = make_kde_funtion(grid=grid_psr, draws = boot_obs_psr[i], min_bandwidth=5, max_bandwidth=15, bandwidth_stepsize=1, cv=5, kernel='gaussian')
    kde_psr.append(kde_func)
kde_psr = np.array(kde_psr)
kde_psr_opt = np.mean(kde_psr, axis=0).reshape(-1)
kde_psr_std = np.std(kde_psr, axis=0).reshape(-1)

np.save('kde_output_data/kde_psr.npy',kde_psr)
np.save('kde_output_data/kde_psr_opt.npy',kde_psr_opt)
np.save('kde_output_data/kde_psr_std.npy',kde_psr_std)
print('KDE done for pulsars.')

kde_frb = []
for i in range(num_resamples):
    kde_func = make_kde_funtion(grid=grid_frb, draws = boot_obs_frb[i], min_bandwidth=50, max_bandwidth=80, bandwidth_stepsize=10, cv=5, kernel='gaussian')
    kde_frb.append(kde_func)
kde_frb = np.array(kde_frb)
kde_frb_opt = np.mean(kde_frb, axis=0).reshape(-1)
kde_frb_std = np.std(kde_frb, axis=0).reshape(-1)

plt.plot(grid_frb, kde_frb_opt,color='crimson')
for i in range(len(kde_frb)):
    plt.plot(grid_frb,kde_frb[i],alpha=.2, color='crimson')
plt.show()

np.save('kde_output_data/kde_frb.npy',kde_frb)
np.save('kde_output_data/kde_frb_opt.npy',kde_frb_opt)
np.save('kde_output_data/kde_frb_std.npy',kde_frb_std)
print('KDE done for FRBs.')

