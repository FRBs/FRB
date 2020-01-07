import sys
sys.path.insert(1, '/eplatts_UCSC_Server/dm_gap/ne2001-master/src')
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw
from percentile_defs import find_quantile
import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 

# Data upload 
grid_frb = np.arange(0,6000,0.01)
grid_psr = np.arange(-500,500,0.01)

dmdiff_psr_optimal = np.load('obs_outputs/dmdiff_psr_optimal.npy')
dmdiff_frb_optimal = np.load('obs_outputs/dmdiff_frb_optimal.npy')

ci_interval = [.95, .96, .97, .98, .99] #confidence
ci_cut_off = np.sqrt(ci_interval) #for joint probability: probability that <ci_interval% of events will occur simultaneously

frb_obs_dm = []
for i in range(len(ci_cut_off)):
    frb_obs_dm_ = find_quantile(grid=grid_frb,distribution=dmdiff_frb_optimal,quantile=1-ci_cut_off[i])
    frb_obs_dm = np.append(frb_obs_dm,frb_obs_dm_)

psr_optimal_dm = find_quantile(grid=grid_psr,distribution=dmdiff_psr_optimal,quantile=ci_cut_off)

print(frb_obs_dm)
print(psr_optimal_dm)

#!!save to dataframe