import numpy as np
import pandas as pd
import random
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw
from scipy.signal import argrelextrema
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from percentile_defs import find_percentile
from pdf_defs import make_kde_funtion

# Data upload 
psrcat_df = pd.read_csv('transient_data/psrcat_df.csv')
grid_psr = np.load('obs_output_data/grid_psr.npy')
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.load('obs_output_data/grid_frb.npy')

delta_dm_psr = psrcat_df['deltaDM']
delta_dm_frb = frbcat_df['deltaDM']

# Constraints to region of interest for max concavity
dm_min_frb = -100
dm_max_frb = 150
dm_min_frb_idx = (np.abs(grid_frb - dm_min_frb)).argmin()
dm_max_frb_idx = (np.abs(grid_frb - dm_max_frb)).argmin()

dm_min_psr = -10
dm_max_psr = 20
dm_min_psr_idx = (np.abs(grid_psr - dm_min_psr)).argmin()
dm_max_psr_idx = (np.abs(grid_psr - dm_max_psr)).argmin()

kde_psr = np.load('kde_output_data/kde_psr.npy')
kde_psr = kde_psr[:,dm_min_psr_idx:dm_max_psr_idx]
kde_psr_opt = np.load('kde_output_data/kde_psr_opt.npy')
kde_psr_opt = kde_psr_opt[dm_min_psr_idx:dm_max_psr_idx]

kde_frb = np.load('kde_output_data/kde_frb.npy')
kde_frb = kde_frb[:,dm_min_frb_idx:dm_max_frb_idx]
kde_frb_opt = np.load('kde_output_data/kde_frb_opt.npy')
kde_frb_opt = kde_frb_opt[dm_min_frb_idx:dm_max_frb_idx]

grid_psr = grid_psr[dm_min_psr_idx:dm_max_psr_idx] 
grid_frb = grid_frb[dm_min_frb_idx:dm_max_frb_idx]

# Pulsars
grid_psr_pp = grid_psr[2:]
psr_pp_opt = np.diff(np.diff(kde_psr_opt)) #second derive
psr_gap = grid_psr_pp[psr_pp_opt.argmax()]
print('Optimal DeltaDM_psr:',psr_gap)

psr_pp = np.diff(np.diff(kde_psr))
psr_gaps = []
for i in range(len(psr_pp)):
    psr_gaps_ = grid_psr_pp[psr_pp[i].argmax()]
    psr_gaps = np.append(psr_gaps,psr_gaps_)
# psr_gaps = np.append(psr_gaps,psr_gap)
print('Mean DeltaDM_psr:',np.mean(psr_gaps))
print('2 sigma:',2*np.std(psr_gaps))
kde_psr_func = make_kde_funtion(grid=grid_psr_pp, draws=psr_gaps, min_bandwidth=1, max_bandwidth=9, bandwidth_stepsize=1, cv=5, kernel='gaussian')

# FRBs
grid_frb_pp = grid_frb[2:]
frb_pp_opt = np.diff(np.diff(kde_frb_opt)) #second derive
frb_gap = grid_frb_pp[frb_pp_opt.argmax()]
print('Optimal DeltaDM_FRB:',frb_gap)

frb_pp = np.diff(np.diff(kde_frb))
frb_gaps = []
for i in range(len(frb_pp)):
    frb_gaps_ = grid_frb_pp[frb_pp[i].argmax()]
    frb_gaps = np.append(frb_gaps,frb_gaps_)
# frb_gaps = np.append(frb_gaps,frb_gap)
print('Mean DeltaDM_FRB:',np.mean(frb_gaps))
print('2 sigma:',2*np.std(frb_gaps))
kde_frb_func = make_kde_funtion(grid=grid_frb_pp, draws=frb_gaps, min_bandwidth=10, max_bandwidth=14, bandwidth_stepsize=1, cv=5, kernel='gaussian')

plt.plot(grid_psr_pp,kde_psr_func,color='teal',label=r'max(DM$_\mathrm{pulsar}$)')
plt.hist(psr_gaps,density=True, bins=50, histtype='stepfilled', alpha=.2, color='teal')
plt.axvline(x=np.mean(psr_gaps),color='teal',linestyle='--',linewidth=1,label='Mean pulsar DM')
plt.plot(grid_frb_pp,kde_frb_func,color='crimson',label=r'min(DM$_\mathrm{FRB}$)')
plt.hist(frb_gaps,density=True, bins=80, histtype='stepfilled', alpha=.2, color='crimson')
plt.axvline(x=np.mean(frb_gaps),color='red',linestyle='--',linewidth=1,label='Mean FRB DM')
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=12)
plt.ylabel('PDF',fontsize=12)
plt.xlim(-50,100)
plt.ylim(0,0.11)
plt.legend(fontsize=12,loc=1)
plt.tight_layout()
plt.savefig('kde_output_figs/DM_frb_psr.png', dpi=300)
plt.show()

# plt.plot(grid_psr_pp,psr_pp_opt)
# for i in range(len(psr_pp)):
#     plt.plot(grid_psr_pp,psr_pp[i],color='#f46036',linewidth=1.5,alpha=.2)
# plt.show()

# plt.plot(grid_frb_pp,frb_pp_opt)
# for i in range(len(frb_pp)):
#     plt.plot(grid_frb_pp,frb_pp[i],color='#f46036',linewidth=1.5,alpha=.2)
# plt.show()

# # Gap
# gap_obs = []
# for i in range(len(frb_gaps)):
#     gap_obs_ = frb_gaps[i]-psr_gaps
#     gap_obs = np.append(gap_obs,gap_obs_)

# # num_resamples=2
# # with NumpyRNGContext(1):
# #     boot_obs_gap = bootstrap(gap_obs, num_resamples)

# gap_obs = np.sort(np.asarray(random.sample(list(gap_obs),10000)))
# kde_obs_func = make_kde_funtion(grid=grid_frb_pp, draws=gap_obs, min_bandwidth=1, max_bandwidth=9, bandwidth_stepsize=1, cv=5, kernel='gaussian')

# plt.plot(grid_frb_pp,kde_obs_func,color='crimson',label=r'DM$_\mathrm{gap,obs}$')
# plt.hist(gap_obs,density=True, bins=500, histtype='stepfilled', alpha=.2, color='crimson')
# plt.axvline(x=np.mean(frb_gaps)-np.mean(psr_gaps),color='k',linestyle='--',linewidth=1,label='Optimal gap')
# plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=12)
# plt.ylabel('PDF',fontsize=12)
# plt.xlim(-80,150)
# plt.legend(fontsize=12,loc=1)
# plt.tight_layout()
# plt.savefig('kde_output_figs/DM_gap_obs.png', dpi=300)
# plt.show()

# print('Optimal DM_gap is',np.mean(frb_gaps)-np.mean(psr_gaps))
# print('Mean DM_gap is',np.mean(gap_obs))

# p = 95. #percentile
# p_lower = (100-p)/2
# p_upper = p+(100-p)/2

# gap_lower = find_percentile(grid_frb_pp, kde_obs_func, p_lower) # min/max from percemtile
# gap_upper = find_percentile(grid_frb_pp, kde_obs_func, p_upper)

# print('Observed lower bound:', gap_lower)
# print('Observed upper bound:', gap_upper)
# print('Std is:', np.std(gap_obs))

# matplotlib.rc('xtick', labelsize=12) 
# matplotlib.rc('ytick', labelsize=12) 
# plt.figure(figsize=(12,8))
# plt.hist(gap_obs, density=True, bins=200, histtype='stepfilled', alpha=.2, color='k')
# plt.plot(grid_frb_pp,kde_gap_opt,color='k',label=r'DM$_\mathrm{gap,obs}$')
# plt.axvline(x=kde_gap_opt,color='red')
# plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
# plt.ylabel('PDF',fontsize=22)
# plt.xlim(-140,180)
# plt.legend(fontsize=22,loc=1)
# plt.tight_layout()
# plt.savefig('kde_output_figs/DM_gap_obs.png', dpi=300)
# plt.show()