import numpy as np
import pandas as pd
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw
from percentile_defs import find_percentile

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.load('obs_output_data/grid_frb.npy')

psrcat_df = pd.read_csv('transient_data/psrcat_df.csv')
grid_psr = np.load('obs_output_data/grid_psr.npy')


# Constraints to region of interest for max concavity
dm_min_frb = -300
dm_max_frb = 300
dm_min_frb_idx = (np.abs(grid_frb - dm_min_frb)).argmin()
dm_max_frb_idx = (np.abs(grid_frb - dm_max_frb)).argmin()

dm_min_psr = -5
dm_max_psr = 50
dm_min_psr_idx = (np.abs(grid_psr - dm_min_psr)).argmin()
dm_max_psr_idx = (np.abs(grid_psr - dm_max_psr)).argmin()

deltaDM_psr_optimal = np.load('obs_output_data/deltaDM_psr_optimal.npy')
deltaDM_psr_optimal = deltaDM_psr_optimal[dm_min_psr_idx:dm_max_psr_idx]
deltaDM_frb_optimal = np.load('obs_output_data/deltaDM_frb_optimal.npy')
deltaDM_frb_optimal = deltaDM_frb_optimal[dm_min_frb_idx:dm_max_frb_idx] 
deltaDM_psr_ensemble = np.load('obs_output_data/deltaDM_psr_ensemble.npy')
deltaDM_psr_ensemble = deltaDM_psr_ensemble[:,dm_min_psr_idx:dm_max_psr_idx]
deltaDM_frb_ensemble = np.load('obs_output_data/deltaDM_frb_ensemble.npy')
deltaDM_frb_ensemble = deltaDM_frb_ensemble[:,dm_min_frb_idx:dm_max_frb_idx]
grid_psr = grid_psr[dm_min_psr_idx:dm_max_psr_idx] 
grid_frb = grid_frb[dm_min_frb_idx:dm_max_frb_idx]

# Pulsars
psr_pp = np.diff(np.diff(deltaDM_psr_optimal)) #second derive
grid_psr_pp = grid_psr[2:]

psr_gap = grid_psr_pp[psr_pp.argmax()]
print('Optimal DeltaDM_psr is',psr_gap)
psr_pp_ensemble = np.diff(np.diff(deltaDM_psr_ensemble))

psr_gaps = []
for i in range(len(deltaDM_psr_ensemble)):
    psr_gaps_ = grid_psr_pp[psr_pp_ensemble[i].argmax()]
    psr_gaps = np.append(psr_gaps,psr_gaps_)
psr_gaps = np.append(psr_gaps,psr_gap)

# FRBs
frb_pp = np.diff(np.diff(deltaDM_frb_optimal))
grid_frb_pp = grid_frb[2:]

frb_gap = grid_frb_pp[frb_pp.argmax()]
print('Optimal DeltaDM_FRB is',frb_gap)
frb_pp_ensemble = np.diff(np.diff(deltaDM_frb_ensemble))

frb_gaps = []
for i in range(len(deltaDM_psr_ensemble)):
    frb_gaps_ = grid_frb_pp[frb_pp_ensemble[i].argmax()]
    frb_gaps= np.append(frb_gaps,frb_gaps_)
frb_gaps = np.append(frb_gaps,frb_gap)

# Gap
gap_obs = []
for i in range(len(frb_gaps)):
    gap_obs_ = frb_gaps[i]-psr_gaps
    gap_obs = np.append(gap_obs,gap_obs_)

alpha_param = 4
b_box_obs = [min(gap_obs)-100,max(gap_obs)+100]
num_ensembles = 100
density_gap_obs = sw.DensityEstimator(gap_obs, alpha=alpha_param, bounding_box=b_box_obs, num_posterior_samples=num_ensembles)
gap_obs_pdf = density_gap_obs.evaluate(grid_frb_pp)

gap_optimal = frb_gap-psr_gap # optimal gap
gap_std = np.std(gap_obs)

print('Optimal DM_gap is',gap_optimal)

p = 95. #percentile
p_lower = (100-p)/2
p_upper = p+(100-p)/2

gap_lower = find_percentile(grid_frb_pp, gap_obs_pdf, p_lower) # min/max from percemtile
gap_upper = find_percentile(grid_frb_pp, gap_obs_pdf, p_upper)

print('Optimal:',gap_optimal)
print('Observed lower bound:', gap_lower)
print('Observed lower bound:', gap_upper)
print('Std is:', gap_std_100)

# matplotlib.rc('xtick', labelsize=12) 
# matplotlib.rc('ytick', labelsize=12) 
# plt.figure(figsize=(12,8))
# plt.hist(gap_obs, density=True, bins=200, histtype='stepfilled', alpha=.2, color='k')
# plt.plot(grid_frb_pp,gap_obs_pdf,color='k',label=r'DM$_\mathrm{gap,obs}$')
# plt.axvline(x=gap_optimal,linestyle='--',linewidth=2,color='k',label=r'Optimal DM$_\mathrm{gap,obs}$')
# plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
# plt.ylabel('PDF',fontsize=22)
# plt.xlim(-140,180)
# plt.legend(fontsize=22,loc=1)
# plt.tight_layout()
# plt.savefig('obs_output_figs/DM_gap_obs.png', dpi=300)
# plt.show()