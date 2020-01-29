import numpy as np
import pandas as pd
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.load('obs_output_data/grid_frb.npy')

psrcat_df = pd.read_csv('transient_data/psrcat_df.csv')
grid_psr = np.load('obs_output_data/grid_psr.npy')

deltaDM_psr_optimal = np.load('obs_output_data/deltaDM_psr_optimal.npy')
deltaDM_psr_optimal = deltaDM_psr_optimal[480000:600000] #remove weird artifacts at edges
deltaDM_frb_optimal = np.load('obs_output_data/deltaDM_frb_optimal.npy')
deltaDM_frb_optimal = deltaDM_frb_optimal[80000:120000] #remove weird artifacts at edges
deltaDM_psr_ensemble = np.load('obs_output_data/deltaDM_psr_ensemble.npy')
deltaDM_psr_ensemble = deltaDM_psr_ensemble[:,480000:600000]
deltaDM_frb_ensemble = np.load('obs_output_data/deltaDM_frb_ensemble.npy')
deltaDM_frb_ensemble = deltaDM_frb_ensemble[:,80000:120000]
grid_psr = grid_psr[480000:600000] 
grid_frb = grid_frb[80000:120000]

# Pulsars
psr_pp = np.diff(np.diff(deltaDM_psr_optimal)) #second derive
grid_psr_pp = grid_psr[2:]

psr_gap = grid_psr_pp[psr_pp.argmax()]
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

matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
plt.hist(gap_obs, density=True, bins=100, histtype='stepfilled', alpha=.4, color='k')
plt.plot(grid_frb_pp,gap_obs_pdf,color='k',label=r'DM$_\mathrm{gap,obs}$')
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=14)
plt.ylabel('PDF',fontsize=12)
plt.xlim(-200,200)
plt.legend(fontsize=12,loc=1)
plt.tight_layout()
plt.savefig('obs_output_figs/DM_gap_obs.png', dpi=300)
plt.show()