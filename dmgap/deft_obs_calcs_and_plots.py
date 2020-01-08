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

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.arange(0,6000,0.01)

psrcat_df = pd.read_csv('transient_data/psrcat_df.csv')
grid_psr = np.arange(-500,500,0.01)

alpha_param = 3 #DEFT smoothing parameter

fig, ax1 = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(12,8))

# Pulsars
dmdiff_psr_density = sw.DensityEstimator(np.array(psrcat_df['dmdiff']), bounding_box=[-200,200], alpha=alpha_param, num_posterior_samples=len(psrcat_df['dmdiff']))
dmdiff_psr_optimal_nonnorm = dmdiff_psr_density.evaluate(grid_psr)
dmdiff_psr_optimal = dmdiff_psr_optimal_nonnorm/np.sum(dmdiff_psr_optimal_nonnorm) #normalise
dmdiff_psr_ensemble_nonnorm = dmdiff_psr_density.evaluate_samples(grid_psr)
# ax1.hist(np.array(psrcat_df['dmdiff']),density=True,histtype='stepfilled',bins=100,alpha=.2,label=r'Pulsar Histogram',color='black')
dmdiff_psr_ensemble = []
for i in range(len(dmdiff_psr_ensemble_nonnorm[0])):
    dmdiff_psr_ensemble_ = dmdiff_psr_ensemble_nonnorm[:,i]/np.sum(dmdiff_psr_ensemble_nonnorm[:,i])
    dmdiff_psr_ensemble.append(dmdiff_psr_ensemble_)
    # ax1.plot(grid_psr,dmdiff_psr_ensemble_*100,color='#f46036',linewidth=1,alpha=.2)
# ax1.plot(grid_psr,dmdiff_psr_optimal*100,color='#f46036',label=r'$\Delta$ DM$_\mathrm{pulsar}$')
# ax1.axis(xmin=-80,xmax=70,ymin=0)
# ax1.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
# ax1.xaxis.label.set_size(22)
# ax1.yaxis.label.set_size(22)
# ax1.legend(fontsize=22)
# plt.tight_layout()
# plt.savefig('obs_outputs/dm_psr.png', dpi=300)
# plt.show()

print('Pulsar calculations done.')

# FRBs
fig, ax2 = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(12,8))

dmdiff_frb_density = sw.DensityEstimator(np.array(frbcat_df['dmdiff']), bounding_box=[0,4000], alpha=alpha_param, num_posterior_samples=len(frbcat_df['dmdiff']))
dmdiff_frb_optimal_nonnorm = dmdiff_frb_density.evaluate(grid_frb)
dmdiff_frb_optimal = dmdiff_frb_optimal_nonnorm/np.sum(dmdiff_frb_optimal_nonnorm) #normalise
dmdiff_frb_ensemble_nonnorm = dmdiff_frb_density.evaluate_samples(grid_frb)
# ax2.hist(np.array(frbcat_df['dmdiff']),density=True,histtype='stepfilled',bins=60,alpha=.2,label=r'FRB Histogram',color='black')
dmdiff_frb_ensemble = []
for i in range(len(dmdiff_frb_ensemble_nonnorm[0])):
    dmdiff_frb_ensemble_ = dmdiff_frb_ensemble_nonnorm[:,i]/np.sum(dmdiff_frb_ensemble_nonnorm[:,i])
    dmdiff_frb_ensemble.append(dmdiff_frb_ensemble_)
    # ax2.plot(grid_frb,dmdiff_frb_ensemble_*100,color='crimson',linewidth=1,alpha=.2)
# ax2.plot(grid_frb,dmdiff_frb_optimal*100,color='crimson',label=r'$\Delta$ DM$_\mathrm{FRB}$')
# ax2.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
# ax2.axis(xmin=0,xmax=2000,ymin=0)
# ax2.xaxis.label.set_size(22)
# ax2.yaxis.label.set_size(22)
# ax2.legend(fontsize=22)
# plt.tight_layout()
# plt.savefig('obs_outputs/dm_frb.png', dpi=300)
# plt.show()
print('FRB calculations done.')

np.save('obs_outputs/dmdiff_psr_optimal.npy',dmdiff_psr_optimal)
np.save('obs_outputs/dmdiff_psr_ensemble.npy',dmdiff_psr_ensemble)
np.save('obs_outputs/dmdiff_frb_optimal.npy',dmdiff_frb_optimal)
np.save('obs_outputs/dmdiff_frb_ensemble.npy',dmdiff_frb_ensemble)







