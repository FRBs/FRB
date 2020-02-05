import numpy as np
import pandas as pd
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.arange(-1000,6000,0.01)
np.save('obs_output_data/grid_frb.npy', grid_frb)

psrcat_df = pd.read_csv('transient_data/psrcat_df.csv')
grid_psr = np.arange(-500,500,0.001)
np.save('obs_output_data/grid_psr.npy', grid_psr)

alpha_param = 3 #DEFT smoothing parameter

fig, ax1 = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(12,8))

# Pulsars
deltaDM_psr_density = sw.DensityEstimator(np.array(psrcat_df['deltaDM']), bounding_box=[-200,200], alpha=alpha_param, num_posterior_samples=100)
deltaDM_psr_optimal_nonnorm = deltaDM_psr_density.evaluate(grid_psr)
deltaDM_psr_optimal = deltaDM_psr_optimal_nonnorm/np.sum(deltaDM_psr_optimal_nonnorm) #normalise
deltaDM_psr_ensemble_nonnorm = deltaDM_psr_density.evaluate_samples(grid_psr)
ax1.hist(np.array(psrcat_df['deltaDM']),density=True,histtype='stepfilled',bins=100,alpha=.2,label=r'Pulsar Histogram',color='black')
deltaDM_psr_ensemble = []
for i in range(len(deltaDM_psr_ensemble_nonnorm[0])):
    deltaDM_psr_ensemble_ = deltaDM_psr_ensemble_nonnorm[:,i]/np.sum(deltaDM_psr_ensemble_nonnorm[:,i])
    deltaDM_psr_ensemble.append(deltaDM_psr_ensemble_)
    ax1.plot(grid_psr,deltaDM_psr_ensemble_*1000,color='#f46036',linewidth=1,alpha=.1)
ax1.plot(grid_psr,deltaDM_psr_optimal*1000,color='#f46036',linewidth=1.5,label=r'$\Delta$DM$_\mathrm{pulsar}$')
ax1.axis(xmin=-100,xmax=50,ymin=0)
ax1.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
ax1.xaxis.label.set_size(22)
ax1.yaxis.label.set_size(22)
ax1.legend(fontsize=22,loc=1)
plt.tight_layout()
plt.savefig('obs_output_figs/dm_psr.png', dpi=300)
plt.show()

print('There are currently',len(psrcat_df), 'pulsars in the sample.')
print('Pulsar calculations done.')

# FRBs
fig, ax2 = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(12,8))

deltaDM_frb_density = sw.DensityEstimator(np.array(frbcat_df['deltaDM']), bounding_box=[-600,5000], alpha=alpha_param, num_posterior_samples=100)
deltaDM_frb_optimal_nonnorm = deltaDM_frb_density.evaluate(grid_frb)
deltaDM_frb_optimal = deltaDM_frb_optimal_nonnorm/np.sum(deltaDM_frb_optimal_nonnorm) #normalise
deltaDM_frb_ensemble_nonnorm = deltaDM_frb_density.evaluate_samples(grid_frb)
ax2.hist(np.array(frbcat_df['deltaDM']),density=True,histtype='stepfilled',bins=60,alpha=.2,label=r'FRB Histogram',color='black')
deltaDM_frb_ensemble = []
for i in range(len(deltaDM_frb_ensemble_nonnorm[0])):
    deltaDM_frb_ensemble_ = deltaDM_frb_ensemble_nonnorm[:,i]/np.sum(deltaDM_frb_ensemble_nonnorm[:,i])
    deltaDM_frb_ensemble.append(deltaDM_frb_ensemble_)
    ax2.plot(grid_frb,deltaDM_frb_ensemble_*100,color='crimson',linewidth=1,alpha=.1)
ax2.plot(grid_frb,deltaDM_frb_optimal*100,color='crimson',linewidth=1.5,label=r'$\Delta$DM$_\mathrm{FRB}$')
ax2.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
ax2.axis(xmin=-500,xmax=3000,ymin=0)
ax2.xaxis.label.set_size(22)
ax2.yaxis.label.set_size(22)
ax2.legend(fontsize=22,loc=1)
plt.tight_layout()
plt.savefig('obs_output_figs/dm_frb.png', dpi=300)
plt.show()

print('There are currently',len(frbcat_df), 'FRBs in the sample.')
print('FRB calculations done.')

# np.save('obs_output_data/deltaDM_psr_optimal.npy',deltaDM_psr_optimal)
# np.save('obs_output_data/deltaDM_psr_ensemble.npy',deltaDM_psr_ensemble)
# np.save('obs_output_data/deltaDM_frb_optimal.npy',deltaDM_frb_optimal)
# np.save('obs_output_data/deltaDM_frb_ensemble.npy',deltaDM_frb_ensemble)

# print('PDF outputs saved.')




