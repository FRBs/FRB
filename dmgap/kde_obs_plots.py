import numpy as np
import pandas as pd
from pandas import DataFrame as df
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt

# Data upload 
psrcat_df = pd.read_csv('transient_data/psrcat_df.csv')
grid_psr = np.load('obs_output_data/grid_psr.npy')
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.load('obs_output_data/grid_frb.npy')

delta_dm_psr = psrcat_df['deltaDM']
delta_dm_frb = frbcat_df['deltaDM']

kde_psr = np.load('kde_output_data/kde_psr.npy')
kde_psr_opt = np.load('kde_output_data/kde_psr_opt.npy')
kde_frb = np.load('kde_output_data/kde_frb.npy')
kde_frb_opt = np.load('kde_output_data/kde_frb_opt.npy')

# Plots
fig, ax1 = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(12,8))
ax1.hist(delta_dm_psr, density=True, bins=100, histtype='stepfilled', alpha=.2, color='black', label=r'FRB Histogram')
ax1.plot(grid_psr, kde_psr_opt, color='#f46036',linewidth=1.5,label=r'$\Delta$DM$_\mathrm{pulsar}$')
for i in range(100):
    ax1.plot(grid_psr,kde_psr[i],color='#f46036',linewidth=1,alpha=.1)
ax1.legend(fontsize=14,loc=1)
ax1.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
ax1.axis(xmin=-100,xmax=50,ymin=0)
ax1.xaxis.label.set_size(14)
ax1.yaxis.label.set_size(14)
plt.tight_layout()
plt.savefig('kde_output_figs/kde_obs_pulsar.png', dpi=300)
plt.show()

fig, ax2 = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(12,8))
ax2.hist(delta_dm_frb, density=True, bins=100, histtype='stepfilled', alpha=.2, color='black', label=r'Pulsar Histogram')
ax2.plot(grid_frb, kde_frb_opt, color='crimson',linewidth=1.5,label=r'$\Delta$DM$_\mathrm{FRB}$')
for i in range(100):
    ax2.plot(grid_frb,kde_frb[i],color='crimson',linewidth=1,alpha=.1)
ax2.legend(fontsize=14,loc=1)
ax2.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
ax2.axis(xmin=-500,xmax=3000,ymin=0)
ax2.xaxis.label.set_size(14)
ax2.yaxis.label.set_size(14)
plt.tight_layout()
plt.savefig('kde_output_figs/kde_obs_frb.png', dpi=300)
plt.show()
