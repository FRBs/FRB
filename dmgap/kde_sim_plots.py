import numpy as np
import pandas as pd
from pandas import DataFrame as df
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 

# Data upload 
frbcat_df = pd.read_csv('transient_data/frbcat_df.csv')
grid_frb = np.load('obs_output_data/grid_frb.npy')

dm_frb_sim = np.load('sim_output_data/dm_frb_sim.npy')
dm_grid = np.load('sim_output_data/dm_grid.npy')

# KDE
#100
kde_100 = np.load('kde_output_data/kde_100.npy')
kde_100_opt = np.load('kde_output_data/kde_100_opt.npy')
kde_100_std = np.load('kde_output_data/kde_100_std.npy')
#1000
kde_1000 = np.load('kde_output_data/kde_1000.npy')
kde_1000_opt = np.load('kde_output_data/kde_1000_opt.npy')
kde_1000_std = np.load('kde_output_data/kde_1000_std.npy')
#5000
kde_5000 = np.load('kde_output_data/kde_5000.npy')
kde_5000_opt = np.load('kde_output_data/kde_5000_opt.npy')
kde_5000_std = np.load('kde_output_data/kde_5000_std.npy')
#10000
kde_10000 = np.load('kde_output_data/kde_10000.npy')
kde_10000_opt = np.load('kde_output_data/kde_10000_opt.npy')
kde_10000_std = np.load('kde_output_data/kde_10000_std.npy')

# DEFT
#100
deft_100_opt = np.load('sim_output_data/optimal_100.npy')
deft_100_boot = np.load('sim_output_data/sampled_100.npy')
deft_std_100 = np.load('sim_output_data/std_100.npy')
#1000
deft_1000_opt = np.load('sim_output_data/optimal_1000.npy')
deft_1000_boot = np.load('sim_output_data/sampled_1000.npy')
deft_std_1000 = np.load('sim_output_data/std_1000.npy')
#5000
deft_5000_opt = np.load('sim_output_data/optimal_5000.npy')
deft_5000_boot = np.load('sim_output_data/sampled_5000.npy')
deft_std_5000 = np.load('sim_output_data/std_5000.npy')
#10000
deft_10000_opt = np.load('sim_output_data/optimal_10000.npy')
deft_10000_boot = np.load('sim_output_data/sampled_10000.npy')
deft_std_10000 = np.load('sim_output_data/std_10000.npy')

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 

plt.figure(figsize=(12,8))
plt.hist(dm_frb_sim,density=True,bins=1000,color='k',alpha=.2)
for i in range(len(kde_100)):
    plt.plot(dm_grid,kde_100[i],color='blue',linewidth=1,alpha=.2)
plt.plot(dm_grid,kde_100_opt,color='k',linewidth=1.5,label='100')

# for i in range(len(kde_1000)):
#     plt.plot(dm_grid,kde_1000[i],color='teal',linewidth=1,alpha=.4)
# plt.plot(dm_grid,kde_1000_opt,color='teal',linewidth=1.5,label='1000')

# for i in range(len(kde_5000)):
#     plt.plot(dm_grid,kde_5000[i],color='green',linewidth=1,alpha=.4)
# plt.plot(dm_grid,kde_5000_opt,color='green',linewidth=1.5,label='5000')

# for i in range(len(kde_10000)):
#     plt.plot(dm_grid,kde_10000[i],color='crimson',linewidth=1,alpha=.4)
# plt.plot(dm_grid,kde_10000_opt,color='crimson',linewidth=1.5,label='10000')

plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
plt.ylabel('PDF',fontsize=22)
plt.xlim(0,2000)
plt.legend(fontsize=22,loc=1)
plt.tight_layout()
# plt.savefig('kde_output_figs/dm_frb_sim.png', dpi=300)
plt.show()


# Plots
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex='col', sharey='row', figsize=(8,12))
ax1.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax1.plot(dm_grid, kde_100_opt, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=100$)')
ax1.fill_between(dm_grid, kde_100_opt-kde_100_std, kde_100_opt+kde_100_std, alpha=.3,color='teal')
ax1.plot(dm_grid, deft_100_opt, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=100$)')
ax1.fill_between(dm_grid, deft_100_opt-deft_std_100, deft_100_opt+deft_std_100, alpha=.4,color='crimson')
ax1.legend(fontsize=14,loc=1)
ax1.set(ylabel='PDF')
ax1.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax1.xaxis.label.set_size(14)
ax1.yaxis.label.set_size(14)

ax2.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax2.plot(dm_grid, kde_1000_opt, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=1000$)')
ax2.fill_between(dm_grid, kde_1000_opt-kde_1000_std, kde_1000_opt+kde_1000_std, alpha=.3,color='teal')
ax2.plot(dm_grid, deft_1000_opt, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=1000$)')
ax2.fill_between(dm_grid, deft_1000_opt-deft_std_1000, deft_1000_opt+deft_std_1000, alpha=.4,color='crimson')
ax2.set(ylabel='PDF')
ax2.legend(fontsize=14,loc=1)
ax2.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax2.xaxis.label.set_size(14)
ax2.yaxis.label.set_size(14)

ax3.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax3.plot(dm_grid, kde_5000_opt, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=5000$)')
ax3.fill_between(dm_grid, kde_5000_opt-kde_5000_std, kde_5000_opt+kde_5000_std, alpha=.3,color='teal')
ax3.plot(dm_grid, deft_5000_opt, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=5000$)')
ax3.fill_between(dm_grid, deft_5000_opt-deft_std_5000, deft_5000_opt+deft_std_5000, alpha=.4,color='crimson')
ax3.set(ylabel='PDF')
ax3.legend(fontsize=14,loc=1)
ax3.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax3.xaxis.label.set_size(14)
ax3.yaxis.label.set_size(14)

ax4.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
ax4.plot(dm_grid, kde_10000_opt, color='teal',linewidth=1,label=r'DM$_\mathrm{FRB}$ KDE ($n=10000$)')
ax4.fill_between(dm_grid, kde_10000_opt-kde_10000_std, kde_10000_opt+kde_10000_std, alpha=.3,color='teal')
ax4.plot(dm_grid, deft_10000_opt, color='crimson',linewidth=1,label=r'DM$_\mathrm{FRB}$ DEFT ($n=10000$)')
ax4.fill_between(dm_grid, deft_10000_opt-deft_std_10000, deft_10000_opt+deft_std_10000, alpha=.4,color='crimson')
ax4.legend(fontsize=14,loc=1)
ax4.axis(xmin=0,xmax=2000,ymin=0,ymax=0.0025)
ax4.set(xlabel=r'DM (pc cm$^{-3}$)',ylabel='PDF')
ax4.xaxis.label.set_size(14)
ax4.yaxis.label.set_size(14)

plt.tight_layout()
plt.savefig('kde_output_figs/kde_sim_vary_n.png', dpi=300)
plt.show()
