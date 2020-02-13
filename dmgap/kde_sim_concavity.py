import numpy as np
import pandas as pd
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import suftware as sw
from scipy.signal import argrelextrema
from percentile_defs import find_percentile
from pdf_defs import make_kde_funtion

# Data upload 
dm_frb_sim = np.load('sim_output_data/dm_frb_sim.npy')
dm_grid = np.load('sim_output_data/dm_grid.npy')

dm_min = -50
dm_max = 150
dm_min_idx = (np.abs(dm_grid - dm_min)).argmin()
dm_max_idx = (np.abs(dm_grid - dm_max)).argmin()

# Constrain to region of interest for max concavity
kde_100_opt =  np.load('kde_output_data/kde_100_opt.npy')
kde_100_opt = kde_100_opt[dm_min_idx:]
kde_100 = np.load('kde_output_data/kde_100.npy')
kde_100 = kde_100[:,dm_min_idx:]

kde_1000_opt =  np.load('kde_output_data/kde_1000_opt.npy')
kde_1000_opt = kde_1000_opt[dm_min_idx:]
kde_1000 = np.load('kde_output_data/kde_1000.npy')
kde_1000 = kde_1000[:,dm_min_idx:]

kde_5000_opt =  np.load('kde_output_data/kde_5000_opt.npy')
kde_5000_opt = kde_5000_opt[dm_min_idx:]
kde_5000 = np.load('kde_output_data/kde_5000.npy')
kde_5000 = kde_5000[:,dm_min_idx:]

kde_10000_opt =  np.load('kde_output_data/kde_10000_opt.npy')
kde_10000_opt = kde_10000_opt[dm_min_idx:]
kde_10000 = np.load('kde_output_data/kde_10000.npy')
kde_10000 = kde_10000[:,dm_min_idx:]

dm_grid = dm_grid[dm_min_idx:]
dm_grid_pp = dm_grid[2:]

# plt.figure(figsize=(12,8))
# plt.hist(dm_frb_sim,density=True,bins=1000,color='k',alpha=.2)
# for i in range(len(kde_100)):
#     plt.plot(dm_grid,kde_100[i],color='blue',linewidth=1,alpha=.2)
# plt.plot(dm_grid,kde_100_opt,color='k',linewidth=1.5,label='100')
# plt.xlim(-200,200)
# plt.show()

# 100
kde_100_opt_pp = np.diff(np.diff(kde_100_opt))
kde_100_pp = np.diff(np.diff(kde_100))
kde_100_gap_idx = argrelextrema(kde_100_opt_pp,np.greater)[0][0]
kde_100_gap = dm_grid_pp[kde_100_gap_idx]

# plt.figure(figsize=(12,8))
# for i in range(len(kde_100_pp)):
#     plt.plot(dm_grid_pp,kde_100_pp[i],color='blue',linewidth=1,alpha=.2)
# plt.plot(dm_grid_pp,kde_100_opt_pp,color='k',linewidth=1.5,label='100')
# plt.xlim(-200,200)
# plt.ylim(0,10**(-11))
# plt.show()

kde_100_gaps = []
for i in range(len(kde_100_pp)):
    plt.plot(dm_grid_pp,kde_100_pp[i],color='k',alpha=.2)
    kde_100_gaps_idx = argrelextrema(kde_100_pp[i],np.greater)[0][0]
    kde_100_gaps_ = dm_grid_pp[kde_100_gaps_idx]
    kde_100_gaps = np.append(kde_100_gaps,kde_100_gaps_)
kde_100_func = make_kde_funtion(grid=dm_grid_pp, draws = kde_100_gaps, min_bandwidth=1, max_bandwidth=8, bandwidth_stepsize=1, cv=5, kernel='gaussian')

# 1000
kde_1000_opt_pp = np.diff(np.diff(kde_1000_opt))
kde_1000_gap_idx = argrelextrema(kde_1000_opt_pp,np.greater)[0][0]
kde_1000_gap = dm_grid_pp[kde_1000_gap_idx]
kde_1000_pp = np.diff(np.diff(kde_1000))
kde_1000_gaps = []
for i in range(len(kde_1000_pp)):
    kde_1000_gaps_idx = argrelextrema(kde_1000_pp[i],np.greater)[0][0]
    kde_1000_gaps_ = dm_grid_pp[kde_1000_gaps_idx]
    kde_1000_gaps = np.append(kde_1000_gaps,kde_1000_gaps_)
kde_1000_func = make_kde_funtion(grid=dm_grid_pp, draws = kde_1000_gaps, min_bandwidth=1, max_bandwidth=8, bandwidth_stepsize=1, cv=5, kernel='gaussian')

# 5000
kde_5000_opt_pp = np.diff(np.diff(kde_5000_opt))
kde_5000_gap_idx = argrelextrema(kde_5000_opt_pp,np.greater)[0][0]
kde_5000_gap = dm_grid_pp[kde_5000_gap_idx]
kde_5000_pp = np.diff(np.diff(kde_5000))
kde_5000_gaps = []
for i in range(len(kde_5000_pp)):
    kde_5000_gaps_idx = argrelextrema(kde_5000_pp[i],np.greater)[0][0]
    kde_5000_gaps_ = dm_grid[kde_5000_gaps_idx]
    kde_5000_gaps = np.append(kde_5000_gaps,kde_5000_gaps_)
kde_5000_func = make_kde_funtion(grid=dm_grid_pp, draws = kde_5000_gaps, min_bandwidth=1, max_bandwidth=5, bandwidth_stepsize=1, cv=5, kernel='gaussian')

# plt.hist(kde_100_gaps,bins=50,alpha=.2,color='blue')
# plt.plot(dm_grid_pp,kde_100_func,color='blue')
# plt.hist(kde_1000_gaps,bins=50,alpha=.2,color='teal')
# plt.plot(dm_grid_pp,kde_1000_func,color='teal')
# plt.hist(kde_5000_gaps,bins=50,alpha=.2,color='crimson')
# plt.plot(dm_grid_pp,kde_5000_func,color='crimson')
# # plt.plot(dm_grid_pp,kde_100_func)
# plt.xlim(-100,200)
# plt.show()

# 10000
kde_10000_opt_pp = np.diff(np.diff(kde_10000_opt))
kde_10000_gap_idx = argrelextrema(kde_10000_opt_pp,np.greater)[0][0]
kde_10000_gap = dm_grid_pp[kde_10000_gap_idx]
kde_10000_pp = np.diff(np.diff(kde_10000))
kde_10000_gaps = []
for i in range(len(kde_10000_pp)):
    kde_10000_gaps_idx = argrelextrema(kde_10000_pp[i],np.greater)[0][0]
    kde_10000_gaps_ = dm_grid[kde_10000_gaps_idx]
    kde_10000_gaps = np.append(kde_10000_gaps,kde_10000_gaps_)
kde_10000_func = make_kde_funtion(grid=dm_grid_pp, draws = kde_10000_gaps, min_bandwidth=0.1, max_bandwidth=2, bandwidth_stepsize=0.1, cv=5, kernel='gaussian')

# plt.figure(figsize=(12,8))
# plt.hist(dm_frb_sim,density=True,bins=1000,color='k',alpha=.2)
# for i in range(50):
#     plt.plot(dm_grid,kde_10000[i],color='blue',linewidth=1,alpha=.2)
# plt.plot(dm_grid,kde_10000_opt,color='k',linewidth=1.5,label='100')
# plt.xlim(-200,1000)
# plt.show()

# Upper and lower gap constraints
p = 95. #percentile
p_lower = (100-p)/2
p_upper = p+(100-p)/2

gap_100_lower = find_percentile(dm_grid_pp, kde_100_func, p_lower) 
gap_100_upper = find_percentile(dm_grid_pp, kde_100_func, p_upper)
gap_1000_lower = find_percentile(dm_grid_pp, kde_1000_func, p_lower) 
gap_1000_upper = find_percentile(dm_grid_pp, kde_1000_func, p_upper)
gap_5000_lower = find_percentile(dm_grid_pp, kde_5000_func, p_lower) 
gap_5000_upper = find_percentile(dm_grid_pp, kde_5000_func, p_upper)
gap_10000_lower = find_percentile(dm_grid_pp, kde_10000_func, p_lower)
gap_10000_upper = find_percentile(dm_grid_pp, kde_10000_func, p_upper)

print('____KDE____')
print('100 Optimal:',kde_100_gap)
print('100 mean:',np.mean(kde_100_gaps))
print('100 lower bound:', gap_100_lower) 
print('100 upper bound:', gap_100_upper) 
print('2 sigma:', 2*np.std(kde_100_gaps))
print('_________')
print('1000 Optimal:',kde_1000_gap)
print('1000 mean:',np.mean(kde_1000_gaps))
print('1000 lower bound:', gap_1000_lower) 
print('1000 upper bound:', gap_1000_upper) 
print('2 sigma:', 2*np.std(kde_1000_gaps))
print('_________')
print('5000 Optimal:',kde_5000_gap)
print('5000 mean:',np.mean(kde_5000_gaps))
print('5000 lower bound:', gap_5000_lower) 
print('5000 upper bound:', gap_5000_upper) 
print('2 sigma:', 2*np.std(kde_5000_gaps))
print('_________')
print('10000 Optimal:',kde_10000_gap)
print('10000 mean:',np.mean(kde_10000_gaps))
print('10000 lower bound:', gap_10000_lower) 
print('10000 upper bound:', gap_10000_upper) 
print('2 sigma:', 2*np.std(kde_10000_gaps))
print('_________')

# for i in range(len(kde_100)):
#     plt.plot(dm_grid_pp,kde_100_pp[i],color='k',alpha=.2)
# plt.plot(dm_grid_pp,kde_100_opt_pp,color='black',linewidth=1.5)
# for i in range(len(kde_1000)):
#     plt.plot(dm_grid_pp,kde_1000_pp[i],color='teal',alpha=.2)
# plt.plot(dm_grid_pp,kde_1000_opt_pp,color='teal',linewidth=1.5)
# for i in range(len(kde_1000)):
#     plt.plot(dm_grid_pp,kde_10000_pp[i],color='crimson',alpha=.2)
# plt.plot(dm_grid_pp,kde_10000_opt_pp,color='crimson',linewidth=1.5)
# plt.xlim(-100,100)
# plt.ylim(0,10**(-11))
# plt.show()

# for i in range(len(kde_100)):
#     plt.plot(dm_grid,kde_100[i],color='k',alpha=.2)
# plt.plot(dm_grid,kde_100_opt,color='k')
# plt.xlim(-100,300)
# plt.show()

# for i in range(len(kde_1000)):
#     plt.plot(dm_grid,kde_1000[i],color='k',alpha=.2)
# plt.plot(dm_grid,kde_1000_opt,color='k')
# plt.xlim(-100,300)
# plt.show()

# for i in range(len(kde_5000)):
#     plt.plot(dm_grid,kde_5000[i],color='k',alpha=.2)
# plt.plot(dm_grid,kde_5000_opt,color='k')
# plt.xlim(-100,300)
# plt.show()

# for i in range(len(kde_10000)):
#     plt.plot(dm_grid,kde_10000[i],color='k',alpha=.1)
# plt.hist(dm_frb_sim,density=True,bins=500,color='k',alpha=.2)
# plt.plot(dm_grid,kde_10000_opt,color='k')
# plt.xlim(-100,300)
# plt.show()

plt.hist(kde_100_gaps,density=True,bins=200,alpha=.2,color='blue')
plt.plot(dm_grid_pp,kde_100_func,color='blue',label=r'$n=100$')
plt.hist(kde_1000_gaps,density=True,bins=100,alpha=.2,color='teal')
plt.plot(dm_grid_pp,kde_1000_func,color='teal',label=r'$n=1000$')
plt.hist(kde_5000_gaps,density=True,bins=100,alpha=.2,color='green')
plt.plot(dm_grid_pp,kde_5000_func,color='green',label=r'$n=5000$')
plt.hist(kde_10000_gaps,density=True,bins=50,alpha=.3,color='crimson')
plt.plot(dm_grid_pp,kde_10000_func,color='crimson',label=r'$n=10000$')
plt.axvline(x=min(dm_frb_sim),color='k',linestyle='--',linewidth=1,label=r'$\Delta$DM$_\mathrm{FRB}$')
plt.xlim(-50,150)
# plt.ylim(0,0.2)
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=12)
plt.ylabel('PDF',fontsize=12)
plt.legend(loc=1)
plt.savefig('kde_output_figs/kde_sim_gaps.png', dpi=300)
plt.show()

# plt.plot(dm_grid,kde_100_opt,label='KDE 100',linestyle='-.')
# # plt.plot(dm_grid,optimal_100,label='DEFT 100',color='teal')
# plt.plot(dm_grid,kde_1000_opt,label='KDE 1000',linestyle='--')
# # plt.plot(dm_grid,optimal_1000,label='DEFT 1000',color='teal',linestyle='--')
# plt.plot(dm_grid,kde_5000_opt,label='KDE 5000',linestyle=':')
# # plt.plot(dm_grid,optimal_5000,label='DEFT 5000',color='teal',linestyle=':')
# plt.plot(dm_grid,kde_10000_opt,label='KDE 10000',linewidth=1)
# for i in range(len(kde_10000)):
#     plt.plot(dm_grid,kde_10000[i],color='k',alpha=.2)
# # plt.plot(dm_grid,optimal_10000,label='DEFT 5000',color='teal',linewidth=1)
# plt.hist(dm_frb_sim,density=True,bins=500,color='k',alpha=.2)
# plt.xlim(-100,1000)
# plt.legend(loc=1)
# plt.show()