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

# Data upload 
dm_frb_sim = np.load('sim_outputs/dm_frb_sim.npy')
dm_grid = np.load('sim_outputs/dm_grid.npy')

deft_optimal_frb_len = np.load('kde_and_deft_data/deft_frb_len_optimal.npy')
deft_optimal_100 = np.load('kde_and_deft_data/deft_100_optimal.npy')
deft_optimal_1000 = np.load('kde_and_deft_data/deft_1000_optimal.npy')
deft_optimal_10000 = np.load('kde_and_deft_data/deft_10000_optimal.npy')

deft_sampled_frb_len = np.load('kde_and_deft_data/deft_frb_len_sampled.npy')
deft_sampled_100 = np.load('kde_and_deft_data/deft_100_sampled.npy')
deft_sampled_100 = deft_sampled_100/sum(deft_sampled_100)
deft_sampled_1000 = np.load('kde_and_deft_data/deft_1000_sampled.npy')
deft_sampled_1000 = deft_sampled_1000/sum(deft_sampled_1000)
deft_sampled_10000 = np.load('kde_and_deft_data/deft_10000_sampled.npy')
deft_sampled_10000 = deft_sampled_10000/sum(deft_sampled_10000)

# DEFT with 1000 draws
# matplotlib.rc('xtick', labelsize=10) 
# matplotlib.rc('ytick', labelsize=10) 
# plt.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
# plt.plot(dm_grid,deft_optimal_1000,color='crimson',label=r'DEFT DM$_\mathrm{FRB}$ (n=1000)')
# plt.plot(dm_grid,deft_sampled_1000,color='crimson',linewidth=0.5,alpha=.1)
# plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
# plt.ylabel('PDF',fontsize=14)
# plt.xlim(0,2000)
# plt.legend(fontsize=14)
# plt.tight_layout()
# plt.savefig('sim_outputs/DM_1000_samples.png', dpi=300)
# plt.show()

# # # DEFT with 100 draws
# matplotlib.rc('xtick', labelsize=10) 
# matplotlib.rc('ytick', labelsize=10) 
# plt.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
# plt.plot(dm_grid,deft_optimal_100,color='crimson',label=r'DEFT DM$_\mathrm{FRB}$ (n=100)')
# plt.plot(dm_grid,deft_sampled_100,color='crimson',linewidth=0.5,alpha=.1)
# plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=14)
# plt.ylabel('PDF',fontsize=14)
# plt.xlim(0,2000)
# plt.legend(fontsize=14)
# plt.tight_layout()
# plt.savefig('sim_outputs/DM_100_samples.png', dpi=300)
# plt.show()

# 10000 Samples
frb_sim_10000_995 = []
for i in range(len(deft_sampled_10000[1])):
    frb_sim_10000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_10000[0:,i],quantile=0.005)
    frb_sim_10000_995 = np.append(frb_sim_10000_995,frb_sim_10000_)
frb_995_density_10000 = sw.DensityEstimator(frb_sim_10000_995, bounding_box=[min(frb_sim_10000_995)-20,max(frb_sim_10000_995)+20], alpha=3, num_posterior_samples=100)
frb_995_optimal_10000 = frb_995_density_10000.evaluate(dm_grid)

frb_sim_10000_998 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_10000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_10000[0:,i],quantile=0.002)
    frb_sim_10000_998 = np.append(frb_sim_10000_998,frb_sim_10000_)
frb_998_density_10000 = sw.DensityEstimator(frb_sim_10000_998, bounding_box=[min(frb_sim_10000_998)-20,max(frb_sim_10000_998)+20], alpha=3, num_posterior_samples=100)
frb_998_optimal_10000 = frb_998_density_10000.evaluate(dm_grid)

frb_sim_10000_999 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_10000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_10000[0:,i],quantile=0.001)
    frb_sim_10000_999 = np.append(frb_sim_10000_999,frb_sim_10000_)
frb_999_density_10000 = sw.DensityEstimator(frb_sim_10000_999, bounding_box=[min(frb_sim_10000_999)-20,max(frb_sim_10000_999)+20], alpha=3, num_posterior_samples=100)
frb_999_optimal_10000 = frb_999_density_10000.evaluate(dm_grid)

# 1000 Samples
frb_sim_1000_98 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_1000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_1000[0:,i],quantile=0.02)
    frb_sim_1000_98 = np.append(frb_sim_1000_98,frb_sim_1000_)
# print(frb_sim_1000_98)
frb_98_density_1000 = sw.DensityEstimator(frb_sim_1000_98, bounding_box=[min(frb_sim_1000_98)-20,max(frb_sim_1000_98)+20], alpha=3, num_posterior_samples=100)
frb_98_optimal_1000 = frb_98_density_1000.evaluate(dm_grid)

frb_sim_1000_99 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_1000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_1000[0:,i],quantile=0.01)
    frb_sim_1000_99 = np.append(frb_sim_1000_99,frb_sim_1000_)
# print(frb_sim_1000_99)
frb_99_density_1000 = sw.DensityEstimator(frb_sim_1000_99, bounding_box=[min(frb_sim_1000_99)-20,max(frb_sim_1000_99)+20], alpha=3, num_posterior_samples=100)
frb_99_optimal_1000 = frb_99_density_1000.evaluate(dm_grid)

frb_sim_1000_995 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_1000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_1000[0:,i],quantile=0.005)
    frb_sim_1000_995 = np.append(frb_sim_1000_995,frb_sim_1000_)
# print(frb_sim_1000_995)
frb_995_density_1000 = sw.DensityEstimator(frb_sim_1000_995, bounding_box=[min(frb_sim_1000_995)-20,max(frb_sim_1000_995)+20], alpha=3, num_posterior_samples=100)
frb_995_optimal_1000 = frb_995_density_1000.evaluate(dm_grid)

# 100 Samples
frb_sim_100_97 = []
for i in range(len(deft_sampled_100[1])):
    frb_sim_100_ = find_quantile(grid=dm_grid,distribution=deft_sampled_100[0:,i],quantile=0.03)
    frb_sim_100_97 = np.append(frb_sim_100_97,frb_sim_100_)
# print(frb_sim_100_97)
frb_97_density_100 = sw.DensityEstimator(frb_sim_100_97, bounding_box=[min(frb_sim_100_97)-20,max(frb_sim_100_97)+20], alpha=4, num_posterior_samples=100)
frb_97_optimal_100 = frb_97_density_100.evaluate(dm_grid)

frb_sim_100_98 = []
for i in range(len(deft_sampled_100[1])):
    frb_sim_100_ = find_quantile(grid=dm_grid,distribution=deft_sampled_100[0:,i],quantile=0.02)
    frb_sim_100_98 = np.append(frb_sim_100_98,frb_sim_100_)
# print(frb_sim_100_98)
frb_98_density_100 = sw.DensityEstimator(frb_sim_100_98, bounding_box=[min(frb_sim_100_98)-20,max(frb_sim_100_98)+20], alpha=4, num_posterior_samples=100)
frb_98_optimal_100 = frb_98_density_100.evaluate(dm_grid)

frb_sim_100_99 = []
for i in range(len(deft_sampled_100[1])):
    frb_sim_100_ = find_quantile(grid=dm_grid,distribution=deft_sampled_100[0:,i],quantile=0.01)
    frb_sim_100_99 = np.append(frb_sim_100_99,frb_sim_100_)
# print(frb_sim_100_99)
frb_99_density_100 = sw.DensityEstimator(frb_sim_100_99, bounding_box=[min(frb_sim_100_99)-20,max(frb_sim_100_99)+20], alpha=4, num_posterior_samples=100)
frb_99_optimal_100 = frb_99_density_100.evaluate(dm_grid)

# # FRB num samples
# frb_sim_frb_len_90 = []
# for i in range(len(deft_sampled_frb_len[1])):
#     frb_sim_frb_len_ = find_quantile(grid=dm_grid,distribution=deft_sampled_frb_len[0:,i],quantile=0.1)
#     frb_sim_frb_len_90 = np.append(frb_sim_frb_len_90,frb_sim_frb_len_)
# print(frb_sim_frb_len_90)
# frb_90_density_frb_len = sw.DensityEstimator(frb_sim_frb_len_90, bounding_box=[0,150], alpha=4, num_posterior_samples=len(deft_sampled_frb_len[1]))
# frb_90_optimal_frb_len = frb_90_density_frb_len.evaluate(dm_grid)

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row', figsize=(10,12))

ax1.plot(dm_grid,frb_97_optimal_100,color='blue',label=r'DM$_\mathrm{FRB}^\mathrm{min}$ ($p=97$%)',linewidth=1)
ax1.plot(dm_grid,frb_98_optimal_100,color='red',label=r'DM$_\mathrm{FRB}^\mathrm{min}$ ($p=98$%)',linewidth=1)
ax1.plot(dm_grid,frb_99_optimal_100,color='darkgreen',label=r'DM$_\mathrm{FRB}^\mathrm{min}$ ($p=99$%)',linewidth=1)
ax1.hist(frb_sim_100_97,density=True,histtype='stepfilled',bins=20,alpha=.2,color='blue')
ax1.hist(frb_sim_100_98,density=True,histtype='stepfilled',bins=20,alpha=.2,color='red')
ax1.hist(frb_sim_100_99,density=True,histtype='stepfilled',bins=20,alpha=.2,color='darkgreen')
ax1.axvline(x=90,linestyle='--',color='k',linewidth=1) #,label=r'DM=90 pc cm$^{-3}$'
ax1.legend(fontsize=14,loc=1)
ax1.set(ylabel='PDF',title='100 Samples')
ax1.axis(xmin=0,xmax=250)
ax1.xaxis.label.set_size(14)
ax1.yaxis.label.set_size(14)

ax2.plot(dm_grid,frb_98_optimal_1000,color='blue',label=r'DM$_\mathrm{FRB}^\mathrm{min}$ ($p=98$%)',linewidth=1)
ax2.plot(dm_grid,frb_99_optimal_1000,color='red',label=r'DM$_\mathrm{FRB}^\mathrm{min}$ ($p=99$%)',linewidth=1)
ax2.plot(dm_grid,frb_995_optimal_1000,color='darkgreen',label=r'DM$_\mathrm{FRB}^\mathrm{min}$ ($p=99.5$%)',linewidth=1)
ax2.hist(frb_sim_1000_98,density=True,histtype='stepfilled',bins=20,alpha=.2,color='blue')
ax2.hist(frb_sim_1000_99,density=True,histtype='stepfilled',bins=20,alpha=.2,color='red')
ax2.hist(frb_sim_1000_995,density=True,histtype='stepfilled',bins=20,alpha=.2,color='darkgreen')
ax2.axvline(x=90,linestyle='--',color='k',linewidth=1)
ax2.legend(fontsize=14,loc=1)
ax2.set(ylabel='PDF',title='1000 Samples')
ax2.axis(xmin=0,xmax=250)
ax2.xaxis.label.set_size(14)
ax2.yaxis.label.set_size(14)

ax3.plot(dm_grid,frb_995_optimal_10000,color='blue',label=r'DM$_\mathrm{FRB}^\mathrm{min}$ ($p=99.5$%)',linewidth=1)
ax3.plot(dm_grid,frb_998_optimal_10000,color='red',label=r'DM$_\mathrm{FRB}^\mathrm{min}$ ($p=99.8$%)',linewidth=1)
ax3.plot(dm_grid,frb_999_optimal_10000,color='darkgreen',label=r'DM$_\mathrm{FRB}^\mathrm{min}$ ($p=99.9$%)',linewidth=1)
ax3.hist(frb_sim_10000_995,density=True,histtype='stepfilled',bins=20,alpha=.2,color='blue')
ax3.hist(frb_sim_10000_998,density=True,histtype='stepfilled',bins=20,alpha=.2,color='red')
ax3.hist(frb_sim_10000_999,density=True,histtype='stepfilled',bins=20,alpha=.2,color='darkgreen')
ax3.axvline(x=90,linestyle='--',color='k',linewidth=1)
ax3.legend(fontsize=14,loc=1)
ax3.set(xlabel=r'DM (pc cm$^{-3}$)', ylabel='PDF',title='10000 Samples')
ax3.axis(xmin=0,xmax=250,ymin=0,ymax=0.3)
ax3.xaxis.label.set_size(14)
ax3.yaxis.label.set_size(14)

plt.tight_layout()
plt.savefig('sim_outputs/dm_gaps_vary_n.png', dpi=300)
plt.show()