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
deft_optimal_2000 = np.load('kde_and_deft_data/deft_2000_optimal.npy')
deft_optimal_4000 = np.load('kde_and_deft_data/deft_4000_optimal.npy')
deft_optimal_10000 = np.load('kde_and_deft_data/deft_10000_optimal.npy')

deft_sampled_frb_len = np.load('kde_and_deft_data/deft_frb_len_sampled.npy')
deft_sampled_100 = np.load('kde_and_deft_data/deft_100_sampled.npy')
deft_sampled_100 = deft_sampled_100/sum(deft_sampled_100)
deft_sampled_1000 = np.load('kde_and_deft_data/deft_1000_sampled.npy')
deft_sampled_1000 = deft_sampled_1000/sum(deft_sampled_1000)
deft_sampled_2000 = np.load('kde_and_deft_data/deft_2000_sampled.npy')
deft_sampled_4000 = np.load('kde_and_deft_data/deft_4000_sampled.npy')
deft_sampled_10000 = np.load('kde_and_deft_data/deft_10000_sampled.npy')

# DEFT with 1000 draws
# matplotlib.rc('xtick', labelsize=10) 
# matplotlib.rc('ytick', labelsize=10) 
# plt.hist(dm_frb_sim, density=True, bins=1000, histtype='stepfilled', alpha=.2, color='black', label=r'Simulated DM$_\mathrm{FRB}$')
# plt.plot(dm_grid,deft_optimal_1000,color='crimson',label=r'DEFT DM$_\mathrm{FRB}$ (n=1000)')
# plt.plot(dm_grid,deft_sampled_1000,color='crimson',linewidth=0.5,alpha=.1)
# plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
# plt.ylabel('PDF',fontsize=14)
# plt.xlim(0,1000)
# plt.legend(fontsize=14)
# plt.tight_layout()
# plt.savefig('sim_outputs/DM_1000_samples.png', dpi=300)
# plt.show()

# # DEFT with 100 draws
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

# 1000 Samples
frb_sim_1000_95 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_1000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_1000[0:,i],quantile=0.05)
    frb_sim_1000_95 = np.append(frb_sim_1000_95,frb_sim_1000_)
# print(frb_sim_1000_95)
frb_95_density_1000 = sw.DensityEstimator(frb_sim_1000_95, bounding_box=[100,150], alpha=3, num_posterior_samples=100)
frb_95_optimal_1000 = frb_95_density_1000.evaluate(dm_grid)

frb_sim_1000_98 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_1000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_1000[0:,i],quantile=0.02)
    frb_sim_1000_98 = np.append(frb_sim_1000_98,frb_sim_1000_)
# print(frb_sim_1000_98)
frb_98_density_1000 = sw.DensityEstimator(frb_sim_1000_98, bounding_box=[90,150], alpha=3, num_posterior_samples=100)
frb_98_optimal_1000 = frb_98_density_1000.evaluate(dm_grid)

frb_sim_1000_99 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_1000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_1000[0:,i],quantile=0.01)
    frb_sim_1000_99 = np.append(frb_sim_1000_99,frb_sim_1000_)
# print(frb_sim_1000_99)
frb_99_density_1000 = sw.DensityEstimator(frb_sim_1000_99, bounding_box=[50,150], alpha=3, num_posterior_samples=100)
frb_99_optimal_1000 = frb_99_density_1000.evaluate(dm_grid)

frb_sim_1000_995 = []
for i in range(len(deft_sampled_1000[1])):
    frb_sim_1000_ = find_quantile(grid=dm_grid,distribution=deft_sampled_1000[0:,i],quantile=0.005)
    frb_sim_1000_995 = np.append(frb_sim_1000_995,frb_sim_1000_)
# print(frb_sim_1000_995)
frb_995_density_1000 = sw.DensityEstimator(frb_sim_1000_995, bounding_box=[30,100], alpha=3, num_posterior_samples=100)
frb_995_optimal_1000 = frb_995_density_1000.evaluate(dm_grid)

# 100 Samples
frb_sim_100_90 = []
for i in range(len(deft_sampled_100[1])):
    frb_sim_100_ = find_quantile(grid=dm_grid,distribution=deft_sampled_100[0:,i],quantile=0.1)
    frb_sim_100_90 = np.append(frb_sim_100_90,frb_sim_100_)
# print(frb_sim_100_90)
frb_90_density_100 = sw.DensityEstimator(frb_sim_100_90, bounding_box=[0,200], alpha=4, num_posterior_samples=100)
frb_90_optimal_100 = frb_90_density_100.evaluate(dm_grid)

frb_sim_100_95 = []
for i in range(len(deft_sampled_100[1])):
    frb_sim_100_ = find_quantile(grid=dm_grid,distribution=deft_sampled_100[0:,i],quantile=0.05)
    frb_sim_100_95 = np.append(frb_sim_100_95,frb_sim_100_)
# print(frb_sim_100_95)
frb_95_density_100 = sw.DensityEstimator(frb_sim_100_95, bounding_box=[0,200], alpha=4, num_posterior_samples=100)
frb_95_optimal_100 = frb_95_density_100.evaluate(dm_grid)

frb_sim_100_98 = []
for i in range(len(deft_sampled_100[1])):
    frb_sim_100_ = find_quantile(grid=dm_grid,distribution=deft_sampled_100[0:,i],quantile=0.02)
    frb_sim_100_98 = np.append(frb_sim_100_98,frb_sim_100_)
# print(frb_sim_100_98)
frb_98_density_100 = sw.DensityEstimator(frb_sim_100_98, bounding_box=[0,150], alpha=4, num_posterior_samples=100)
frb_98_optimal_100 = frb_98_density_100.evaluate(dm_grid)

frb_sim_100_99 = []
for i in range(len(deft_sampled_100[1])):
    frb_sim_100_ = find_quantile(grid=dm_grid,distribution=deft_sampled_100[0:,i],quantile=0.01)
    frb_sim_100_99 = np.append(frb_sim_100_99,frb_sim_100_)
# print(frb_sim_100_99)
frb_99_density_100 = sw.DensityEstimator(frb_sim_100_99, bounding_box=[0,150], alpha=4, num_posterior_samples=100)
frb_99_optimal_100 = frb_99_density_100.evaluate(dm_grid)

# # FRB num samples
# frb_sim_frb_len_90 = []
# for i in range(len(deft_sampled_frb_len[1])):
#     frb_sim_frb_len_ = find_quantile(grid=dm_grid,distribution=deft_sampled_frb_len[0:,i],quantile=0.1)
#     frb_sim_frb_len_90 = np.append(frb_sim_frb_len_90,frb_sim_frb_len_)
# print(frb_sim_frb_len_90)
# frb_90_density_frb_len = sw.DensityEstimator(frb_sim_frb_len_90, bounding_box=[0,150], alpha=4, num_posterior_samples=len(deft_sampled_frb_len[1]))
# frb_90_optimal_frb_len = frb_90_density_frb_len.evaluate(dm_grid)


matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
plt.plot(dm_grid,frb_95_optimal_1000,color='darkgreen',label=r'FRB DM CI=95%',linewidth=1)
plt.plot(dm_grid,frb_98_optimal_1000,color='blue',label=r'FRB DM CI=98%',linewidth=1)
plt.plot(dm_grid,frb_99_optimal_1000,color='red',label=r'FRB DM CI=99%',linewidth=1)
plt.plot(dm_grid,frb_995_optimal_1000,color='orange',label=r'FRB DM CI=99.5%',linewidth=1)
plt.hist(frb_sim_1000_95,density=True,histtype='stepfilled',bins=20,alpha=.2,color='darkgreen')
plt.hist(frb_sim_1000_98,density=True,histtype='stepfilled',bins=20,alpha=.2,color='blue')
plt.hist(frb_sim_1000_99,density=True,histtype='stepfilled',bins=20,alpha=.2,color='red')
plt.hist(frb_sim_1000_995,density=True,histtype='stepfilled',bins=20,alpha=.2,color='orange')
plt.axvline(x=90,linestyle='--',color='k',label='Simulated FRB cut-off')
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=14)
plt.ylabel('PDF',fontsize=14)
plt.title('1000 Draws')
plt.legend(fontsize=14)
plt.xlim([75,145])
plt.tight_layout()
plt.savefig('sim_outputs/dm_gaps_1000.png', dpi=300)
plt.show()

matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
plt.plot(dm_grid,frb_90_optimal_100,color='darkgreen',label=r'FRB DM CI=90%',linewidth=1)
plt.plot(dm_grid,frb_95_optimal_100,color='blue',label=r'FRB DM CI=95%',linewidth=1)
plt.plot(dm_grid,frb_98_optimal_100,color='red',label=r'FRB DM CI=98%',linewidth=1)
plt.plot(dm_grid,frb_99_optimal_100,color='orange',label=r'FRB DM CI=99%',linewidth=1)
plt.hist(frb_sim_100_90,density=True,histtype='stepfilled',bins=20,alpha=.2,color='darkgreen')
plt.hist(frb_sim_100_95,density=True,histtype='stepfilled',bins=20,alpha=.2,color='blue')
plt.hist(frb_sim_100_98,density=True,histtype='stepfilled',bins=20,alpha=.2,color='red')
plt.hist(frb_sim_100_99,density=True,histtype='stepfilled',bins=20,alpha=.2,color='orange')
plt.axvline(x=90,linestyle='--',color='k',label='Simulated FRB cut-off')
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=14)
plt.ylabel('PDF',fontsize=14)
plt.title('100 Draws')
plt.legend(fontsize=14)
plt.xlim([0,200])
plt.tight_layout()
plt.savefig('sim_outputs/dm_gaps_100.png', dpi=300)
plt.show()