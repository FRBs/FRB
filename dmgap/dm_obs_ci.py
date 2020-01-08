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
import seaborn as sns

# Data upload 
grid_frb = np.arange(0,6000,0.01)
grid_psr = np.arange(-500,500,0.01)

dmdiff_psr_optimal = np.load('obs_outputs/dmdiff_psr_optimal.npy')
dmdiff_frb_optimal = np.load('obs_outputs/dmdiff_frb_optimal.npy')
dmdiff_psr_ensemble = np.load('obs_outputs/dmdiff_psr_ensemble.npy')
dmdiff_frb_ensemble = np.load('obs_outputs/dmdiff_frb_ensemble.npy')

# plt.plot(grid_psr,dmdiff_psr_ensemble[0,0:])
# plt.show()

ci_interval = [.95, .96, .98] #confidence
ci_cut_off_psr = np.sqrt(ci_interval) #for joint probability: probability that <ci_interval% of events will occur simultaneously
ci_cut_off_frb = 1-np.sqrt(ci_interval)

# frb_obs_dm = find_quantile(grid=grid_frb,distribution=dmdiff_frb_optimal,quantile=ci_cut_off_frb)
# psr_optimal_dm = find_quantile(grid=grid_psr,distribution=dmdiff_psr_optimal,quantile=ci_cut_off_psr)

# print(frb_obs_dm)
# print(psr_optimal_dm)

# 95% CI for FRBs and pulsars
frb_ensemble_dm_95 = []
for i in range(len(dmdiff_frb_ensemble)):
    frb_ensemble_dm_ = find_quantile(grid=grid_frb,distribution=dmdiff_frb_ensemble[i,0:],quantile=ci_cut_off_frb[0])
    frb_ensemble_dm_95 = np.append(frb_ensemble_dm_95,frb_ensemble_dm_)

psr_ensemble_dm_95 = []
for i in range(len(dmdiff_psr_ensemble)):
    psr_ensemble_dm_ = find_quantile(grid=grid_psr,distribution=dmdiff_psr_ensemble[i,0:],quantile=ci_cut_off_psr[0])
    psr_ensemble_dm_95 = np.append(psr_ensemble_dm_95,psr_ensemble_dm_)

# 96% CI for FRBs and pulsars
frb_ensemble_dm_96 = []
for i in range(len(dmdiff_frb_ensemble)):
    frb_ensemble_dm_ = find_quantile(grid=grid_frb,distribution=dmdiff_frb_ensemble[i,0:],quantile=ci_cut_off_frb[1])
    frb_ensemble_dm_96 = np.append(frb_ensemble_dm_96,frb_ensemble_dm_)

psr_ensemble_dm_96 = []
for i in range(len(dmdiff_psr_ensemble)):
    psr_ensemble_dm_ = find_quantile(grid=grid_psr,distribution=dmdiff_psr_ensemble[i,0:],quantile=ci_cut_off_psr[1])
    psr_ensemble_dm_96 = np.append(psr_ensemble_dm_96,psr_ensemble_dm_)

# 98% CI for FRBs and pulsars
frb_ensemble_dm_98 = []
for i in range(len(dmdiff_frb_ensemble)):
    frb_ensemble_dm_ = find_quantile(grid=grid_frb,distribution=dmdiff_frb_ensemble[i,0:],quantile=ci_cut_off_frb[2])
    frb_ensemble_dm_98 = np.append(frb_ensemble_dm_98,frb_ensemble_dm_)

psr_ensemble_dm_98 = []
for i in range(len(dmdiff_psr_ensemble)):
    psr_ensemble_dm_ = find_quantile(grid=grid_psr,distribution=dmdiff_psr_ensemble[i,0:],quantile=ci_cut_off_psr[2])
    psr_ensemble_dm_98 = np.append(psr_ensemble_dm_98,psr_ensemble_dm_)

gap_95 = []
for i in range(len(frb_ensemble_dm_95)):
    gap_95_ = frb_ensemble_dm_95[i]-psr_ensemble_dm_95
    gap_95 = np.append(gap_95,gap_95_)
gap_95_density = sw.DensityEstimator(gap_95, bounding_box=[-200,200], alpha=3, num_posterior_samples=100)
gap_95_optimal = gap_95_density.evaluate(grid_psr)

gap_96 = []
for i in range(len(frb_ensemble_dm_96)):
    gap_96_ = frb_ensemble_dm_96[i]-psr_ensemble_dm_96
    gap_96 = np.append(gap_96,gap_96_)
gap_96_density = sw.DensityEstimator(gap_96, bounding_box=[-200,200], alpha=3, num_posterior_samples=100)
gap_96_optimal = gap_96_density.evaluate(grid_psr)

gap_98 = []
for i in range(len(frb_ensemble_dm_98)):
    gap_98_ = frb_ensemble_dm_98[i]-psr_ensemble_dm_98
    gap_98 = np.append(gap_98,gap_98_)
gap_98_density = sw.DensityEstimator(gap_98, bounding_box=[-200,200], alpha=3, num_posterior_samples=100)
gap_98_optimal = gap_98_density.evaluate(grid_psr)

matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
plt.plot(grid_psr,gap_95_optimal,color='darkgreen',label=r'DM Gap CI=95%',linewidth=1)
plt.plot(grid_psr,gap_96_optimal,color='blue',label=r'DM Gap CI=96%',linewidth=1)
plt.plot(grid_psr,gap_98_optimal,color='red',label=r'DM Gap CI=98%',linewidth=1)
plt.hist(gap_95,density=True,histtype='stepfilled',bins=100,alpha=.2,color='darkgreen')
plt.hist(gap_96,density=True,histtype='stepfilled',bins=100,alpha=.2,color='blue')
plt.hist(gap_98,density=True,histtype='stepfilled',bins=100,alpha=.2,color='red')
plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=14)
plt.ylabel('PDF',fontsize=14)
plt.legend(fontsize=14)
plt.xlim([-25,80])
plt.tight_layout()
plt.savefig('obs_outputs/dm_gaps.png', dpi=300)
plt.show()


# plt.hist(np.array(psr_ensemble_dm_95),density=True,histtype='stepfilled',bins=50,alpha=.2,label=r'Pulsars',color='darkgreen')
# plt.hist(np.array(frb_ensemble_dm_95),density=True,histtype='stepfilled',bins=50,alpha=.2,label=r'FRB',color='red')
# plt.xlabel(r'DM (pc cm$^{-3}$)',fontsize=22)
# plt.ylabel('PDF',fontsize=22)
# plt.legend(fontsize=22)
# plt.tight_layout()
# plt.show()

# plt.hist(np.array(frb_ensemble_dm_95-psr_ensemble_dm_95),density=True,histtype='stepfilled',bins=50,alpha=.2,label=r'CI=95%')
# plt.hist(np.array(frb_ensemble_dm_98-psr_ensemble_dm_98),density=True,histtype='stepfilled',bins=50,alpha=.2,label=r'CI=98%')
# plt.hist(np.array(frb_ensemble_dm_99-psr_ensemble_dm_99),density=True,histtype='stepfilled',bins=50,alpha=.2,label=r'CI=99%')
