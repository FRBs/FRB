# PACKAGE IMPORT
import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import random
import matplotlib
from matplotlib import rcParams
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from DM_definitions import DM_known_draws, FRB_distribution_from_SFR, make_pdf, make_kde_funtion
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')
rcParams['figure.figsize'] = 22,11
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 
"DM_FRB AND SAMPLES"
num_kde_samples = 1000
DM_frb = np.load('DM_outputs/DM_frb.npy')
DM_grid = np.arange(len(DM_frb))
DM_stepsize = 0.1
num_frb_draws = 1000
DM_rand_draws = np.asarray(random.sample(list(DM_frb),num_frb_draws))

"Make KDE distribution kernel with KNOWN DM values, then take loads of samples of this this distribution to make a PDF"
DM_kde_func_ = make_kde_funtion(grid=DM_grid, draws = DM_known_draws, min_bandwidth=25, max_bandwidth=100, bandwidth_stepsize=1, cv=30, kernel='gaussian')
DM_kde_func = DM_kde_func_/np.sum(DM_kde_func_) #normalised
DM_kde = make_pdf(distribution=DM_kde_func,num_of_draws=num_kde_samples,grid=DM_grid,stepsize=DM_stepsize)
DM_kde_draws = DM_kde[0]
DM_kde_distribution  = DM_kde[1] #DM_kde_func scaled

"Make KDE distribution kernel with RANDOM DM values, then take loads of samples of this this distribution to make a PDF"
DM_kde_rand_func_ = make_kde_funtion(grid=DM_grid, draws = DM_rand_draws, min_bandwidth=20, max_bandwidth=60, bandwidth_stepsize=1, cv=30, kernel='gaussian')
DM_kde_rand_func = DM_kde_rand_func_/np.sum(DM_kde_rand_func_) #normalised
DM_kde_rand = make_pdf(distribution=DM_kde_rand_func,num_of_draws=num_kde_samples,grid=DM_grid,stepsize=DM_stepsize)
DM_kde_draws_rand = DM_kde_rand[0]
DM_kde_rand_distribution  = DM_kde_rand[1]

# np.save('DM_outputs/Final_kernel/Final_kernel_choice_2000.npy', DM_kde_rand_distribution/np.sum(DM_kde_rand_distribution))
# np.save('DM_outputs/Final_kernel/DM_grid.npy', DM_grid)
# np.save('DM_outputs/Final_kernel/DM_kde_rand_draws.npy',DM_kde_draws_rand)



"PLOTS"

"Plot for multiple kernels"
# plt.hist(DM_frb, density=True, bins=1000, histtype='stepfilled', alpha=0.5,label=r'$DM_{FRB}$')
# for i, kernel in enumerate(['tophat', 'linear','epanechnikov']):
#     kde_pdfs = make_kde_funtion(grid=DM_grid, draws=DM_known_draws, bandwidth=80, kernel=kernel)
#     plt.plot(DM_grid, kde_pdfs, linewidth=1.5, alpha=0.8, label='%s' % kernel)
#     plt.hist(DM_known_draws, density=True, bins=100, histtype='stepfilled', alpha=0.2)
#     plt.legend(fontsize=14)
#     plt.xlim(0,2600)
# plt.show()
# plt.hist(DM_frb, density=True, bins=1000, histtype='stepfilled', alpha=0.5,label=r'$DM_{FRB}$')
# for i, kernel in enumerate(['gaussian','exponential', 'cosine']):
#     kde_pdfs = make_kde_funtion(grid=DM_grid, draws=DM_rand_draws, bandwidth=120, kernel=kernel)
#     plt.plot(DM_grid, kde_pdfs, linewidth=1.5, alpha=0.5, label='%s' % kernel)
#     plt.hist(DM_kde_draws_rand, density=True, bins=80, histtype='stepfilled', alpha=0.2)
#     plt.legend(fontsize=14)
#     plt.xlim(0,2000)
# plt.show()

"Plot for different bandwidths"
# bandwidths=np.array([10,20,50,100])
# plt.hist(DM_frb, density=True, bins=100, histtype='stepfilled', alpha=0.5,label=r'$DM_{FRB}$ %s' %np.min(DM_frb))
# for i in bandwidths:
#     kde_pdfs = make_kde_funtion(grid=DM_grid, draws=DM_draws, bandwidth=i, kernel='gaussian')
#     plt.plot(DM_grid, kde_pdfs, linewidth=1.2, alpha=0.5,label='{} {}'.format(i, np.percentile(kde_pdfs[2],0.1)))
#     plt.hist(DM_known_draws, density=True, bins=100, histtype='stepfilled', alpha=0.2)
#     plt.title('Gaussian')
#     plt.legend(fontsize=14)
#     plt.xlim(0,2800)
# plt.show()

"Plot DM_FRB and KDE"
# plt.plot(DM_grid, DM_kde_distribution/np.sum(DM_kde_distribution), linewidth=1, alpha=1, label=r'DM observed KDE', color='purple')
# plt.plot(DM_grid, DM_kde_rand_distribution/np.sum(DM_kde_rand_distribution), linewidth=1, alpha=1, label=r'DM random KDE', color='grey')
# plt.hist(DM_known_draws, density=True, bins=100, histtype='stepfilled', alpha=0.5, label=r'DM observed',color='purple')
# plt.hist(DM_rand_draws, density=True, bins=100, histtype='stepfilled', alpha=0.8, label=r'DM random',color='k')
# plt.hist(DM_frb, density=True, bins=500, histtype='stepfilled', alpha=0.5,label=r'DM FRB')
# # for xc in DM_known_draws:
# #     plt.scatter(xc,0,marker='o',s=2,c='k')
# plt.xlabel('$DM$', fontsize=14)
# plt.ylabel('PDF', fontsize=14)
# plt.legend(fontsize=14)
# plt.xlim(0,2200)
# plt.savefig('DM_outputs/DM_kde'+str(num_frb_draws)+'.png')
# plt.show()


"Plot DM_FRB rand draws"
# plt.plot(DM_grid, DM_kde_rand_distribution/np.sum(DM_kde_rand_distribution), linewidth=1, alpha=1, label=r'DM KDE', color='purple')
# # plt.hist(DM_kde_draws_rand, density=True, bins=100, histtype='stepfilled', alpha=0.5, label=r'DM draws from KDE',color='purple')
# plt.hist(DM_frb, density=True, bins=500, histtype='stepfilled', alpha=0.5,label=r'DM simulated')
# # for xc in DM_rand_draws:
# #     plt.scatter(xc,0,marker='o',s=2,c='k')
# plt.xlabel('$DM$', fontsize=14)
# plt.ylabel('PDF', fontsize=14)
# plt.legend(fontsize=14)
# plt.xlim(0,2500)
# plt.savefig('DM_outputs/DM_kde_choice.png')
# plt.show()

# count_obs = 0
# for i in DM_kde_draws: 
#     if i < 100: 
#         count_obs = count_obs + 1
# print(count_obs)

# count_rand = 0
# for i in DM_kde_draws_rand: 
#     if i < 100: 
#         count_rand = count_rand + 1
# print(count_rand)

"Seaborn plot"
sns.distplot(DM_frb, hist = True, kde = False, rug = False, bins=500,
             color = 'darkblue', 
             kde_kws={'linewidth': 2, 'bw':99},
             rug_kws={'color': 'red'})
# sns.distplot(DM_kde_draws, hist = False, kde = False, rug = True, bins=100,
#              color = 'darkblue', 
#              kde_kws={'linewidth': 2, 'bw':99},
#              rug_kws={'color': 'red'})
plt.axvline(x=100,linewidth=2,linestyle='dashed',label='Simulated gap value',color='k')
# plt.title('{}/1000 draws fall below true gap'.format(count_obs),fontsize=40)
plt.xlabel('$DM_{Simulated}$', fontsize=28)
plt.ylabel('PDF', fontsize=28)
plt.legend(fontsize=28)
plt.xlim(0,1500)
plt.tight_layout()
plt.savefig('DM_outputs/KDE_FRB.png')
plt.show()

"Seaborn plot"
# sns.distplot(DM_rand_draws, hist = True, kde = False, rug = False, bins=100,
#              color = 'darkblue', 
#              kde_kws={'linewidth': 2, 'bw':25},
#              rug_kws={'color': 'black'})
# sns.distplot(DM_kde_draws_rand, hist = False, kde = False, rug = True, bins=100,
#              color = 'darkblue', 
#              kde_kws={'linewidth': 2, 'bw':99},
#              rug_kws={'color': 'red'})
# plt.axvline(x=100,linewidth=2,linestyle='dashed',label='Simulated gap value',color='k')
# plt.title('{}/1000 draws fall below simulated gap'.format(count_rand),fontsize=40)
# plt.xlabel('$DM_{Simulated}$', fontsize=28)
# plt.ylabel('PDF', fontsize=28)
# plt.legend(fontsize=28)
# plt.xlim(0,1500)
# # plt.ylim(0,0.01)
# plt.tight_layout()
# plt.savefig('DM_outputs/KDE_rand_counts.png')
# plt.show()