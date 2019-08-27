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
from DM_definitions import DM_known_draws, make_pdf, make_kde_funtion
from DM_kde_gen import DM_grid, DM_frb, DM_rand_draws
import suftware as sw
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

matplotlib.rc('xtick', labelsize=30) 
matplotlib.rc('ytick', labelsize=30) 
# rcParams['figure.figsize'] = 15,8

num_draws_from_kde = 1000
DM_stepsize = 0.001

kde_original = np.load('DM_outputs/Final_kernel/Final_kernel_choice_2000.npy')

#make kde from function and take draws
kde_original_pdf = make_pdf(distribution=kde_original,num_of_draws=num_draws_from_kde,grid=DM_grid,stepsize=DM_stepsize)
DM_kde_draws_rand = kde_original_pdf[0]
DM_kde_rand_distribution  = kde_original_pdf[1]

plt.plot(DM_grid,DM_kde_rand_distribution/np.sum(DM_kde_rand_distribution))
sns.distplot(DM_kde_draws_rand, hist = True, kde = False, rug = True, bins=200,
             color = 'darkblue', 
             kde_kws={'linewidth': 2, 'bw':25},
             rug_kws={'color': 'black'})
plt.xlabel('$DM_{Simulated}$', fontsize=18)
plt.ylabel('PDF', fontsize=18)
plt.xlim(0,1500)
# plt.ylim(0,0.01)
plt.tight_layout()
plt.savefig('DM_outputs/DM_kde_counts.png')
plt.show()

count = 0
for i in DM_kde_draws_rand: 
    if i < 100: 
        count = count + 1
print(count)



"Plot DM_FRB rand draws"
# plt.plot(DM_grid, DM_kde_rand_distribution/np.sum(DM_kde_rand_distribution), linewidth=1, alpha=1, label=r'DM KDE', color='purple')
# plt.hist(DM_kde_draws_rand, density=True, bins=1000, histtype='stepfilled', alpha=0.5, label=r'DM draws from KDE',color='purple')
# # plt.hist(DM_frb, density=True, bins=500, histtype='stepfilled', alpha=0.5,label=r'DM simulated')
# # for xc in DM_rand_draws:
# #     plt.scatter(xc,0,marker='o',s=2,c='k')
# plt.xlabel('$DM$', fontsize=30)
# plt.ylabel('PDF', fontsize=30)
# plt.legend(fontsize=14)
# plt.xlim(0,2500)
# plt.savefig('DM_outputs/DM_kde_choice.png')
# plt.show()



"Seaborn plot"
# sns.distplot( DM_known_draws, hist = True, kde = True, rug = True, bins=100,
#              color = 'darkblue', 
#              kde_kws={'linewidth': 2, 'bw':99},
#              rug_kws={'color': 'black'})
# plt.xlabel('$DM_{Observed}$', fontsize=28)
# plt.ylabel('PDF', fontsize=28)
# plt.legend(fontsize=28)
# plt.xlim(0,2500)
# plt.tight_layout()
# plt.savefig('DM_outputs/DM_kde_known_sns.png')
# plt.show()

# "Seaborn plot"
# sns.distplot(DM_rand_draws, hist = True, kde = True, rug = True, bins=100,
#              color = 'darkblue', 
#              kde_kws={'linewidth': 2, 'bw':25},
#              rug_kws={'color': 'black'})
# plt.xlabel('$DM_{Simulated}$', fontsize=18)
# plt.ylabel('PDF', fontsize=18)
# plt.xlim(0,1500)
# # plt.ylim(0,0.01)
# plt.tight_layout()
# plt.savefig('DM_outputs/DM_kde_choice_sns_1000.png')
# plt.show()

# print(kde_original)