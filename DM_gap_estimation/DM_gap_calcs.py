import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import sklearn
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import matplotlib
import matplotlib.pyplot as plt
from importlib import reload
from DM_definitions import DM_known_draws, FRB_distribution_from_SFR, make_pdf, make_kde_funtion
import DM_kde_gen
from DM_kde_gen import DM_frb, DM_grid, DM_kde_rand_distribution, DM_kde_distribution, DM_kde_draws_rand, DM_kde_draws
import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

num_runs = 50
percentile_val = 2

"FIND GAP OF THE SIMULATED DATA"
DM_frb_runs = []
for _ in range(num_runs):
    reload(DM_kde_gen)
    DM_frb_runs_ = np.percentile(DM_kde_gen.DM_frb,percentile_val)
    DM_frb_runs = np.append(DM_frb_runs,DM_frb_runs_)
DM_frb_uncertainty = np.std(DM_frb_runs)
print('The std on the simulated data is: ', DM_frb_uncertainty)
print('The average difference from true gap value on the simulated data is: ', np.mean(np.ones(len(DM_frb_runs))*100-DM_frb_runs))

"FIND GAP USING PDF MADE BY A KDE USING SMALL SET OF THE SIMULATED DATA"
DM_sim_draw_runs = []
for _ in range(num_runs):
    reload(DM_kde_gen)
    DM_sim_draw_runs_ = np.percentile(DM_kde_gen.DM_kde_draws_rand,percentile_val)
    DM_sim_draw_runs = np.append(DM_sim_draw_runs,DM_sim_draw_runs_)
print(len(DM_sim_draw_runs))
print('The lowest gap value is: ', min(DM_sim_draw_runs))
print('The highest gap value is: ', max(DM_sim_draw_runs))
DM_sim_draw_uncertainty = np.std(DM_sim_draw_runs)
print('The std on the KDE made by subset of simulated data is: ', DM_sim_draw_uncertainty)
print('The average difference from true gap value on the simulated data is: ', np.mean(np.ones(len(DM_sim_draw_runs))*100-DM_sim_draw_runs))

"FIND GAP USING PDF MADE BY KDE USING REAL DATA"
DM_kde_runs = []
for _ in range(num_runs):
    reload(DM_kde_gen)
    DM_kde_runs_ = np.percentile(DM_kde_gen.DM_kde_draws,percentile_val)
    DM_kde_runs = np.append(DM_kde_runs,DM_kde_runs_)
DM_kde_uncertainty = np.var(DM_kde_runs)

"FIND DIFFERENCE BETWEEN GAPS OF SIMULATED DATA AND PDF MADE USING KDE ON SMALL SUBSET OF SIMULATED DATA"
# DM_gap_uncertainty_runs = []
# for _ in range(num_runs):
#     reload(DM_kde_gen)
#     DM_gap_uncertainty_runs_ = min(DM_kde_gen.DM_frb)-min(DM_kde_gen.DM_kde_draws)
#     DM_gap_uncertainty_runs = np.append(DM_gap_uncertainty_runs,DM_gap_uncertainty_runs_)

# DM_gap_uncertainty = np.var(DM_gap_uncertainty_runs) #min(DM_gap_uncertainty_runs) #
# print(DM_gap_uncertainty_runs)
# print(DM_gap_uncertainty)

fig = plt.figure()

for _ in range(num_runs):
    reload(DM_kde_gen)
    plt.axvline(x=min(DM_kde_gen.DM_kde_draws_rand),linewidth=1,linestyle='dashed')
plt.plot(DM_grid, DM_kde_rand_distribution/np.sum(DM_kde_rand_distribution), linewidth=1, alpha=1, label=r'DM simulated KDE',c='purple')
"Plot DM_FRB and KDE FROM SIM"
plt.hist(DM_frb, density=True, bins=100, histtype='stepfilled', alpha=0.5,label=r'DM True simulated')
plt.hist(DM_kde_draws_rand, density=True, bins=len(DM_known_draws), histtype='stepfilled', alpha=0.4, label=r'DM simulated KDE',color='purple')
plt.xlabel('$DM$', fontsize=14)
plt.ylabel('PDF', fontsize=14)
plt.legend(fontsize=14)
plt.xlim(0,1500)
plt.show(block=True)

"SAVE OUTPUT GAPS"
np.save('DM_gap_outputs/DM_many_gaps_output.npy',DM_sim_draw_runs)