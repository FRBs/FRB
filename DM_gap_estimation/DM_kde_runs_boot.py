import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import sklearn
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
# from DM_definitions_boot import DM_known_draws, FRB_distribution_from_SFR, make_pdf, make_kde_funtion
import DM_kde_gen_boot
import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

"SAVE MULTIPLE OUTPUTS"


#timestamp
ts = str(int(time.time()))

# np.save('DM_gap_outputs/DM_frb_sim_npy/DM_frb_sim_hist'+ts+'.npy',DM_kde_gen.DM_frb)
#np.save('DM_gap_outputs/DM_rand_npy/rand_hist/DM_rand_hist'+ts+'.npy',DM_kde_gen.DM_kde_draws_rand)
# np.save('DM_gap_outputs/DM_known_npy/known_hist/DM_known_hist'+ts+'.npy',DM_kde_gen_boot.DM_draws)

# np.save('DM_gap_outputs/DM_rand_boot/100/DM_rand_boot_dist'+ts+'.npy',DM_kde_gen_boot.DM_kde_rand_distribution)
# np.save('DM_gap_outputs/DM_known_npy/known_distri/DM_known_dist'+ts+'.npy',DM_kde_gen.DM_kde_distribution)
np.save('DM_gap_outputs/DM_known_npy/known_distri/1000/'+ts+'.npy',DM_kde_gen_boot.DM_known_rand_draws)

