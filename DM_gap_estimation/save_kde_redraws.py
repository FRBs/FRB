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
import DM_redraw_kde_gen
import seaborn as sns
sns.set()
sns.set_style('darkgrid')
sns.set_context('poster')

"SAVE MULTIPLE OUTPUTS"

#timestamp
ts = str(int(time.time()))

"Change number of draws in DM_redraw_kde_gen and then uncomment relevant path here (sorry!)"
np.save('DM_gap_outputs/DM_redraws/100/DM_redraw'+ts+'.npy',DM_redraw_kde_gen.DM_kde_rand_distribution)
# np.save('DM_gap_outputs/DM_redraws/1000/DM_redraw'+ts+'.npy',DM_redraw_kde_gen.DM_kde_rand_distribution)


