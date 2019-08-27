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
from DM_static_definitions import DM_known_draws, FRB_distribution_from_SFR, make_pdf, make_kde_funtion
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
num_kde_samples = 10*4
DM_frb = np.load('DM_outputs/DM_frb.npy')
DM_grid = np.arange(len(DM_frb))
DM_stepsize = 0.1
num_frb_draws = 1000
DM_rand_draws = np.asarray(random.sample(list(DM_frb),num_frb_draws))

"Make KDE distribution kernel with KNOWN DM values, then take loads of samples of this this distribution to make a PDF"
DM_kde_func_ = make_kde_funtion(grid=DM_grid, draws = DM_known_draws, bandwidth=99, kernel='gaussian')
DM_kde_func = DM_kde_func_/np.sum(DM_kde_func_) #normalised
DM_kde = make_pdf(distribution=DM_kde_func,num_of_draws=num_kde_samples,grid=DM_grid,stepsize=DM_stepsize)
DM_kde_draws = DM_kde[0]
DM_kde_distribution  = DM_kde[1] #DM_kde_func scaled

"Make KDE distribution kernel with RANDOM DM values, then take loads of samples of this this distribution to make a PDF"
DM_kde_rand_func_ = make_kde_funtion(grid=DM_grid, draws = DM_rand_draws, bandwidth=60, kernel='gaussian')
DM_kde_rand_func = DM_kde_rand_func_/np.sum(DM_kde_rand_func_) #normalised
DM_kde_rand = make_pdf(distribution=DM_kde_rand_func,num_of_draws=num_kde_samples,grid=DM_grid,stepsize=DM_stepsize)
DM_kde_draws_rand = DM_kde_rand[0]
DM_kde_rand_distribution  = DM_kde_rand[1]