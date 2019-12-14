import numpy as np
import scipy as sp
import pandas as pd
import sklearn
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity

import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
 
rv_amount = 10**6

def make_pdf(distribution,num_of_draws,grid,stepsize):
    # Makes PDF of given distribution
    x_grid = np.arange(len(grid)) #rv_discrete only accepts interger values
    values = grid
    pdf = sp.stats.rv_discrete(values=(x_grid, distribution))
    draws_ = pdf.rvs(size=num_of_draws)
    draws = values[draws_] #rescale to floats
    distribution_scaled = distribution/stepsize
    return draws, distribution_scaled

def make_kde_funtion(grid, draws, min_bandwidth, max_bandwidth, bandwidth_stepsize, cv, kernel):
    # Returns KDE distribution
    draws = np.asarray(draws)
    params = {'bandwidth': np.arange(min_bandwidth, max_bandwidth, bandwidth_stepsize)}
    grid_cv = GridSearchCV(KernelDensity(kernel=kernel), params, cv=cv)
    grid_cv.fit(draws.reshape(-1,1))
    bandwidth_opt= grid_cv.best_estimator_.bandwidth
    kde_skl = KernelDensity(kernel=kernel,bandwidth=bandwidth_opt)
    kde_skl.fit(draws.reshape(-1,1))
    log_kde = kde_skl.score_samples(grid.reshape(-1,1))
    kde = np.exp(log_kde)
    return kde

def find_kde_ensemble(grid, draws, min_bandwidth, max_bandwidth, bandwidth_stepsize, cv, kernel, n, reps):
    # Find ensemble via bootstrapping
    # n=number of samples, reps=number of ensembles to take 
    resample_set = np.sort(np.random.choice(draws, (n, reps),replace=True))
    kde_bootstrapped = []
    for i in range(len(draws)):
        kde_predictions = make_kde_funtion(grid=grid, draws=resample_set[i], min_bandwidth=min_bandwidth, max_bandwidth=max_bandwidth, bandwidth_stepsize=bandwidth_stepsize, cv=cv, kernel=kernel)
        kde_bootstrapped.append([kde_predictions])
    return kde_bootstrapped