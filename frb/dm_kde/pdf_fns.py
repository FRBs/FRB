""" Functions to for making PDFs"""

import numpy as np
import scipy as sp
import pandas as pd
import sklearn
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity

rv_amount = 10**6

def make_pdf(distribution, num_of_draws, grid, stepsize):
    """
    Makes PDF of given distribution
    
    Arguments:
        distribution (array):
            Array of values describing PDF.
        num_of_draws (int):
            Number of samples for PDF.
        grid (array):
            Desired grid for PDF.
        stepsize (float):
            Stepsize of desired grid.

    Outputs:
        draws (array):
            Samples drawn from PDF.
        distribution_scaled (array):
            PDF corresponding to input grid.

    """
    x_grid = np.arange(len(grid)) #rv_discrete only accepts interger values
    values = grid
    pdf = sp.stats.rv_discrete(values=(x_grid, distribution))
    draws_ = pdf.rvs(size=num_of_draws)
    draws = values[draws_] #rescale to floats
    distribution_scaled = distribution/stepsize
    return draws, distribution_scaled

def make_kde_funtion(grid, draws, min_bandwidth, max_bandwidth, bandwidth_stepsize, cv, kernel):
    " cv is number of cross-validation folds "
    """
       Returns KDE distribution
    
    Arguments:
        grid (array):
            Grid for PDF.
        draws (array):
            Sample from which to approximate PDF
        min_bandwidth (float):
            Start of bandwidth search range.
        max_bandwidth (float):
            End of bandwidth search range.
        bandwidth_stepsize (float):
            Stepsize for bandwidth search.
        cv (int):
            Number of folds for cross-validation
        kernel (str):
            Kernel to use. Valid kernels are 'gaussian', 'tophat', 'epanechnikov', 'exponential', 'linear', 'cosine'.
    Outputs:
        kde (array):
            PDF approximated by KDE
    """
    draws = np.asarray(draws)
    params = {'bandwidth': np.arange(min_bandwidth, max_bandwidth, bandwidth_stepsize)}
    grid_cv = GridSearchCV(KernelDensity(kernel=kernel), params, cv=cv)
    grid_cv.fit(draws.reshape(-1,1))
    bandwidth_opt= grid_cv.best_estimator_.bandwidth
    # print('Optimal bandwidth is:',bandwidth_opt)
    kde_skl = KernelDensity(kernel=kernel,bandwidth=bandwidth_opt)
    kde_skl.fit(draws.reshape(-1,1))
    log_kde = kde_skl.score_samples(grid.reshape(-1,1))
    kde = np.exp(log_kde)
    return kde