""" Definitions relating to the observed DM distributions of pulsars and FRBs.
Standard KDE is used for pulsars and varying asymmetric KDE for FRBs.
To use this code, please download the asymmetric_kde package from
https://github.com/tillahoffmann/asymmetric_kde [Hoffman and Jones, 2015] """
from pkg_resources import resource_filename
from asymmetric_kde import ProperGammaEstimator
from pdf_fns import make_kde_funtion
from data import FRBs
import numpy as np
import pandas as pd
from scipy.signal import argrelextrema
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
import random
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def kde_for_frb(grid, dm_frb, num_samples, num_resamples=100, save_to_path=None):
    """
    Find FRB KDE and minima
    
    Arguments:
        dm_grid (array):
            DM values.
        dm_frb (array):
            PDF of DM_FRB.
        num_samples (int):
            Number of draws from dm_frb.
        num_resamples (int, optional):
            Number of times to resample the draws. Default=100.
        save_to_path (str, optional):
            To save data outputs (as .npy), specify path. Default=None.

    Outputs:
        dm_kde (array):
            KDEs for observed FRB DM distribution.
            A sample of size=num_samples is taken from the FRB PDF and then resampled num_resamples times.
        frb_lim (array):
            The max gradients of the FRB DM KDEs, corresponding to an lower limit on the MW halo DM

    """
    dm_grid = grid
    dm_grid_p = dm_grid[1:] #grid for the derivative 
    logging.info('Finding optimal KDE...')
    ige = ProperGammaEstimator(dm_frb, 'plugin')
    dm_kde_opt = ige(dm_grid) #KDEs
    frb_p_opt = np.diff(dm_kde_opt)
    frb_lim_opt = dm_grid_p[argrelextrema(frb_p_opt, np.greater)[0][0]] 
    frb_lim = []
    dm_kde = []
    logging.info('Bootstrapping KDE...')
    for _ in range(num_resamples):
        dm_sample = np.sort(np.asarray(random.sample(list(dm_frb),num_samples)))
        ige_sample = ProperGammaEstimator(dm_sample, 'plugin')
        ige_func = ige_sample(dm_grid) #KDEs
        frb_p = np.diff(ige_func)
        frb_lim_sample = dm_grid_p[argrelextrema(frb_p, np.greater)[0][0]] #max gradient at front of PDF
        frb_lim.append(frb_lim_sample)
        dm_kde.append(ige_func)

    if save_to_path==None:
        pass
    else:
        np.save(save_to_path+'/kde_frb_lim.npy', frb_lim)
        np.save(save_to_path+'/kde_frb_lim_opt.npy', frb_lim_opt)
        np.save(save_to_path+'/kde_frb.npy', dm_kde)
        np.save(save_to_path+'/kde_frb_opt.npy', dm_kde_opt)
        logging.info('Output saved to '+str(save_to_path))

    return dm_kde, frb_lim

def kde_for_psr(grid, dm_psr, num_samples, num_resamples=100, min_bandwidth=8, max_bandwidth=15, bandwidth_stepsize=0.1, cv=100, kernel='gaussian', save_to_path=None):
    """
    Find pulsar KDE and maxima
    
    Arguments:
        dm_grid (array):
            Value of DM values.
        dm_psr (array):
            PDF of DM_psr.
        num_samples (int):
            Number of draws from dm_psr
        num_resamples (int, optional):
            Number of times to resample the draws. Default=100.
        min_bandwidth (int, optional):
            Minimum bandwidth for cross-validation. Default=8.
        max_bandwidth (int, optional):
            Maximum bandwidth for cross-validation. Default=15.
        bandwidth_stepsize (int, optional):
            Bandwidth stepsize. Default=1.
        cv (int, optional):
            Number of folds for cross-validation. Default=5.
        kernel (str, optional):
            Kernel to use. Default='gaussian'
        save_to_path (str, optional):
            To save data outputs (as .npy), specify path. Default=None.

    Outputs:
        dm_kde (array):
            KDEs for observed pulsar DM distribution.
            A sample of size=num_samples is taken from the FRB PDF and then resampled num_resamples times.
        psr_lim (array):
            The max gradients of the pulsar DM KDEs, corresponding to an lower limit on the MW halo DM

    """
    dm_grid = grid
    dm_grid_p = dm_grid[1:] #grid for the derivative
    logging.info('Finding optimal KDE...')
    dm_kde_opt = make_kde_funtion(grid=dm_grid, draws=dm_psr, min_bandwidth=min_bandwidth, max_bandwidth=max_bandwidth, bandwidth_stepsize=bandwidth_stepsize, cv=cv, kernel=kernel)
    psr_p_opt = np.diff(dm_kde_opt)
    psr_lim_opt = dm_grid_p[psr_p_opt.argmin()]
    psr_lim = []
    dm_kde = []
    with NumpyRNGContext(1):
        boot_obs_psr = bootstrap(dm_psr, num_resamples)
    logging.info('Bootstrapping KDE...')
    for i in range(num_resamples):
        kde_func = make_kde_funtion(grid=dm_grid, draws=boot_obs_psr[i], min_bandwidth=min_bandwidth, max_bandwidth=max_bandwidth, bandwidth_stepsize=bandwidth_stepsize, cv=cv, kernel=kernel)
        psr_p = np.diff(kde_func)
        psr_lim_sample = dm_grid_p[psr_p.argmin()]
        psr_lim.append(psr_lim_sample)
        dm_kde.append(kde_func)

    if save_to_path==None:
        pass
    else:
        np.save(save_to_path+'/kde_psr_lim.npy', psr_lim)
        np.save(save_to_path+'/kde_psr_lim_opt.npy', psr_lim_opt)
        np.save(save_to_path+'/kde_psr.npy', dm_kde)
        np.save(save_to_path+'/kde_psr_opt.npy', dm_kde_opt)
        logging.info('Output saved to '+str(save_to_path))

    return dm_kde, psr_lim