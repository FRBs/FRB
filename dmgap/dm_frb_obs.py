""" Definitions relating to the observed DM distributions of pulsars and FRBs.
Standard KDE is used for pulsars and varying assymetric KDE for FRBs.
To use this code, please download the asymmetric_kde package from
https://github.com/tillahoffmann/asymmetric_kde [Hoffman and Jones, 2015] """

from asymmetric_kde import ProperGammaEstimator
from frb import igm
import numpy as np
import scipy as sp
import pandas as pd
from pandas import DataFrame as df
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import lognorm
from scipy.signal import argrelextrema
from pdf_fns import make_kde_funtion, make_pdf, rv_amount

from astropy.cosmology import WMAP9 as cosmo
from astropy import units
import random

def kde_for_frb(grid, dm_frb, num_samples, num_resamples=100, save_to_path=None):
    """
    Find FRB KDE and minima
    
    Arguments:
        dm_grid (array):
            Value of DM values.
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
    frb_lim = []
    dm_kde = []
    for _ in range(num_resamples):
        dm_sample = np.sort(np.asarray(random.sample(list(dm_frb),num_samples)))
        ige_sample = ProperGammaEstimator(dm_sample, 'plugin')
        ige_func = ige_sample(dm_grid) #KDEs
        frb_p = np.diff(ige_func)
        frb_lim_sample = dm_grid_p[argrelextrema(frb_p, np.greater)[0][0]] #max gradient at front of PDF
        frb_lim = np.append(frb_lim,frb_lim_sample)
        dm_kde.append(ige_func)

    if save_to_path==None:
        pass
    else:
        np.save(save_to_path, frb_lim)
        np.save(save_to_path, dm_kde)

    return dm_kde, frb_lim

def kde_for_psr(grid, dm_psr, num_samples, num_resamples=100, save_to_path=None):
    """
    Find pulsar KDE and maxima
    
    Arguments:
        dm_grid (array):
            Value of DM values.
        dm_psr (array):
            PDF of DM_psr.
        num_samples (int):
            Number of draws from dm_frb_obs.
        num_resamples (int, optional):
            Number of times to resample the draws. Default=100.
        save_to_path (str, optional):
            To save data outputs (as .npy), specify path. Default=None.

    Outputs:
        dm_kde (array):
            KDEs for observed pulsar DM distribution.
            A sample of size=num_samples is taken from the FRB PDF and then resampled num_resamples times.
        psr_lim (array):
            The max gradients of the pulsar DM KDEs, corresponding to an lower limit on the MW halo DM

    """
    dm_frb = dm_obs
    dm_grid = grid
    dm_grid_p = dm_grid[1:] #grid for the derivative 
    frb_lim = []
    dm_kde = []
    for _ in range(num_resamples):
        dm_sample = np.sort(np.asarray(random.sample(list(dm_frb),num_samples)))
        ige_sample = ProperGammaEstimator(dm_sample, 'plugin')
        ige_func = ige_sample(dm_grid) #KDEs
        frb_p = np.diff(ige_func)
        frb_lim_sample = dm_grid_p[argrelextrema(frb_p, np.greater)[0][0]] #max gradient at front of PDF
        frb_lim = np.append(frb_lim,frb_lim_sample)
        dm_kde.append(ige_func)

    if save_to_path==None:
        pass
    else:
        np.save(save_to_path, frb_lim)
        np.save(save_to_path, dm_kde)

    return dm_kde, frb_lim