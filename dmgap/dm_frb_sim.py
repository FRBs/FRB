""" Definitions relating to the simulated DM distribution of FRBs.
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

def simulate_frb_dm(frb_data, z_stepsize=10**-4, z_max=3, dm_min=-10, dm_max=5000, dm_stepsize=10**-2, save_to_path=None):
    """
    Simulate DM_FRB by constructing PDFs for DM_cosmic, DM_host and DM_halo.
    Please check your frbcat_df matches the template used here.
    Change the halo and host distributions as you feel appropriate.

    Arguments:
        transient_data (str):
            Path to data (in .csv format).
        z_stepsize (float, optional):
            Redshift stepsize. Default=10**-4.
        z_max (int, optimal):
            Maximum redshift value. Default=3.
        dm_min (int, optional):
            Minimum value of DM grid. Default=-10.
        dm_max (int, optional):
            Minimum value of DM grid. Default=5000.
        dm_stepsize (float, optional):
            DM_grid stepsize. Default=10**-2.
        save_to_path (str, optional):
            To save data outputs (as .npy), specify path. Default=None.

    Outputs:
        dm_grid (array):
            DM grid
        dm_frb_sim (array):
            PDF of DM_FRB

    """
    # Estimate redshift distribution from observed FRB DMs (corrected for ISM)
    frbcat_df = pd.read_csv(frb_data)
    dm_frb = np.array(frbcat_df['deltaDM'])*units.pc/units.cm**3
    z_grid = np.arange(0,z_max,z_stepsize)
    z = []
    for i in range(len(dm_frb)):
        z_ = igm.z_from_DM(dm_frb[i], corr_nuisance=False)
        z = np.append(z,z_)

    # Kernel density estimation
    z_func_ = make_kde_funtion(grid=z_grid, draws=z, min_bandwidth=0.01, max_bandwidth=0.5, bandwidth_stepsize=0.1, cv=5, kernel='gaussian')
    z_func = z_func_/np.sum(z_func_)
    z = make_pdf(distribution=z_func,num_of_draws=rv_amount,grid=z_grid,stepsize=z_stepsize)
    z_draws = z[0]

    # Calculate average cosmic DM
    dm_ave_fn = igm.average_DM(z_max, cosmo=None, cumul=True, neval=rv_amount+1, mu=1.3)
    dm_ave_ = dm_ave_fn[0].value
    z_ave = dm_ave_fn[1]
    dm_int_ave = sp.interpolate.CubicSpline(z_ave,dm_ave_)
    dm_ave = np.sort(dm_int_ave(z_draws))
    
    # Find cosmic DM
    F = 0.2
    sigma_dm = np.minimum(F*z_draws**(-0.5),0.5)
    sigma_dm = [i for i in sigma_dm if i<0.5]
    sigma_gauss = 1
    rand_draws_gauss = np.random.normal(loc=0, scale=sigma_gauss, size=rv_amount)
    rand_draws_gauss = [i for i in rand_draws_gauss if i>-sigma_gauss and i<sigma_gauss]
    sigma_dm = random.sample(sigma_dm,len(rand_draws_gauss))
    dm_ave = np.array(random.sample(list(dm_ave),len(rand_draws_gauss)))
    dm_cos = dm_ave + rand_draws_gauss*dm_ave*sigma_dm
    dm_cos = [i for i in dm_cos if i > 0]

    # Estimate halo
    dm_halo  = 30. #some assigned MW halo DM 

    # Host DM - adjust this to your prefered (estimated) host DM
    mu_host = 40.
    sigma_host = 0.5
    dm_host = lognorm(s=sigma_host, loc=-4, scale=mu_host).rvs(size=len(dm_cos))

    dm_cos = np.array(random.sample(list(dm_cos),len(dm_host)))

    dm_frb_sim = np.sort(dm_cos)+dm_halo+np.sort(dm_host)
    dm_grid = np.arange(dm_min,dm_max,dm_stepsize)

    if save_to_path==None:
        pass
    else:
        np.save(save_to_path+'/dm_frb_sim.npy', dm_frb_sim)
        np.save(save_to_path+'/dm_grid.npy', dm_grid)

    return dm_grid, dm_frb_sim

def kde_for_frb_sim(grid, dm_sim, num_samples, num_resamples=100, save_to_path=None):
    """
    Estimate simulated DM_FRB using a random draws from the simulation. Used to analyse KDE effectiveness
    for different sample sizes (num_samples).
    
    Arguments:
        dm_grid (array):
            Value of DM values.
        dm_frb_sim (array):
            PDF of DM_FRB
        num_samples (int):
            Number of draws from dm_frb_sim.
        num_resamples (int, optional):
            Number of times to resample the draws. Default=100.
        save_to_path (str, optional):
            To save data outputs (as .npy), specify path. Default=None.

    Outputs:
        dm_kde (array):
            KDEs for simutated FRB DM distributions.
            A sample of size=num_samples is taken from the simulated PDF and then resampled num_resamples times.
        frb_lim (array):
            The max gradients of the FRB DM KDEs, corresponding to an upper limit on the MW halo DM

    """
    dm_frb_sim = np.load(dm_sim)
    dm_grid = np.load(grid)
    dm_grid_p = dm_grid[1:] #grid for the derivative 
    frb_lim = []
    dm_kde = []
    for _ in range(num_resamples):
        dm_sample = np.sort(np.asarray(random.sample(list(dm_frb_sim),num_samples)))
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