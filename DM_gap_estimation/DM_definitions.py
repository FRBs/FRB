import numpy as np
import scipy as sp
import pandas as pd
import sklearn
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from astropy.cosmology import WMAP9 as cosmo

"VARIABLES"
rv_amount = 10**5

"DATA UPLOAD"
z_known_draws = np.asarray([0.1927, 0.3212, 0.4755, 0.29, 0.66, 0.1178, 0.378]) #observed redshift values
data_FRB_raw = pd.read_csv('FRB_data/frbcat_20190722.csv')
DM_known_draws = data_FRB_raw['rmp_dm'].str.split('&').str[0].astype(float).values #observed DM values
# print('The minimum DM is: ',min(DM_known_draws))
# print('The maxiumu DM is: ',max(DM_known_draws))

def FRB_distribution_from_SFR(z_grid):
    Mpc = cosmo.comoving_distance(z_grid).value #redshift conversion to Mpc
    psi = 0.015*(1+z_grid)**2.7/(1+((1+z_grid)/2.9)**5.6) #star formation rate density
    FRB_rate_density = psi*np.exp(-z_grid/0.5) #scale by exponential cut off
    V_sphere = 4/3*np.pi*Mpc**3
    V_shell = [y - x for x,y in zip(V_sphere,V_sphere[1:])]
    z_FRB_likelihood_SFR_ = FRB_rate_density[1:]*V_shell #FRB distribution density to distance
    z_FRB_likelihood_SFR = z_FRB_likelihood_SFR_/sum(z_FRB_likelihood_SFR_) #normalise
    return z_FRB_likelihood_SFR

def make_pdf(distribution,num_of_draws,grid,stepsize):
    "Makes PDF of given distribution"
    x_grid = np.arange(len(grid)) #rv_discrete only accepts interger values
    values = grid
    pdf = sp.stats.rv_discrete(values=(x_grid, distribution))
    draws_ = pdf.rvs(size=num_of_draws)
    draws = values[draws_] #rescale to floats
    distribution_scaled = distribution/stepsize
    return draws, distribution_scaled

def make_kde_funtion(grid, draws, min_bandwidth, max_bandwidth, bandwidth_stepsize, cv, kernel):
    "Returns KDE distribution"
    draws = np.asarray(draws)
    params = {'bandwidth': np.arange(min_bandwidth, max_bandwidth, bandwidth_stepsize)}
    grid_cv = GridSearchCV(KernelDensity(kernel=kernel), params, cv=cv)
    grid_cv.fit(draws.reshape(-1,1))
    bandwidth_opt= grid_cv.best_estimator_.bandwidth
    kde_skl = KernelDensity(kernel=kernel,bandwidth=bandwidth_opt)
    kde_skl.fit(draws.reshape(-1,1))
    log_kde = kde_skl.score_samples(grid.reshape(-1,1))
    kde = np.exp(log_kde)
    # print('best bandwidth: {0}'.format(bandwidth_opt),'for a {0}'.format(kernel), 'kernel')
    return kde