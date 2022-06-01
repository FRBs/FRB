import os
import numpy as np
from frb.galaxies import prospector

from IPython import embed

datafile = os.path.join(
    os.getenv('FRB_GDB'), 'F4', 'gordon2022',
    'FRB20180924_NP_VISTA_mask_mcmc.h5')

# read in h5 file
results, obs_dict, model = prospector.results_from(datafile,
                                               dangerous=False)

# load parameter names
parnames = np.array(results['theta_labels'], dtype='U20')

# pull the chains and probability weights
samples = results['chain']
weights = results.get('weights', None)

# get representative sample of the posterior by sampling 100,000 times
theta_samp = prospector.sample_posterior(
    samples, weights=weights, nsample=100000)
embed(header='22 of test')

########################################################################################
# theta_samp returns a 100000x22 numpy array where the 0th axis is the sampled posterior
# and the 1st axis is each stellar population parameter. The parnames array lists which
# index is which parameter and is in the same order as the chain and theta_samp. 
# theta_samp is the element that we need to do statistics on to derive the median and one
# sigma uncertainties to report.
########################################################################################