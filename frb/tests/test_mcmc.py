"""Test MCMC code"""
import os
import numpy as np
from pkg_resources import resource_filename


from matplotlib import pyplot as plt

from frb.dm import mcmc
from frb.dm import cosmic
from frb import utils
from frb import frb

import pytest

pm_required = pytest.mark.skipif(not hasattr(mcmc, 'pm'),
                                 reason='test requires pymc3')

# FRBs -- Init for the golden 5
frb180924 = frb.FRB.by_name('FRB20180924B')
frb181112 = frb.FRB.by_name('FRB20181112A')
frb190102 = frb.FRB.by_name('FRB20190102C')
frb190608 = frb.FRB.by_name('FRB20190608B')
frb190711 = frb.FRB.by_name('FRB20190711A')
mcmc.frbs = [frb180924, frb181112, frb190102, frb190608, frb190711]
mcmc.frb_DMs = np.array([frb.DM.value-frb.DMISM.value for frb in mcmc.frbs])
mcmc.frb_zs = np.array([frb.z for frb in mcmc.frbs])
DM_FRBp = mcmc.frb_DMs - mcmc.DM_MWhalo
mcmc.DM_FRBp_grid = np.outer(np.ones(mcmc.DM_values.size), DM_FRBp)
mcmc.DMhost_grid = np.outer(mcmc.DM_values, (1+mcmc.frb_zs))  # Host rest-frame DMs
mcmc.DMvalues_grid =  np.outer(mcmc.DM_values, np.ones(DM_FRBp.size))
mcmc.Deltavalues_grid = np.outer(mcmc.Delta_values, np.ones(DM_FRBp.size))

def test_pdf():
    # Mainly testing the jit aspect
    F = 0.32
    f_C0_3 = cosmic.grab_C0_spline()

    # Randomish values
    nFRB = 1000
    z_FRB = np.random.uniform(low=0.1, high=0.7, size=nFRB)
    sigma = F / np.sqrt(z_FRB)
    Delta = np.random.uniform(low=0.7, high=1.5, size=nFRB)
    C0 = f_C0_3(sigma)

    # Run it
    PDF_Cosmic = mcmc.mcquinn_DM_PDF_grid(Delta, C0, sigma)

    # Once more
    nFRB = 1000
    z_FRB = np.random.uniform(low=0.1, high=0.7, size=nFRB)
    sigma = F / np.sqrt(z_FRB)
    Delta = np.random.uniform(low=0.7, high=1.5, size=nFRB)
    PDF_Cosmic = mcmc.mcquinn_DM_PDF_grid(Delta, C0, sigma)

def test_allprob():
    F=0.32
    # All
    like = mcmc.all_prob(mcmc.cosmo_Obh70, F, None,
                            mcmc.frb_zs)
    # One by one
    ln_like = 0.
    probs = []
    for frb in mcmc.frbs:
        prob = mcmc.one_prob(mcmc.cosmo_Obh70, F, 
                            frb.DM.value - frb.DMISM.value, frb.z,
                mu=150., lognorm_s=1., lognorm_floor=0.,
                beta=3., orig=False)
        ln_like += np.log(prob)
        probs.append(prob)

    assert np.isclose(like, ln_like)


@pm_required
def test_pm():
    # This takes 5min to run
    # Hiding this import in here
    import pymc3 as pm

    parm_dict = mcmc.grab_parmdict()
    outroot = os.path.join(resource_filename('frb', 'tests'), 
                           'files', 'mcmc')

    with mcmc.pm_four_parameter_model(parm_dict, beta=3.):
        # Sample
        #trace = pm.sample(40000, tune=2000) # This defaults to 4 chains
        trace = pm.sample(1000, tune=500) # This defaults to 4 chains
        # Save the traces -- Needs to be done before the plot
        pm.save_trace(trace, directory=outroot, overwrite=True)
        print("All done with the 4 parameter, beta=3 run ")
        # Save a plot
        plt.clf()
        _ = pm.plot_trace(trace)
        #plt.savefig(os.path.join(outroot, 'traceplot.png'))
        # Parameters
        jdict = utils.jsonify(parm_dict)
        utils.savejson(os.path.join(outroot, 'parms.json'), jdict, easy_to_read=True)