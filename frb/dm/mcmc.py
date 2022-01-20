""" Methods for MCMC analysis of the Macquart relation """
import numpy as np
from numba import njit 

from scipy.stats import lognorm
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from frb import defs

import warnings

try:
    import pymc3 as pm
except ImportError:
    warnings.warn("You need to install pymc3 to run mcmc codes")
else:
    import theano.tensor as tt
    from theano.compile.ops import as_op

from frb.dm import cosmic, igm

from IPython import embed

DM_values = np.linspace(1., 7000., 7000)
Delta_values = np.linspace(1./400., 20., 20000)# / 400.  # Fiducial
DM_MWhalo = 50.  # Fixed value for Macquart+20

# Init as blank
frbs = []
frb_DMs = None #np.array([frb.DM.value-frb.DMISM.value for frb in frbs])
frb_zs = None #np.array([frb.z for frb in frbs])
DM_FRBp_grid = None
DMhost_grid = None #np.outer(DM_values, (1+z_FRB))  # Host rest-frame DMs
DMvalues_grid = None # np.outer(DM_values, np.ones(DM_FRBp.size))
Deltavalues_grid = None  #np.outer(Delta_values, np.ones_like(C0)), 

#

# DM Cosmic
DM_cosmic, zeval = igm.average_DM(1., cosmo=defs.frb_cosmo, cumul=True)
spl_DMc = IUS(zeval, DM_cosmic.value)
cosmo_ObH0 = defs.frb_cosmo.Ob0 * defs.frb_cosmo.H0.value
cosmo_Obh70 = defs.frb_cosmo.Ob0 * (defs.frb_cosmo.H0.value/70.)

# Load splines
spl_sigma = cosmic.grab_sigma_spline()
f_C0_3 = cosmic.grab_C0_spline()


# Parameters
def grab_parmdict(tight_ObH=False):
    """ Generate the parameter dict for the MCMC run

    Args:
        tight_ObH (bool, optional): [description]. Defaults to False.

    Raises:
        IOError: [description]

    Returns:
        dict: [description]
    """
    parm_dict = {}
    parm_dict['F'] = dict(dist='Uniform', lower=0.011, upper=0.5, latex='F', unit='None')
    parm_dict['mu'] = dict(dist='Uniform', lower=20., upper=200, latex='\exp(\mu)', unit='pc cm$^{-3}$')
    parm_dict['lognorm_s'] = dict(dist='Uniform', lower=0.2, upper=2, latex='\sigma_{\\rm host}', unit='pc $^{-3}$')
    if tight_ObH:
        # NEED TO FIX THIS for h70!
        raise IOError
        parm_dict['ObH0'] = dict(dist='Normal', mu=cosmo_ObH0,
                                 sigma=cosmo_ObH0*0.03)
    else:
        parm_dict['Obh70'] = dict(dist='Uniform', lower=0.015, upper=0.095)
    parm_dict['Obh70']['latex'] = '\\Omega_b h_{70}'
    parm_dict['Obh70']['unit'] = 'None' #$\\rm km \, s^{-1} Mpc^{-1}$'
    # Other info
    parm_dict['Other'] = dict(DM_MWhalo=50., zscale_host=True, floor=0.)
    # Return    
    return parm_dict

def one_prob(Obh70, F, DM_FRBp, z_FRB, mu=150., lognorm_s=1.,
             lognorm_floor=0., orig=False, beta=4.):
    """
    Calculate the probability for a single FRB

    Args:
        Obh70 (float): Value of Omega_b * H_0 
        F (float): Feedback parameter
        DM_FRBp (np.ndarray): Values of DM_FRBp for analysis
        z_FRB (np.ndarray): z values for evaluation
        mu (float, optional):
            Mean of log-normal PDF
        lognorm_s (float, optional):
            Sigma of log-normal PDF
        lognorm_floor (float, optional):
            Floor to the log-normal PDF
        orig (bool, optional):
            if True (not recommended!), use the original approach to 
            calculating sigma
        beta (float, optional):
            Parameter for DM PDF

    Returns:
        float: Likelihood probability

    """
    # Update DM_FRBp with MW halo (global)
    DM_FRBp = DM_FRBp - DM_MWhalo

    # DM values that matter
    keep = DM_values <= DM_FRBp
    sub_DMvalues = DM_values[keep]

    # PDF Nuisance
    # Scale by redshift of the FRB
    PDF_Nuisance = lognorm.pdf(sub_DMvalues*(1+z_FRB), lognorm_s, lognorm_floor, mu)

    # Sigma
    if orig:
        sigma_DMp = F / np.sqrt(z_FRB)
        sigma = tt_spl_sigma(sigma_DMp)
    else:
        sigma = F / np.sqrt(z_FRB)
    # C0
    if beta == 4.:
        raise ValueError("Bad beta value")
        C0 = tt_spl_C0(sigma)
    elif beta == 3.:
        C0 = tt_spl_C0_3(sigma)
    else:
        raise IOError

    # Delta
    #avgDM = spl_DMc(z_FRB) * (ObH0 / cosmo_ObH0)
    avgDM = spl_DMc(z_FRB) * (Obh70 / cosmo_Obh70)
    sub_DMcosmic = DM_FRBp - sub_DMvalues
    Delta = sub_DMcosmic / avgDM

    # PDF DM_cosmic
    PDF_Cosmic = cosmic.DMcosmic_PDF(Delta, C0, sigma, beta=beta)
    all_PDF_Cosmic = cosmic.DMcosmic_PDF(Delta_values, C0, sigma, beta=beta)

    # Integrate
    Prob = np.sum(PDF_Nuisance * PDF_Cosmic) / np.sum(all_PDF_Cosmic)

    # Return
    return Prob#, sigma, C0


# C code or jit
@njit(parallel=False)
def mcquinn_DM_PDF_grid(Delta_values, C0, sigma, alpha=3., beta=3.):
    """
    PDF(Delta) for the McQuinn formalism describing the DM_cosmic PDF

    Args:
        Delta (2D ndarray):
            DM / averageDM values
        C0 (np.ndarray):
            C0 values
        sigma (np.ndarray):
            sigma values
        alpha (float, optional):
        beta (float, optional):

    Returns:

    """
    #
    Dalpha_grid = Delta_values**(-1*alpha)
    Dbeta_grid = Delta_values**(-1*beta)
    C0_grid = np.outer(np.ones(Delta_values.shape[0]), C0)
    sigma_grid = np.outer(np.ones(Delta_values.shape[0]), sigma)
    # Evaluate
    PDF_grid = np.exp(-(Dalpha_grid-C0_grid)**2 / (
            2*alpha**2 * sigma_grid**2))*Dbeta_grid
    # Return
    return PDF_grid


def all_prob(Obh70, F, in_DM_FRBp, z_FRB, mu=150., lognorm_s=1.,
             lognorm_floor=0., beta=3.):
    """
    Calculate the probability for a set of FRBs

    Args:
        Obh70 (float): Value of Omega_b * H_0 
        F (float): Feedback parameter
        in_DM_FRBp (np.ndarray): Values of DM_FRBp for analysis
            Not used?!
        z_FRB (np.ndarray): z values for evaluation
        mu (float, optional):
            Mean of log-normal PDF
        lognorm_s (float, optional):
            Sigma of log-normal PDF
        lognorm_floor (float, optional):
            Floor to the log-normal PDF
        beta (float, optional):
            Parameter for DM PDF

    Returns:
        float:  Log like-lihood

    """
    '''
    # Update DM_FRBp with MW halo (global)
    DM_FRBp = in_DM_FRBp - DM_MWhalo

    # DM values that matter
    DM_FRBp_grid = np.outer(np.ones(DM_values.size), DM_FRBp)
    DMhost_grid = np.outer(DM_values, (1+z_FRB))  # Host rest-frame DMs
    '''

    # PDF Nuisance
    # Scale by redshift of the FRB
    PDF_Nuisance_grid = lognorm.pdf(DMhost_grid, lognorm_s, 
                                    lognorm_floor, mu)

    # Sigma
    sigma = F / np.sqrt(z_FRB)
    # C0
    if beta == 4.:
        raise ValueError("Bad beta value in all_prob")
        C0 = tt_spl_C0(sigma)
    elif beta == 3.:
        #C0 = tt_spl_C0_3(sigma)
        C0 = f_C0_3(sigma)
    else:
        raise IOError

    # Delta
    avgDM = spl_DMc(z_FRB) * (Obh70 / cosmo_Obh70)
    DMcosmic = DM_FRBp_grid - DMvalues_grid #np.outer(DM_values, np.ones(DM_FRBp.size))
    Delta = DMcosmic / avgDM
    goodD = Delta > 0.

    # PDF DM_cosmic
    Delta[~goodD] = 1.  # To avoid bad things
    PDF_Cosmic = mcquinn_DM_PDF_grid(Delta, C0, sigma, beta=beta)
    PDF_Cosmic[~goodD] = 0.

    all_PDF_Cosmic = mcquinn_DM_PDF_grid(Deltavalues_grid, C0, sigma, beta=beta)

    # Integrate
    Prob = np.sum(PDF_Nuisance_grid * PDF_Cosmic, axis=0) / np.sum(
        all_PDF_Cosmic, axis=0)
    log_like = np.sum(np.log(Prob))

    # Return
    return log_like



try:
    @as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar], otypes=[tt.dvector])
    def calc_likelihood_four_beta3(Obh70, F, mu, lognorm_s):
        """
        Calculate likelihood for the real data

        Args:
            Obh70 (float): Value of Omega_b * H_0 
            F (float): Feedback parameter
            mu (float): Mean of log-normal PDF
            lognorm_s (float): Sigma of log-normal PDF

        Returns:
            np.ndarray:  Array of log likelihood values, one per FRB
                in the global variable frbs

        """

        lognorm_floor=0.
        # Loop on the FRBs
        ln_like = 0.
        one_by_one = False
        if one_by_one:
            for frb in frbs:
                prob = one_prob(Obh70, F, frb.DM.value - frb.DMISM.value, frb.z,
                            mu=mu, lognorm_s=lognorm_s, lognorm_floor=lognorm_floor,
                            beta=3., orig=False)
                ln_like += np.log(prob)
        else:
            # The loop above is much faster for only ~10 events
            ln_like = all_prob(Obh70, F, frb_DMs, frb_zs, #frb.DM.value - frb.DMISM.value, frb.z,
                        mu=mu, lognorm_s=lognorm_s, lognorm_floor=lognorm_floor, beta=3.)

        # Return
        return np.array([ln_like])   # Should be log
except: # Hiding this theano method in a try/except
    pass


def pm_four_parameter_model(parm_dict:dict, tight_ObH=False, beta=3.):
    """ Builds a pymc3 model for the 4-parameter MCMC

    Args:
        parm_dict (dict): dict with the pymc3 parameters
        tight_ObH (bool, optional): If True, restrict the ObH0 value based on CMB. Defaults to False.
        beta (float, optional): PDF parameter. Defaults to 3..

    Raises:
        IOError: [description]

    Returns:
        pm.Model: pymc3 model
    """
    # Load the model
    with pm.Model() as model:
        # Variables
        if tight_ObH:
            assert parm_dict['Obh70']['dist'] == 'Normal'
            Obh70 = pm.Normal('Obh70', mu=cosmo_Obh70, sigma=cosmo_Obh70*0.03)
        else:
            assert parm_dict['Obh70']['dist'] == 'Uniform'
            Obh70 = pm.Uniform('Obh70', lower=parm_dict['Obh70']['lower'],
                               upper=parm_dict['Obh70']['upper'])
        assert parm_dict['F']['dist'] == 'Uniform'
        F = pm.Uniform('F', lower=parm_dict['F']['lower'], upper=parm_dict['F']['upper'])
        mu = pm.Uniform('mu', lower=parm_dict['mu']['lower'], upper=parm_dict['mu']['upper'])
        lognorm_s = pm.Uniform('lognorm_s', lower=parm_dict['lognorm_s']['lower'], upper=parm_dict['lognorm_s']['upper'])
        # Potential
        if beta == 4.:
            like = pm.Potential('like', calc_likelihood_four(Obh70, F, mu, lognorm_s))
        elif beta == 3.:
            like = pm.Potential('like', calc_likelihood_four_beta3(Obh70, F, mu, lognorm_s))
        else:
            raise IOError
    return model


#@as_op(itypes=[tt.dscalar], otypes=[tt.dscalar])
def tt_spl_sigma(value):
    return float(spl_sigma(value))

def tt_spl_C0_3(value):
    return float(f_C0_3(value))
