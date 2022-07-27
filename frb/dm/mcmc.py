""" Methods for MCMC analysis of the Macquart relation """
import numpy as np
import warnings

from numba import njit
from scipy.stats import lognorm
from scipy.special import hyp1f1, gamma
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from frb import defs
from frb.dm import cosmic, igm

try:
    import pymc3 as pm
except ImportError:
    warnings.warn("You need to install pymc3 to run mcmc codes")
else:
    import theano.tensor as tt
    from theano.compile.ops import as_op

DM_values = np.linspace(1., 7000., 7000)  # Only used in one_prob
Delta_values = np.linspace(1./400., 20., 20000)# / 400.  # Fiducial
DM_MWhalo = 50.  # Fixed value for Macquart+20

# Init as blank
frbs = []
frb_DMs = None #np.array([frb.DM.value-frb.DMISM.value for frb in frbs])
frb_zs = None #np.array([frb.z for frb in frbs])

# The following are nomore used in the new likelihood function.
DM_FRBp_grid = None
DMhost_grid = None #np.outer(DM_values, (1+z_FRB))  # Host rest-frame DMs
DMvalues_grid = None # np.outer(DM_values, np.ones(DM_FRBp.size))
Deltavalues_grid = None  #np.outer(Delta_values, np.ones_like(C0)),

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
        #C0 = tt_spl_C0(sigma)
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


def log_likelihood(Obh70, F, in_DM_FRBp, z_FRB, mu=100., lognorm_s=1.,
                   lognorm_floor=0., beta=3., step=1.):
    """Calculate the probability for a set of FRBs

    Compared to the previously used "all_prob", this function includes a
    factor of <DM_cosmic> that comes from changeing the integration to Delta,
    i.e. from dDM = dDelta * <DM_cosmic>. On top, this function has a
    significant speed up due to looping instead of leaving parts of an
    array empty.

    Args:
        Obh70 (float): Value of Omega_b * h_70
        F (float): Feedback parameter
        in_DM_FRBp (np.ndarray): Measured DMs with DM_MW already subtracted.
        z_FRB (np.ndarray): FRB redshifts.
        mu (float, optional):
            Mean of log-normal PDF of DM_host (in DM units).
        lognorm_s (float, optional):
            Sigma of log-normal PDF of DM_host (in log space).
        beta (float, optional):
            Parameter for DM PDF.
        step (float, optional):
            Step size in DM units for the integral over the DM.

    Returns:
        float:  Log likelihood
    """
    # Sigma
    sigma = F / np.sqrt(z_FRB)
    # C0
    if beta == 4.:
        raise ValueError("Bad beta value in all_prob")
        #C0 = tt_spl_C0(sigma)
    elif beta == 3.:
        #C0 = tt_spl_C0_3(sigma)
        C0 = f_C0_3(sigma)
    else:
        raise IOError

    likelihoods = np.zeros(z_FRB.shape)

    avgDM = spl_DMc(z_FRB) * (Obh70 / cosmo_Obh70)

    for i, in_DM in enumerate(in_DM_FRBp):
        DM_values = np.arange(step/2, in_DM, step)
        DMcosmic = in_DM - DM_values
        Delta = DMcosmic / avgDM[i]
        PDF_Cosmic = cosmic.DMcosmic_PDF(Delta, C0[i], sigma[i], beta=beta)
        PDF_host = lognorm.pdf(DM_values*(1+z_FRB[i]), lognorm_s, lognorm_floor, mu)
        likelihoods[i] = step*np.sum(PDF_Cosmic * PDF_host)

    if beta == 3.:
        # Calculate the normalization "analytically" to be fast.
        hyp_x = -C0**2/18/sigma**2
        normalizations = 3*(12*sigma)**(1/3)/(gamma(1/3)*3*sigma*hyp1f1(1/6, 1/2, hyp_x)
                                             + 2**(1/2)*C0*gamma(5/6)*hyp1f1(2/3, 3/2, hyp_x))
    else:
        # Integrate numerically. This is slower by a factor 146 (with 20000 samples).
        step = 20/20000
        Delta = np.linspace(step, 20.-step, 20000)
        normalizations = cosmic.DMcosmic_PDF(Delta, C0[:, np.newaxis], sigma[:, np.newaxis], beta=beta)
        normalizations = 1 / (step * normalizations.sum(axis=-1))

    log_like = np.sum(np.log(likelihoods*normalizations/avgDM))
    return log_like


def log_likelihood_variable_step(Obh70, F, in_DM_FRBp, z_FRB, mu=100.,
                                 lognorm_floor=0, lognorm_s=1., beta=3.,
                                 res=500):
    """Calculate the log likelihood for a set of FRBs.

    Compared to the previously used "all_prob", this function includes a
    factor of <DM_cosmic> that comes from changeing the integration to Delta,
    i.e. from dDM = dDelta * <DM_cosmic>. On top, this function has a
    significant speed up due to a variable stepsize and the use of
    broadcasting. If an FRBs DM is too high you should increase the
    resolution.

    Args:
        Obh70 (float): Value of Omega_b * h_70
        F (float): Feedback parameter
        in_DM_FRBp (np.ndarray): Measured DMs with DM_MW already subtracted.
        z_FRB (np.ndarray): FRB redshifts.
        mu (float, optional):
            Mean of log-normal PDF of DM_host (in DM units).
        lognorm_s (float, optional):
            Sigma of log-normal PDF of DM_host (in log space).
        beta (float, optional):
            Parameter for DM PDF.
        res (int, optional):
            Number of steps to use for the integral over the DM.
            In a test with 200 FRBs the difference between res=100 and
            res=10000 was only 0.3%.

    Returns:
        float:  Log likelihood
    """
    # Sigma for each FRB.
    sigma = F / np.sqrt(z_FRB)
    # C0 for each FRB.
    if beta == 4.:
        raise ValueError("Bad beta value in all_prob")
        #C0 = tt_spl_C0(sigma)
    elif beta == 3.:
        #C0 = tt_spl_C0_3(sigma)
        C0 = f_C0_3(sigma)
    else:
        raise IOError

    # Get the average DM for each z from the globally created spline.
    avgDM = spl_DMc(z_FRB) * (Obh70 / cosmo_Obh70)

    # Integrate over P_host and P_cosmic in eq. 7 of Macquart et al. 2020 using the rectangle rule.
    steps = in_DM_FRBp/(res+1)  # Integration steps
    DM_values = np.linspace(steps/2, in_DM_FRBp-steps/2, res, axis=-1)  # 0th axis are the FRBs.

    Delta = (in_DM_FRBp[:, np.newaxis] - DM_values) / avgDM[:, np.newaxis]
    PDF_Cosmic = cosmic.DMcosmic_PDF(Delta, C0[:, np.newaxis], sigma[:, np.newaxis], beta=beta)
    PDF_host = lognorm.pdf(DM_values*(1+z_FRB[:, np.newaxis]), s=lognorm_s, scale=mu)
    likelihoods = steps*np.sum(PDF_Cosmic * PDF_host, axis=-1)

    if beta == 3.:
        # Calculate the normalization "analytically" to be fast.
        hyp_x = -C0**2/18/sigma**2
        normalizations = 3*(12*sigma)**(1/3)/(gamma(1/3)*3*sigma*hyp1f1(1/6, 1/2, hyp_x)
                                             + 2**(1/2)*C0*gamma(5/6)*hyp1f1(2/3, 3/2, hyp_x))
    else:
        # Integrate numerically. This is slower by a factor 146 (with 20000 samples).
        step = 20/20000
        Delta = np.linspace(step, 20.-step, 20000)
        normalizations = cosmic.DMcosmic_PDF(Delta, C0[:, np.newaxis], sigma[:, np.newaxis], beta=beta)
        normalizations = 1 / (step * normalizations.sum(axis=-1))

    # Normalization matters because it is different for each FRB. The factor avgDM comes because
    # normalizations is the integral over Delta instead of DM. avgDM was missing in previous
    # versions.
    log_like = np.sum(np.log(likelihoods*normalizations/avgDM))
    return log_like


try:
    @as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar], otypes=[tt.dvector])
    def calc_likelihood_four_beta3(Obh70, F, mu, lognorm_s):
        """Calculate likelihood for the real data

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
            ln_like = log_likelihood_variable_step(Obh70, F, frb_DMs, frb_zs, #frb.DM.value - frb.DMISM.value, frb.z,
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
    if (DM_FRBp_grid is not None
        or DMhost_grid is not None
        or DMvalues_grid is not None
        or Deltavalues_grid is not None):
        warnings.warn("You set one of the grid variables. That has no effect anymore.")
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
