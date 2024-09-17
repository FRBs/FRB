# Imports

import matplotlib.pyplot as plt
import numpy as np
import os
from pandas import read_csv, read_table
from prospect.models import priors, SedModel
from prospect.models.sedmodel import PolySedModel
from prospect.models.templates import TemplateLibrary
from prospect.sources import CSPSpecBasis
from sedpy.observate import load_filters
from prospect.utils.obsutils import fix_obs
from scipy.stats import truncnorm

#####################################################################

# Adding mass-metallicity relation

class MassMet(priors.Prior):
    """
    A Gaussian prior designed to approximate the Gallazzi et al. 2005 
    stellar mass--stellar metallicity relationship.
    """

    prior_params = ['mass_mini', 'mass_maxi', 'z_mini', 'z_maxi']
    distribution = truncnorm
    massmet = np.loadtxt('gallazzi_05_massmet.txt')
    def __len__(self):
        """ Hack to work with Prospector 0.3
        """
        return 2

    def scale(self,mass):
        upper_84 = np.interp(mass, self.massmet[:,0], self.massmet[:,3]) 
        lower_16 = np.interp(mass, self.massmet[:,0], self.massmet[:,2])
        return (upper_84-lower_16)

    def loc(self,mass):
        return np.interp(mass, self.massmet[:,0], self.massmet[:,1])

    def get_args(self,mass):
        a = (self.params['z_mini'] - self.loc(mass)) / self.scale(mass)
        b = (self.params['z_maxi'] - self.loc(mass)) / self.scale(mass)
        return [a, b]

    @property
    def range(self):
        return ((self.params['mass_mini'], self.params['mass_maxi']),\
                (self.params['z_mini'], self.params['z_maxi']))

    def bounds(self, **kwargs):
        if len(kwargs) > 0:
            self.update(**kwargs)
        return self.range

    def __call__(self, x, **kwargs):
        """Compute the value of the probability density function at x and
        return the ln of that.

        :params x:
            x[0] = mass, x[1] = metallicity. Used to calculate the prior

        :param kwargs: optional
            All extra keyword arguments are used to update the `prior_params`.

        :returns lnp:
            The natural log of the prior probability at x, scalar or ndarray of
            same length as the prior object.
        """
        if len(kwargs) > 0:
            self.update(**kwargs)
        p = np.atleast_2d(np.zeros_like(x))
        a, b = self.get_args(x[...,0])
        p[...,1] = self.distribution.pdf(x[...,1], a, b, loc=self.loc(x[...,0]), scale=self.scale(x[...,0]))
        with np.errstate(invalid='ignore'):
            p[...,1] = np.log(p[...,1])
        return p

    def sample(self, nsample=None, **kwargs):
        """Draw a sample from the prior distribution.

        :param nsample: (optional)
            Unused
        """
        if len(kwargs) > 0:
            self.update(**kwargs)
        mass = np.random.uniform(low=self.params['mass_mini'],high=self.params['mass_maxi'],size=nsample)
        a, b = self.get_args(mass)
        met = self.distribution.rvs(a, b, loc=self.loc(mass), scale=self.scale(mass), size=nsample)

        return np.array([mass, met])

    def unit_transform(self, x, **kwargs):
        """Go from a value of the CDF (between 0 and 1) to the corresponding
        parameter value.

        :param x:
            A scalar or vector of same length as the Prior with values between
            zero and one corresponding to the value of the CDF.

        :returns theta:
            The parameter value corresponding to the value of the CDF given by
            `x`.
        """
        if len(kwargs) > 0:
            self.update(**kwargs)
        mass = x[0]*(self.params['mass_maxi'] - self.params['mass_mini']) + self.params['mass_mini']
        a, b = self.get_args(mass)
        met = self.distribution.ppf(x[1], a, b, loc=self.loc(mass), scale=self.scale(mass))
        return np.array([mass,met])
    
def massmet_to_mass(massmet=None, **extras):
    return 10**massmet[0]
def massmet_to_logzol(massmet=None,**extras):
    return massmet[1]

#####################################################################################################

# Adding dust1 and gas_logz functions

def dust2_to_dust1(dust2=None, **kwargs):
    return dust2

def gas_logz(gas_logz=None, **kwargs):
    return gas_logz

######################################################################################################

# Run parameters

run_params = {'verbose':True,
              'debug':False,
              'outfile':'FRB_190608_full_UV_newspec2',
              'output_pickles': False,
              # Optimization parameters
              'do_powell': False,
              'ftol':3e-16, 'maxfev': 5000,
              'do_levenburg': True,
              'nmin': 5,
              # dynesty Fitter parameters
              'nested_bound': 'multi', # bounding method
              'nested_sample': 'rwalk', # sampling method
              'nested_nlive_init': 500, # higher number increases speed of run
              'nested_nlive_batch': 500, # higher number increases speed of run
              'nested_bootstrap': 0,
              'nested_dlogz_init': 0.05,
              'nested_weight_kwargs': {"pfrac": 1.0},
              'nested_stop_kwargs': {"post_thresh": 0.1},
              # Obs data parameters
              'objid':190608,
              'logify_spectrum':False,
              'normalize_spectrum':False,
              'luminosity_distance': None,  # in Mpc
              # Model parameters
              'add_neb': True,
              'add_duste': True,
              # SPS parameters
              'zcontinuous': 1,
              }

# --------------
# OBS
# --------------
def load_spec_ascii(filename):
    f = read_table(filename, delimiter = '|', header=0, names=['NA', 'wave', 'flux', 'sig', 'NA2'], 
                   usecols =['wave', 'flux', 'sig'])
    wave = np.array(f['wave'])
    flux = np.array(f['flux'])
    err = np.array(f['sig'])
    return wave, flux, err

def load_spec_csv(filename):
    f = read_csv(filename, header=None, names=['wave', 'flux', 'err'], delimiter=',')
    wave = np.array(f['wave'])
    flux = np.array(f['flux'])
    err = np.array(f['err'])
    return wave, flux, err

# Convert from flux per Ang to flux per Hz
def convert(flux, wave, micro=True):
    """
    Convert from flux in per Angstrom to
    flux in per Hz (Jy). Note: wavelength
    vector must be in Ang.
    """
    new_flux = 3.34e4*(wave**2)*flux
    if micro == True:
        new_flux = new_flux*1e6
    return new_flux

def chop_wave(wave, w1, w2):
    cond_1 = np.where(w1 > wave)
    cond_2 = np.where(w2 < wave)
    return cond_1, cond_2

def load_obs(objid=run_params["objid"], phottable=None,
             luminosity_distance=0, snr=10, **kwargs):
    
    sdss = ['sdss_{}0'.format(b) for b in ['u', 'g', 'r', 'i', 'z']]
    hst = ['wfc3_ir_f160w', 'wfc3_uvis_f300x']
    wise = ['wise_w1', 'wise_w2', 'wise_w3', 'wise_w4']

    filternames = sdss + hst + wise
        
    # these are apparent magnitudes
    
    M_AB = np.array([18.9867385, 18.0239936, 17.5538606, 17.21585955, 17.087944049999997, 16.693-0.02, 19.765-0.30, 14.372155755376001+2.699, 13.834173935194169+3.339, 10.762158747454539+5.174, 8.647572481509055+6.620])


    # uncertainties in apparent magnitudes
    M_AB_unc = np.array([0.08660249, 0.01794288, 0.01485198, 0.01709021, 0.04910865, 0.001, 0.014, 0.03, 0.04, 0.115, 0.415]) 

    
    # convert from apparent magnitude to flux in maggies
    mags = 10**(-0.4*M_AB)

    # Find lower limits of mag uncertainty
    mag_down = [x-y for (x,y) in zip(M_AB, M_AB_unc)]

    # convert mag uncertainties to flux uncertainties
    flux_down = [10**(-0.4*x) for x in mag_down]
    flux_uncertainty = [y-x for (x,y) in zip(mags, flux_down)]

    # Build output dictionary.
    obs = {}
    # This is a list of sedpy filter objects.    See the
    # sedpy.observate.load_filters command for more details on its syntax.
    obs['filters'] = load_filters(filternames)
    # This is a list of maggies, converted from mags.  It should have the same
    # order as `filters` above.
    obs["phot_wave"] = np.array([f.wave_effective for f in load_filters(filternames)])
    obs['maggies'] = np.array(mags)
    obs['maggies_unc'] = np.array(flux_uncertainty)
    # Here we mask out any NaNs or infs
    obs['phot_mask'] = np.isfinite(np.squeeze(mags))
    
    # Adding in spectrum
    spec_wave, spec_fd, spec_fd_unc = load_spec_csv('FRB190608_corrected_spec.csv')
    
#    spec_fd = convert(spec_fd, spec_wave, micro=True)
#    spec_fd_unc = convert(spec_fd_unc, spec_wave, micro=True)

    
#    c1, c2 = chop_wave(spec_wave, 5575, 5581)
    
#    spec_wave = np.append(spec_wave[c1], spec_wave[c2])
#    spec_fd = np.append(spec_fd[c1], spec_fd[c2])
#    spec_fd_unc = np.append(spec_fd_unc[c1], spec_fd_unc[c2])
    
    obs['wavelength'] = spec_wave
    obs['spectrum'] = spec_fd * 1e-10
    obs['unc'] = spec_fd_unc * 1e-10

    # Add unessential bonus info.  This will be stored in output
    obs['objid'] = objid
   
    # Adding in mask of bad spectrum range
    obs['mask'] = ([])
    for i in range(len(obs['wavelength'])):
        obs['mask'].append(True)
    for j in range(1662, 1665):
        obs['mask'][j] = not obs['mask'][j]
    obs = fix_obs(obs)

    return obs

# --------------
# SPS Object
# --------------

def load_sps(zcontinuous=1, compute_vega_mags=False, **extras):
    sps = CSPSpecBasis(zcontinuous=zcontinuous,
                       compute_vega_mags=compute_vega_mags)
    return sps


# -----------------
# Gaussian Process
# ------------------

def load_gp(**extras):
    return None, None


# Load model function

def load_model(object_redshift=0.11778, add_duste=True, opt_spec=True, 
               add_dust1 = True, massmet = True, agn = True,
               add_neb=True, luminosity_distance=None, **extras):
    model_params = TemplateLibrary["parametric_sfh"]

    #fixed values
    model_params["imf_type"]["init"] = 1 # Chabrier
    model_params["dust_type"]["init"] = 1 # Milky Way extinction law
    model_params["sfh"]["init"] = 4 # non delayed-tau 
    model_params["logzsol"]["isfree"] = True
    model_params["tau"]["isfree"] = True
    model_params["dust2"]["isfree"] = True
    
    
    # adjust priors
    model_params["tau"]["prior"] = priors.LogUniform(mini=0.1, maxi=10)
    model_params["mass"]["prior"] = priors.LogUniform(mini=1e8, maxi=1e12)
    model_params["tage"]["prior"] = priors.TopHat(mini=0.0, maxi=12.196)
    model_params["logzsol"]["prior"] = priors.TopHat(mini=-0.5,maxi=0.19)
    model_params["dust2"]["prior"] = priors.TopHat(mini=0.0,maxi=2.0)
    #model_params["zred"]["prior"] = priors.TopHat(mini=1.0, maxi=1.5)

    # add AGN
    if agn:
        model_params.update(TemplateLibrary["agn"])
        model_params['agn_tau']['isfree'] = True # optical depth
        model_params['agn_tau']['prior'] = priors.LogUniform(mini=10.0, maxi=90.0)
        model_params['fagn']['isfree'] = True
        model_params['fagn']['prior'] = priors.LogUniform(mini=1e-5, maxi=2)
        model_params['add_dust_agn'] = {'N':1, 'init':True, 'isfree':False, "units":" ", 'prior':None}

    # Add dust1 ONLY if you have a star-forming galaxy
    if add_dust1:
        model_params['dust1'] = {'N':1, 'init':0.5, 'isfree':False,
                                'depends_on': dust2_to_dust1}
    # Setting redshift
    if object_redshift is not None:

        model_params["zred"]['isfree'] = False
        model_params["zred"]['init'] = object_redshift

    # Add nebular emission parameters and turn nebular emission on
    # Add gas_loz and gas_logu as free parameters if you have a spectrum
    if add_neb:
        model_params.update(TemplateLibrary["nebular"])
        model_params['gas_logu'] = {'N':1, 'init': -2, 'isfree':True,
                                    'prior': priors.TopHat(mini=-4, maxi=-1), 'units': 'Q_H/N_H'}
        model_params['gas_logz'] = {'N':1, 'init': 0.0, 'units': 'log Z/Z_\\odot', 'depends_on': gas_logz,
                                    'isfree':True, 'prior': priors.TopHat(mini=-2.0, maxi=0.5)}
        
    # Add dust emission in FIR 
    if add_duste:
        model_params.update(TemplateLibrary["dust_emission"])
        model_params["duste_gamma"]["isfree"] = True
        model_params["duste_gamma"]["prior"] = priors.LogUniform(mini=0.001, maxi=0.15)

        model_params["duste_qpah"]["isfree"] = True
        model_params["duste_qpah"]["prior"] = priors.TopHat(mini=0.5, maxi=7.0)
        
        model_params["duste_umin"]["isfree"] = True
        model_params["duste_umin"]["prior"] = priors.TopHat(mini=0.1,maxi=25)
        
    # Adding massmet param   
    if massmet:
        model_params['massmet'] = {"name": "massmet", "N": 2, "isfree": True, "init": [8.0, 0.0],
                     "prior": MassMet(z_mini=-0.5, z_maxi=0.10, mass_mini=9, mass_maxi=11)}
        model_params['mass']['isfree']=False
        model_params['mass']['depends_on']= massmet_to_mass
        model_params['logzsol']['isfree'] =False
        model_params['logzsol']['depends_on']=massmet_to_logzol

    if opt_spec:
        model_params.update(TemplateLibrary["optimize_speccal"])
        # fit for normalization of spectrum
        model_params['spec_norm'] = {'N': 1,'init': 1.0,'isfree': True,'prior': 
                                     priors.Normal(sigma=0.2, mean=1.0), 'units': 'f_true/f_obs'}
        # Increase the polynomial size to 12
        model_params['polyorder'] = {'N': 1, 'init': 10,'isfree': False}
        
        run_params["opt_spec"] = True
    
        # Now instantiate the model using this new dictionary of parameter specifications
        model = PolySedModel(model_params)
        
    elif opt_spec == False:
        model = SedModel(model_params)
        run_params["opt_spec"] = False
        
    return model
