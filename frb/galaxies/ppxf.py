""" Module for running pPXF analyses"""

from pkg_resources import resource_filename

import numpy as np

from matplotlib import pyplot as plt
#from goodies import closest
from astropy import constants
from astropy import units
from astropy.table import Table

c = constants.c.to(units.km / units.s).value

from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra.io import readspec

from ppxf import ppxf
from ppxf import ppxf_util as util
from ppxf import miles_util as lib
import time

from frb.defs import frb_cosmo as cosmo 

from IPython import embed

def run(spec_file, R, zgal, results_file=None, spec_fit='tmp.fits', chk=True,
        flux_scale=1., atmos=[], gaps=[], wvmnx=(0.,1e9)):
    """
    Wrapper for running and handling outputs

    Outputs are written to disk

    Args:
        spec_file (str or XSpectrum1D):
        R (float):
        zgal (float):
        results_file (str, optional):
        spec_fit:
        chk:
        flux_scale:
        atmos (list of tuple):
            List of (wvmin,wvmax) regions to mask during the analysis
            This is a list of lists, e.g.  [[7150., 7300.]]
        gaps (list):
            Regions to ignore due to detector gaps or any other bad regions
            This is a list of lists, e.g.  [[6675., 6725.]]
        wvmnx:
    """
    # Init

    # Load spectrum
    if isinstance(spec_file, XSpectrum1D):
        spec = spec_file
    else:
        spec = readspec(spec_file)

    if chk:
        spec.plot()

    # Rebin
    wave = spec.wavelength.value

    diff = wave[1:] - wave[0:-1]
    meddiff = np.median(diff)
    print(meddiff)
    newwave = np.arange(wave[0], wave[-2], meddiff) * units.angstrom
    newspec = spec.rebin(newwave, do_sig=True, grow_bad_sig=True)

    # Scale to MUSE flux units
    newspec = XSpectrum1D.from_tuple((newspec.wavelength.value,
                                      newspec.flux * flux_scale,
                                      newspec.sig * flux_scale))

    # Mask
    wave = newspec.wavelength.value
    goodpixels = np.where((wave >= wvmnx[0]) & (wave <= wvmnx[1]))[0]

    if atmos is not None:
        mask_lam = atmos
        mask_lam.extend(gaps)
        for i, mrange in enumerate(mask_lam):
            theseidxs = np.where((wave < mrange[0]) | (wave > mrange[1]))[0]
            goodpixels = np.intersect1d(theseidxs, goodpixels)
    goodpixels = np.unique(goodpixels)
    goodpixels.sort()

    if chk:
        plt.clf()
        plt.plot(newspec.wavelength[goodpixels], newspec.flux[goodpixels])
        plt.show()

    # Run it
    ppfit, miles, star_weights = fit_spectrum(
        newspec, zgal, R, degree_mult=0, degree_add=3,
        goodpixels=goodpixels, reddening=1., rebin=False)

    # Age
    age, metals = miles.mean_age_metal(star_weights)
    print('Age = {} Gyr'.format(age))
    print('Metals = {}'.format(metals))

    # Mass -- This is a bit approximate as Dwv is a guess for now
    actualflux = ppfit.bestfit * constants.L_sun.cgs / units.angstrom / (
            4 * np.pi * (cosmo.luminosity_distance(zgal).to(units.cm)) ** 2 / (1 + zgal))
    # When fitting, the routine thought our data and model spectra had same units...
    Dwv = 1700.  # Ang, width of the band pass
    scfactor = np.median(ppfit.bestfit * (units.erg / units.s / units.cm ** 2 / units.angstrom) / actualflux) * Dwv
    # To get the actual model mass required to fit spectrum, scale by this ratio
    massmodels = scfactor * miles.total_mass(star_weights)
    print('log10 M* = {}'.format(np.log10(massmodels)))

    # Reddening
    print('E(B-V) = {}'.format(ppfit.reddening))

    # Write?
    if results_file is not None:
        dump_ppxf_results(ppfit, miles, zgal, results_file)
    if spec_fit is not None:
        bestfit = dump_bestfit(ppfit, outfile=spec_fit, z=zgal)

    # Final check
    if chk:
        bestfit = dump_bestfit(ppfit, z=zgal)
        plt.clf()
        plt.plot(newspec.wavelength, newspec.flux)
        plt.plot(bestfit.wavelength, bestfit.flux)
        plt.show()


def fit_spectrum(spec, zgal, specresolution, tie_balmer=False,
                 miles_dir=None, rebin=True,
                 limit_doublets=False, degree_add=None,degree_mult=5,
                 **kwargs):
    """This is a wrapper for pPXF to fit stellar population models as well as
    emission lines to galaxy spectra.  Although pPXF allows otherwise, the
    emission lines are kinematically fixed to one another as are the stellar
    models, and the stars and gas independently vary with one another.

    Please see the pPXF documentation for more details on the vast array of
    parameters and options afforded by the software.

    The pPXF software may be downloaded at
    http://www-astro.physics.ox.ac.uk/~mxc/software/

    Parameters
    ----------
    spec : XSpectrum1D
        Spectrum to be fitted
    zgal : float
        Redshift of galaxy
    specresolution : float
        Spectral resolution (R) of the data
    tie_balmer : bool, optional
        Assume intrinsic Balmer decrement.  See documentation in ppxf_util.py,
        as this has implications for the derived reddening.
    limit_doublets : bool, optional
        Limit the ratios of [OII] and [SII] lines to ranges allowed by atomic
        physics.  See documentation in ppxf_util.py, as this has implications
        for getting the true flux values from those reported.
    degree_add : int, optional
        Degree of the additive Legendre polynomial used to modify the template
        continuum in the fit.
    degree_mult : int,optional
    miles_dir: str, optional
      Location of MILES models

    Returns
    -------
    ppfit : ppxf
        Object returned by pPXF; attributes are data pertaining to the fit
    miles : miles
        Contains information about the stellar templates used in the fit.
        See the documentation in miles_util.py for full details
    weights : 1d numpy vector
        Weights of the *stellar* template components.  Equivalent to the
        first N elements of ppfit.weights where N is the number of stellar
        templates used in the fit.
    """


    if spec is None:
        print('Galaxy has no spectrum associated')
        return None

    ### Spectra must be rebinned to a log scale (in wavelength).
    ### pPXF provides a routine to do this, but the input spectra
    ### must be sampled uniformly linear. linetools to the rescue!
    if rebin:
        meddiff = np.median(spec.wavelength.value[1:] - spec.wavelength.value[0:-1])
        newwave = np.arange(spec.wavelength.value[0], spec.wavelength.value[-1], meddiff)
        spec = spec.rebin(newwave * units.AA, do_sig=True, grow_bad_sig=True)

    # get data and transform wavelength to rest frame for template fitting
    wave = spec.wavelength.to(units.Angstrom).value
    flux = spec.flux.value
    flux = flux * (1. + zgal)
    noise = spec.sig.value
    noise = noise * (1. + zgal)
    wave = wave / (1. + zgal)
    # transform to air wavelengths for MILES templates and get approx. FWHM
    wave *= np.median(util.vac_to_air(wave) / wave)

    # use only wavelength range covered by templates
    mask = (wave > 3540) & (wave < 7409)
    #mask = (wave > 1682) & (wave < 10000.)
    maskidx = np.where(mask)[0]
    # also deal with declared good regions of the spectrum
    if 'goodpixels' in kwargs:
        goodpix = kwargs['goodpixels']
        pixmask = np.in1d(maskidx, goodpix)
        newgoodpix = np.where(pixmask)[0]
        kwargs['goodpixels'] = newgoodpix

    wave = wave[mask]
    flux = flux[mask]
    noise = noise[mask]

    # Nonnegative noise values are not allowed
    noise[noise <= 0] = np.max(noise)

    # pPXF requires the spectra to be log rebinned, so do it
    flux, logLam, velscale = util.log_rebin(np.array([wave[0], wave[-1]]), flux)
    noise, junk1, junk2 = util.log_rebin(np.array([wave[0], wave[-1]]), noise)

    ### The following lines unnecessary for DEIMOS/Hecto spectra due to their
    ### units but rescaling may be necessary for some
    # galaxy = flux / np.median(flux)  # Normalize spectrum to avoid numerical issues
    # print 'Scale flux by', round(np.median(flux),2)
    galaxy = flux  # use the native units

    # pPXF wants the spectral resolution in units of wavelength
    FWHM_gal = wave/ specresolution

    ### Set up stellar templates
    #miles_dir = resource_filename('ppxf', '/miles_models/')
    #miles_dir = resource_filename('ppxf', '/emiles_padova_chabrier/')
    if miles_dir is None:
        miles_dir = resource_filename('ppxf', '/miles_padova_chabrier/')
    #path4libcall = miles_dir + 'Mun1.30*.fits'
    #path4libcall = miles_dir + 'Ech1.30*.fits'
    path4libcall = miles_dir + 'Mch1.30*.fits'
    miles = lib.miles(path4libcall, velscale, FWHM_gal, wave_gal=wave)

    ### Stuff for regularization dimensions
    reg_dim = miles.templates.shape[1:]
    stars_templates = miles.templates.reshape(miles.templates.shape[0], -1)

    # See the pPXF documentation for the keyword REGUL
    regul_err = 0.01  # Desired regularization error

    ### Now the emission lines!  Only include lines in fit region.
    if 'goodpixels' in kwargs:
        gal_lam = wave[newgoodpix]
        # Also, log rebinning the spectrum change which pixels are 'goodpixels'
        newnewgoodpix = np.searchsorted(np.exp(logLam),gal_lam,side='left')
        uqnewnewgoodpix = np.unique(newnewgoodpix)
        if uqnewnewgoodpix[-1] == len(wave):
            uqnewnewgoodpix =uqnewnewgoodpix[:-1]
        kwargs['goodpixels'] = uqnewnewgoodpix

    else:
        gal_lam = wave

    def FWHM_func(wave):  # passed to generate emission line templates
        return wave / specresolution

    gas_templates, gas_names, line_wave = util.emission_lines_mask(
        miles.log_lam_temp, gal_lam, FWHM_func,
        tie_balmer=tie_balmer, limit_doublets=limit_doublets)

    # How many gas components do we have?
    balmerlines = [ll for ll in gas_names if ll[0] == 'H']
    numbalm = len(balmerlines)
    numforbid = len(gas_names) - numbalm

    # Stack all templates
    templates = np.column_stack([stars_templates, gas_templates])

    # other needed quantities
    dv = c * (miles.log_lam_temp[0] - logLam[0])
    ### use the following line if not transforming to z=0 first
    # vel = c * np.log(1 + zgal)  # eq.(8) of Cappellari (2017)
    vel = 0.  # We already transformed to the restframe!
    start = [vel, 25.]  # (km/s), starting guess for [V, sigma]

    ### Set up combination of templates
    n_temps = stars_templates.shape[1]
    n_balmer = 1 if tie_balmer else numbalm  # Number of Balmer lines included in the fit
    n_forbidden = numforbid  # Number of other lines included in the fit

    # Assign component=0 to the stellar templates, component=1 to the Balmer
    # emission lines templates, and component=2 to the forbidden lines.
    # component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
    component = [0] * n_temps + [1] * (n_balmer + n_forbidden)  # tie the gas lines together
    gas_component = np.array(component) > 0  # gas_component=True for gas templates

    # Fit (V, sig, h3, h4) moments=4 for the stars
    # and (V, sig) moments=2 for the two gas kinematic components
    if len(gas_names) > 0:
        moments = [2, 2]  # fix the gas kinematic components to one another
        start = [[vel, 50.], start]  # Adopt different gas/stars starting values
    else:
        moments = [2]  # only stars to be fit
        start = [vel, 50.]

        # If the Balmer lines are tied one should allow for gas reddening.
    # The gas_reddening can be different from the stellar one, if both are fitted.
    gas_reddening = 0 if tie_balmer else None


    if degree_add is None:
        degree_add = -1
    t = time.time()
    ppfit = ppxf.ppxf(templates, galaxy, noise, velscale, start,
                      plot=False, moments=moments, degree=degree_add, vsyst=dv,
                      lam=np.exp(logLam), clean=False, regul=1. / regul_err,
                      reg_dim=reg_dim,component=component, gas_component=gas_component,
                      gas_names=gas_names, gas_reddening=gas_reddening, mdegree=degree_mult,
                      **kwargs)

    print('Desired Delta Chi^2: %.4g' % np.sqrt(2 * galaxy.size))
    print('Current Delta Chi^2: %.4g' % ((ppfit.chi2 - 1) * galaxy.size))
    print('Elapsed time in PPXF: %.2f s' % (time.time() - t))

    weights = ppfit.weights[~gas_component]  # Exclude weights of the gas templates
    weights = weights.reshape(reg_dim)  # Normalized

    return ppfit, miles, weights

def total_mass(miles, weights, quiet=False):
    """
    Computes the total mass of living stars and stellar remnants
    in models fitted, given the weights produced and output by pPXF.

    A Salpeter IMF is assumed (slope=1.3) initially.
        -   TODO: Employ Chabrier models
    The returned mass excludes the gas lost during stellar evolution.

    This procedure uses the mass predictions
    from Vazdekis+12 and Ricciardelli+12
    http://adsabs.harvard.edu/abs/2012MNRAS.424..157V
    http://adsabs.harvard.edu/abs/2012MNRAS.424..172R
    they were downloaded in December 2016 below and are included in pPXF with permission
    http://www.iac.es/proyecto/miles/pages/photometric-predictions/based-on-miuscat-seds.php

    Parameters
    ----------
    miles : miles
        Miles object output from fit_spectrum()
    weights : 1d numpy array
        Weights vector corresponding to stellar templates.
        Output from fit_spectrum()
    quiet : bool, optional
        If True, do not print stellar mass result

    Returns
    -------
    mass_no_gas : float
        Total stellar mass of templates fitted.
        NOTE: this value will generally need to be scaled if there is any
        mismatch in units between data and models.  The MILES model spectra
        are generally given in units of L_sun/M_sun/Angstrom
    """
    from astropy.table import Table

    assert miles.age_grid.shape == miles.metal_grid.shape == weights.shape, \
        "Input weight dimensions do not match"

    #file_dir = path.dirname(path.realpath(__file__))  # path of this procedure
    #miles_dir = resource_filename('ppxf', '/miles_models/')
    miles_dir = resource_filename('ppxf', '/miles_padova_chabrier/')

    #file1 = miles_dir + "/Vazdekis2012_ssp_mass_Padova00_UN_baseFe_v10.0.txt"
    file1 = miles_dir + "out_mass_CH_PADOVA00.txt"

    colnames = ['IMF','slope','[M/H]','Age','Mtotal','M(*+remn)','M*','Mremn',
                'Mgas','M(*+remn)/Lv','M*/Lv','Mv','unknown']
    tab = Table.read(file1,format='ascii',names=colnames)

    slope1 = tab['slope']
    MH1 = tab['[M/H]']
    Age1 = tab['Age']
    m_no_gas = tab['Mtotal']

    # The following loop is a brute force but very safe and general
    # way of matching the photometric quantities to the SSP spectra.
    # It makes no assumption on the sorting and dimensions of the files
    mass_no_gas_grid = np.empty_like(weights)
    for j in range(miles.n_ages):
        for k in range(miles.n_metal):
            p1 = (np.abs(miles.age_grid[j, k] - Age1) < 0.001) & \
                 (np.abs(miles.metal_grid[j, k] - MH1) < 0.01) & \
                 (np.abs(1.30 - slope1) < 0.01)
            mass_no_gas_grid[j, k] = m_no_gas[p1]

    mass_no_gas = np.sum(weights*mass_no_gas_grid)

    if not quiet:
        print('Total mass: %.4g' % mass_no_gas)

    return mass_no_gas

def dump_bestfit(ppfit, outfile=None,  z=0.):
    """
    Create the bestfit in the observer frame and with vacuum wavelengths

    Parameters
    ----------
    ppfit
    outfile

    Returns
    -------
    bestfit: XSpectrum1D

    """
    meta = dict(airvac='air', headers=[None])
    # Spectrum
    bestfit = XSpectrum1D.from_tuple((ppfit.lam*(1+z), ppfit.bestfit/(1+z)), meta=meta)
    # Convert to vacuum
    bestfit.airtovac()
    # Write
    if outfile is not None:
        bestfit.write(outfile)
    # Return
    return bestfit


def dump_ppxf_results(ppfit, miles, z, outfile):
    """
    Write the stnadard results and the
    gas_component measurements to a simple ASCII file

    Parameters
    ----------
    ppfit: ppxf
    outfile: str

    Returns
    -------

    """
    # Get the lines (air)
    emission_lines, line_names, line_wave = util.emission_lines(
        np.array([0.1, 0.2]), [1000., 1e5], 0, limit_doublets=False, vacuum=True)
    # Construct a simple Table
    gas_tbl = Table()

    # Standard pPXF fit results
    meta = {}
    meta['EBV'] = ppfit.reddening

    star_weights = ppfit.weights[~ppfit.gas_component]
    star_weights = star_weights.reshape(ppfit.reg_dim)
    age, metals = miles.mean_age_metal(star_weights)
    meta['AGE'] = age
    meta['METALS'] = metals

    # Mass -- Approximate
    # Mass -- This is a bit approximate as Dwv is a guess for now
    actualflux = ppfit.bestfit * constants.L_sun.cgs / units.angstrom / (
            4 * np.pi * (cosmo.luminosity_distance(z).to(units.cm)) ** 2 / (1 + z))
    # When fitting, the routine thought our data and model spectra had same units...
    Dwv = 1700.  # Ang, width of the band pass
    scfactor = np.median(ppfit.bestfit * (units.erg / units.s / units.cm ** 2 / units.angstrom) / actualflux) * Dwv
    # To get the actual model mass required to fit spectrum, scale by this ratio
    massmodels = scfactor * miles.total_mass(star_weights)
    meta['LOGMSTAR'] = np.log10(massmodels.value)

    gas_tbl.meta = meta

    gas = ppfit.gas_component
    comp = ppfit.component[gas]
    gas_tbl['comp'] = comp
    gas_tbl['name'] = ppfit.gas_names
    gas_tbl['flux'] = ppfit.gas_flux
    gas_tbl['err'] = ppfit.gas_flux_error

    # Wavelengths
    waves = []
    for name in ppfit.gas_names:
        idx = np.where(line_names == name)[0][0]
        waves.append(line_wave[idx])
    gas_tbl['wave'] = waves

    vs = [ppfit.sol[icomp][0] for icomp in comp]
    sigs = [ppfit.sol[icomp][1] for icomp in comp]
    gas_tbl['v'] = vs
    gas_tbl['sig'] = sigs
    # Write
    gas_tbl.write(outfile, format='ascii.ecsv', overwrite=True)
    print("Wrote: {:s}".format(outfile))


class ppxfFit(object):
    def __init__(self,ppfit,miles,weights):
        """ Stripped down class of pPXf fit attributes to improve load times
        and data management.

        Parameters
        ----------
        ppfit : ppxf
            pPXF object output from fit_spectrum()
        miles : miles
            Miles object output from fit_spectrum()
        weights : array
            Weights on stellar pop models in fit; output from fit_spectrum()

        Returns
        -------
        mass_no_gas : float
            Total stellar mass of templates fitted.
            NOTE: this value will generally need to be scaled if there is any
            mismatch in units between data and models.  The MILES model spectra
            are generally given in units of L_sun/M_sun/Angstrom
        """

        self.galaxy = ppfit.galaxy
        self.nspec = ppfit.nspec  # nspec=2 for reflection-symmetric LOSVD
        self.npix = ppfit.npix  # total pixels in the galaxy spectrum
        self.noise = ppfit.noise
        self.clean = ppfit.clean
        self.fraction = ppfit.fraction
        self.gas_reddening = ppfit.gas_reddening
        self.degree = ppfit.degree
        self.mdegree = ppfit.mdegree
        self.method = ppfit.method
        self.quiet = ppfit.quiet
        self.sky = ppfit.sky
        self.vsyst = ppfit.vsyst
        self.regul = ppfit.regul
        self.lam = ppfit.lam
        self.nfev = ppfit.nfev
        self.reddening = ppfit.reddening
        self.reg_dim = ppfit.reg_dim
        self.reg_ord = ppfit.reg_ord
        self.star = ppfit.star
        self.npix_temp = ppfit.npix_temp
        self.ntemp = ppfit.ntemp
        self.factor = ppfit.factor
        self.sigma_diff = ppfit.sigma_diff
        self.status = ppfit.status  # Initialize status as failed
        self.velscale = ppfit.velscale

        self.bestfit = ppfit.bestfit
        self.chi2 = ppfit.chi2
        self.templates_rfft = ppfit.templates_rfft
        self.goodpixels = ppfit.goodpixels
        self.moments = ppfit.moments
        self.error = ppfit.error
        self.factor = ppfit.factor
        self.npad = ppfit.npad
        self.sol = ppfit.sol
        self.apoly = ppfit.apoly
        self.mpoly = ppfit.mpoly
        self.weights = ppfit.weights
        self.weights_ppxf = self.weights

        self.gas_component = ppfit.gas_component
        self.gas_bestfit = ppfit.gas_bestfit
        self.gas_flux = ppfit.gas_flux
        self.gas_flux_error = ppfit.gas_flux_error
        self.gas_names = ppfit.gas_names
        self.gas_mpoly = ppfit.gas_mpoly
        self.weights_ppxf = self.weights

        ### Now for MILES stuff
        self.age_grid = miles.age_grid
        self.metal_grid = miles.metal_grid
        self.n_ages = miles.n_ages
        self.n_metal = miles.n_metal

        self.mean_log_age, self.mean_metal = miles.mean_age_metal(weights, quiet=True)
        self.normfactor = miles.normfactor





