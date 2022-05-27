"""
A module to automate CIGALE. Currently works for a single galaxy.
It generates a configuration file and runs the standard pcigale
script. Requires pcigale already installed on the system. 
"""

import numpy as np
import sys, os, glob, multiprocessing, warnings
from collections import OrderedDict

from astropy.table import Table

try:
    from pcigale.session.configuration import Configuration
except ImportError:
    print("You will need to install pcigale to use the cigale.py module")
else:
    from pcigale.analysis_modules import get_module
    from pcigale.data import Database

from frb.surveys.catalog_utils import _detect_mag_cols, convert_mags_to_flux


from IPython import embed

# Default list of SED modules for CIGALE
_DEFAULT_SED_MODULES = ("sfhdelayed", "bc03", "nebular", "dustatt_calzleit", "dale2014",
                        "restframe_parameters", "redshifting")

#TODO Create a function to check the input filters
#Or create a translation file like eazy's.
#def check_filters(data_file):

def _sed_default_params(module):
    """
    Set the default parameters for CIGALE

    Args:
        module (str):
            Specify the SED using the CIGALE standard names, e.g. sfhdelayed, bc03, etc.

    Returns:
        params (dict): the default dict of SED modules
        and their initial parameters.
    """
    params = {}
    if module == "sfhdelayed":
        params['tau_main'] = (10**np.linspace(1,3,10)).tolist() #e-folding time of main population (Myr)
        params['age_main'] = (10**np.linspace(3,4,10)).tolist() #age (Myr)
        params['tau_burst'] = 50.0 #burst e-folding time (Myr)
        params['age_burst'] = 20.0
        params['f_burst'] = 0.0 #burst fraction by mass
        params['sfr_A'] = 0.1 #SFR at t = 0 (Msun/yr)
        params['normalise'] = False # Normalise SFH to produce one solar mass
    elif module == "bc03":
        params['imf'] = 1 #0: Salpeter 1: Chabrier
        params['metallicity'] = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05] 
        params['separation_age'] = 10 # Separation between yound and old stellar population (Myr)
    elif module == 'nebular':
        params['logU'] = -2.0 # Ionization parameter
        params['f_esc'] = 0.0 # Escape fraction of Ly continuum photons
        params['f_dust'] = 0.0 # Fraction of Ly continuum photons absorbed
        params['lines_width'] = 300.0
        params['emission'] = True
    elif module == 'dustatt_calzleit':
        params['E_BVs_young'] = [0.12, 0.25, 0.37, 0.5, 0.62, 0.74, 0.86] #Stellar color excess for young continuum
        params['E_BVs_old_factor'] = 1.0 # Reduction of E(B-V) for the old population w.r.t. young
        params['uv_bump_wavelength'] = 217.5 #central wavelength of UV bump (nm)
        params['uv_bump_width'] = 35.6 #UV bump FWHM (nm)
        params['uv_bump_amplitude'] = 1.3 # Amplitude of the UV bump. For the Milky Way: 3.
        # The following parameter can have a significant affect on stellar mass
        #  We use the recommendation in Lo Faro+2017
        params['powerlaw_slope'] = -0.13  # Slope delta of the power law modifying the attenuation curve.
        # These filters have no effect
        params['filters'] = 'B_B90 & V_B90 & FUV'
    elif module == 'dale2014':
        params['fracAGN'] = [0.0,0.05,0.1,0.2]
        params['alpha'] = 2.0
    elif module == 'restframe_parameters':
        params['beta_calz94'] = False
        params['D4000'] = False
        params['IRX'] = False
        params['EW_lines'] = '500.7/1.0 & 656.3/1.0'
        params['luminosity_filters'] = 'u_prime & r_prime'
        params['colours_filters'] = 'u_prime-r_prime'
    elif module == 'redshifting':
        params['redshift'] = '' #Use input redshifts
    return params


def gen_cigale_in(photometry_table, zcol, idcol=None, infile="cigale_in.fits",
                  overwrite=True, **kwargs):
    """
    Generates the input catalog from
    a photometric catalog.

    Args:
        photometry_table (astropy Table):
            A table from some photometric
            catalog with magnitudes and
            error measurements. Currently supports
            DES, DECaLS, SDSS, Pan-STARRS and WISE

            The naming convention follows those specified in frb.galaxies.defs
            with the exception of WISE which use WISE-1, etc. although the code
            also handles WISE-W1, etc.
        zcol (str):
            Name of the column with redshift estimates
        idcol (str, optional):
            Name of the column with object IDs. By default,
            the code looks for the first column with "ID" in
            its name. If that's not present, it creates a
            column with row numbers for IDs.
        infile (str, optional):
            Output name + path for the CIGALE input file generated
        overwrite (bool, optional):
            If true, overwrites file if it already exists
        kwargs: only here to catch extras
    """
    #Table must have a column with redshift estimates
    if not isinstance(zcol, str):
        raise IOError("zcol must be a column name. i.e. a string")
    assert zcol in photometry_table.colnames, "{} not found in the table. Please check".format(zcol)

    magcols, mag_errcols = _detect_mag_cols(photometry_table)
    cigtab = photometry_table.copy()
    cigtab.rename_column(zcol,"redshift")
    photom_cols = magcols+mag_errcols

    # Rename any column with "ID" in it to "id"
    if idcol is None:
        try:
            idcol = [col for col in cigtab.colnames if "ID" in col.upper()][0]
        except IndexError:
            print("No column with 'ID' in name. Adding a column.")
            idcol = 'id'
            cigtab[idcol] = np.arange(len(cigtab))+1 

    cigtab.rename_column(idcol,"id")
    
    #First round of renaming
    cigtab = convert_mags_to_flux(cigtab)
    cigtab = cigtab[['id','redshift']+photom_cols]

    # Rename our filters to CIGALE names, as needed
    new_names = {
        'SDSS_u': 'sdss.up',
        'SDSS_g': 'sdss.gp',
        'SDSS_r': 'sdss.rp',
        'SDSS_i': 'sdss.ip',
        'SDSS_z': 'sdss.zp',
        'VLT_u': 'VLT_FORS2_u',
        'VLT_g': 'VLT_FORS2_g',
        'VLT_I': 'VLT_FORS2_I',
        'VLT_z': 'VLT_FORS2_z',
        'WISE_W1': 'WISE1',
        'WISE_W2': 'WISE2',
        'WISE_W3': 'WISE3',
        'WISE_W4': 'WISE4',
        'VISTA_Y': 'vista.vircam.Y',
        'VISTA_J': 'vista.vircam.J',
        'VISTA_H': 'vista.vircam.H',
        'VISTA_Ks': 'vista.vircam.Ks',
        'LRISr_I': 'LRIS_I',
        'LRISb_V': 'LRIS_V',
        'WFC3_F160W': 'hst.wfc3.F160W',
        'WFC3_F300X': 'WFC3_F300X',
        'Spitzer_3.6': 'spitzer.irac.ch1',
        'Spitzer_4.5': 'spitzer.irac.ch2',
        'NSC_g': 'DES_g',
        'NSC_r': 'DES_r',
        'NSC_i': 'DES_i',
        'NSC_z': 'DES_z',
        'NSC_Y': 'DES_Y',
        'DECam_r': 'DES_r'

    }
    for key in new_names:
        if key in photom_cols:
            cigtab.rename_column(key, new_names[key])
            # Try Error
            if key+'_err' in photom_cols:
                cigtab.rename_column(key+'_err', new_names[key]+'_err')

    cigtab.write(infile,overwrite=overwrite)
    return


def _initialise(data_file, config_file="pcigale.ini",
                cores=None, save_sed=False, variables="", sed_modules=_DEFAULT_SED_MODULES,
                sed_modules_params=None, **kwargs):
    """
    Initialise a CIGALE configuration file and write to disk.
    
    Args:
        data_file (str):
            Path to the input photometry data file.
        config_file (str, optional):
            Path to the file where CIGALE's configuration
            is stored.
        cores (int, optional):
            Number of CPU cores to be used. Defaults
            to all cores on the system.
        save_sed (bool, optional):
            Save the best fit SEDs to disk for each galaxy.
        variables (str or list, optional):
            A single galaxy property name to save to results
            or a list of variable names. Names must belong
            to the list defined in the CIGALE documentation.
        sed_modules (list or tuple, optional):
            A list of SED modules to be used in the 
            PDF analysis. If this is being input, there
            should be a corresponding correct dict
            for sed_modules_params.
        sed_module_params (dict, optional):
            A dict containing parameter values for
            the input SED modules. Better not use this
            unless you know exactly what you're doing.
        kwargs: only here to catch extras
    Returns:
        cigconf (pcigale.session.configuration.Configuration):
                CIGALE Configuration object
    """
    # Check
    if sed_modules !=_DEFAULT_SED_MODULES:
        assert sed_modules_params is not None,\
             "If you're not using the default modules, you'll have to input SED parameters"
    # Init
    cigconf = Configuration(config_file) #a set of dicts, mostly
    cigconf.create_blank_conf() #Initialises a pcigale.ini file

    # fill in initial values
    cigconf.pcigaleini_exists = True
    cigconf.config['data_file'] = data_file
    cigconf.config['param_file'] = ""
    cigconf.config['sed_modules'] = sed_modules
    cigconf.config['analysis_method'] = 'pdf_analysis'
    if cores is None:
        cores = multiprocessing.cpu_count() #Use all cores
    cigconf.config['cores'] = cores
    cigconf.generate_conf() #Writes defaults to config_file
    cigconf.config['analysis_params']['variables'] = variables
    cigconf.config['analysis_params']['save_best_sed'] = save_sed
    cigconf.config['analysis_params']['lim_flag'] = True

    # Change the default values to new defaults:
    if sed_modules_params is None:
        sed_modules_params = {}
        for module in sed_modules:
            sed_modules_params[module] = _sed_default_params(module)
    cigconf.config['sed_modules_params'] = sed_modules_params

    # Overwrites the config file
    cigconf.config.write()

    # Return
    return cigconf


def run(photometry_table, zcol, data_file="cigale_in.fits", config_file="pcigale.ini",
        wait_for_input=False, plot=True, outdir='out', compare_obs_model=False, **kwargs):
    """
    Input parameters and then run CIGALE.

    Args:
        photometry_table (astropy Table):
            A table from some photometric catalog with magnitudes and
            error measurements. Currently supports
            DES, DECaLS, SDSS, Pan-STARRS and WISE
        zcol (str):
            Name of the column with redshift estimates.
        data_file (str, optional):
            Root name for the photometry data file generated used as input to CIGALE
        config_file (str, optional):
            Root name for the file where CIGALE's configuration is generated
        wait_for_input (bool, optional):
            If true, waits for the user to finish editing the auto-generated config file
            before running.
        plot (bool, optional):
            Plots the best fit SED if true
        cores (int, optional):
            Number of CPU cores to be used. Defaults
            to all cores on the system.
        outdir (str, optional):
            Path to the many outputs of CIGALE
            If not supplied, the outputs will appear in a folder named out/
        compare_obs_model (bool, optional):
            If True compare the input observed fluxes with the model fluxes
            This writes a Table to outdir named 'photo_observed_model.dat'

    kwargs:  These are passed into gen_cigale_in() and _initialise()
        save_sed (bool, optional):
            Save the best fit SEDs to disk for each galaxy.
        variables (str or list, optional):
            A single galaxy property name to save to results
            or a list of variable names. Names must belong
            to the list defined in the CIGALE documentation.
        sed_modules (list of 'str', optional):
            A list of SED modules to be used in the 
            PDF analysis. If this is being input, there
            should be a corresponding correct dict
            for sed_modules_params.
        sed_module_params (dict, optional):
            A dict containing parameter values for
            the input SED modules. Better not use this
            unless you know exactly what you're doing.

    """
    gen_cigale_in(photometry_table,zcol,infile=data_file,overwrite=True, **kwargs)
    _initialise(data_file, config_file=config_file,**kwargs)
    if wait_for_input:
        input("Edit the generated config file {:s} and press any key to run.".format(config_file))
    cigconf = Configuration(config_file)
    analysis_module = get_module(cigconf.configuration['analysis_method'])
    analysis_module.process(cigconf.configuration)
    if plot:
        try:
            from pcigale_plots import sed  # This modifies the backend to Agg so I hide it here
            old_version = True
        except ImportError:
            from pcigale_plots.plot_types.sed import sed
            old_version = False
        
        if old_version:
            import pcigale
            #warnings.warn("You are using CIGALE version {:s}, for which support is deprecated. Please update to 2020.0 or higher.".format(pcigale.__version__))
            sed(cigconf,"mJy",True)
        else:
            # TODO: Let the user customize the plot.
            series = ['stellar_attenuated', 'stellar_unattenuated', 'dust', 'agn', 'model']
            sed(cigconf,"mJy",True, (False, False), (False, False), series, "pdf", "out")
        # Set back to a GUI
        import matplotlib
        matplotlib.use('TkAgg')

    # Rename the default output directory?
    if outdir != 'out':
        try:
            os.system("rm -rf {}".format(outdir))
            os.system("mv out {:s}".format(outdir))
        except:
            print("Invalid output directory path. Output stored in out/")

    # Move input files into outdir too
    os.system("mv {:s} {:s}".format(data_file, outdir))
    os.system("mv {:s} {:s}".format(config_file, outdir))
    os.system("mv {:s}.spec {:s}".format(config_file, outdir))

    # Compare?
    if compare_obs_model:
        #Generate an observation/model flux comparison table.
        with Database() as base:
            filters = OrderedDict([(name, base.get_filter(name))
                                for name in cigconf.configuration['bands']
                                if not (name.endswith('_err') or name.startswith('line')) ])
            filters_wl = np.array([filt.pivot_wavelength
                                    for filt in filters.values()])
            mods = Table.read(outdir+'/results.fits')

            try:
                obs = Table.read(os.path.join(outdir, cigconf.configuration['data_file']))
            except:
                print("Something went wrong here. Astropy was unable to read the observations table. Please ensure it is in the fits format.")
                return
            for model, obj in zip(mods, obs):
                photo_obs_model = Table()
                photo_obs_model['lambda_filter'] = [wl/1000 for wl in filters_wl]
                photo_obs_model['model_flux'] = np.array([model["best."+filt] for filt in filters.keys()])
                photo_obs_model['observed_flux'] = np.array([obj[filt] for filt in filters.keys()])
                photo_obs_model['observed_flux_err'] = np.array([obj[filt+'_err'] for filt in filters.keys()])
                photo_obs_model.write(outdir+"/photo_observed_model_"+str(model['id'])+".dat",format="ascii",overwrite=True)
            #import pdb; pdb.set_trace()
            
    return

def host_run(host, cut_photom=None, cigale_file=None):
    """
    Run CIGALE on an FRBGalaxy's photometry
    and store results in a folder with the
    FRBGalaxy's name.
    Args
    ----
    photom (astropy Table): Table containing
        galaxy photometry. Table columns
        must be in the format '<SOURCE>_<BAND>'
        and '<SOURCE>_<BAND>_err'.
        e.g. SDSS_u, SDSS_u_err, Pan-STARRS_g
    host (FRBGalaxy): A host galaxy.
    cigale_file (str, optional): Name of main
        CIGALE output file. Must be in the format
        `<something>_CIGALE.fits`. No file is
        renamed if nothing is provided.
    """
    cigale_tbl = Table()
    cigale_tbl['z'] = [host.z]
    cigale_tbl['ID'] = host.name

    # Deal with photometry
    if cut_photom is not None:
        photom_obj = cut_photom
    else:
        photom_obj = host.photom
    for key in photom_obj.keys():
        cigale_tbl[key] = photom_obj[key]

    # Run
    run(cigale_tbl, 'z', outdir=host.name, compare_obs_model=True, idcol='ID')

    # Rename/move
    if cigale_file is not None:
        os.system('cp -rp {:s}/results.fits {:s}'.format(host.name, cigale_file))
        model_file = cigale_file.replace('CIGALE', 'CIGALE_model')
        os.system('cp -rp {:s}/{:s}_best_model.fits {:s}'.format(host.name, host.name, model_file))
        photo_file = cigale_file.replace('CIGALE.fits', 'CIGALE_photo.dat')
        os.system('cp -rp {:s}/photo_observed_model_{:s}.dat {:s}'.format(host.name, host.name, photo_file))
        # SFH
        sfh_file = cigale_file.replace('CIGALE', 'CIGALE_SFH')
        os.system('mv {:s}/{:s}_SFH.fits {:s}'.format(host.name, host.name, sfh_file))
    return
