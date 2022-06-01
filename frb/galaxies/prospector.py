""" WARNING -- WIP -- Methods related to running and reading Prospector outputs """
import os, sys
import pickle
import json

import numpy as np
import h5py

from frb.galaxies import nebular
from frb.galaxies import utils as galaxy_utils
from frb.galaxies import defs

# TODO -- THIS NEEDS MORE FILTERS
prospector_filter_dict = {}
for key in defs.SDSS_bands:
    prospector_filter_dict['SDSS_{:s}'.format(key)] = 'sdss_{:s}'.format(key)

def load_spec_from_specdb(host_obj, INSTR:str, deredden=True,
                          ebv_method = 'SandF'):
    # Grab it
    meta, spectrum = host_obj.get_metaspec(instr=INSTR)

    # De-redden?
    if deredden:
        EBV = nebular.get_ebv(host_obj.coord, definition=ebv_method)['meanValue'] 
        spectrum = galaxy_utils.deredden_spec(spectrum, EBV)
    
    # Return
    return spectrum.wavelength.values, spectrum.flux.values, spectrum.sig.values


def grab_magnitudes(host_obj, FRB_filters:list):
    # Load up magnitudes from a host object

    prospector_filters, m_AB, m_AB_err = [], [], []

    for ifilt in FRB_filters:
        # Vette
        assert ifilt in defs.valid_filters
        # Translate
        prospector_filters.append(prospector_filter_dict[ifilt])
        # Grab data
        m_AB.append(host_obj.photom[ifilt])
        m_AB_err.append(host_obj.photom[ifilt+'_err'])
    
    # Return
    return prospector_filters, m_AB, m_AB_err
    
    
# #######################################################
# CODE TAKEN FROM PROSPECTOR
# #######################################################


def results_from(filename, model_file=None, dangerous=True, **kwargs):
    """Read a results file with stored model and MCMC chains.
    :param filename:
        Name and path to the file holding the results.  If ``filename`` ends in
        "h5" then it is assumed that this is an HDF5 file, otherwise it is
        assumed to be a pickle.
    :param dangerous: (default, True)
        If True, use the stored paramfile text to import the parameter file and
        reconstitute the model object.  This executes code in the stored
        paramfile text during import, and is therefore dangerous.
    :returns results:
        A dictionary of various results including:
          + `"chain"`  - Samples from the posterior probability (ndarray).
          + `"lnprobability"` - The posterior probability of each sample.
          + `"weights"` -  The weight of each sample, if `dynesty` was used.
          + `"theta_labels"` - List of strings describing free parameters.
          + `"bestfit"` - The prediction of the data for the posterior sample with
            the highest `"lnprobability"`, as a dictionary.
          + `"run_params"` - A dictionary of arguments supplied to prospector at
            the time of the fit.
          + `"paramfile_text"` - Text of the file used to run prospector, string
    :returns obs:
        The obs dictionary
    :returns model:
        The models.SedModel() object, if it could be regenerated from the stored
        `"paramfile_text"`.  Otherwise, `None`.
    """
    # Read the basic chain, parameter, and run_params info
    if filename.split('.')[-1] == 'h5':
        res = read_hdf5(filename, **kwargs)
        if "_mcmc.h5" in filename:
            mf_default = filename.replace('_mcmc.h5', '_model')
        else:
            mf_default = "x"
    else:
        with open(filename, 'rb') as rf:
            res = pickle.load(rf)
        mf_default = filename.replace('_mcmc', '_model')

    # Now try to read the model object itself from a pickle
    if model_file is None:
        mname = mf_default
    else:
        mname = model_file
    param_file = (res['run_params'].get('param_file', ''),
                  res.get("paramfile_text", ''))
    model, powell_results = read_model(mname, param_file=param_file,
                                       dangerous=dangerous, **kwargs)
    #if dangerous:
    #    try:
    #        model = get_model(res)
    #    except:
    #        model = None
    res['model'] = model
    if powell_results is not None:
        res["powell_results"] = powell_results

    return res, res["obs"], model


def read_model(model_file, param_file=('', ''), dangerous=False, **extras):
    """Read the model pickle.  This can be difficult if there are user defined
    functions that have to be loaded dynamically.  In that case, import the
    string version of the paramfile and *then* try to unpickle the model
    object.
    :param model_file:
        String, name and path to the model pickle.
    :param dangerous: (default: False)
        If True, try to import the given paramfile.
    :param param_file:
        2-element tuple.  The first element is the name of the paramfile, which
        will be used to set the name of the imported module.  The second
        element is the param_file contents as a string.  The code in this
        string will be imported.
    """
    model = powell_results = None
    if os.path.exists(model_file):
        try:
            with open(model_file, 'rb') as mf:
                mod = pickle.load(mf)
        except(AttributeError):
            # Here one can deal with module and class names that changed
            with open(model_file, 'rb') as mf:
                mod = load(mf)
        except(ImportError, KeyError):
            # here we load the parameter file as a module using the stored
            # source string.  Obviously this is dangerous as it will execute
            # whatever is in the stored source string.  But it can be used to
            # recover functions (especially dependcy functions) that are user
            # defined
            path, filename = os.path.split(param_file[0])
            modname = filename.replace('.py', '')
            if dangerous:
                user_module = import_module_from_string(param_file[1], modname)
            with open(model_file, 'rb') as mf:
                mod = pickle.load(mf)

        model = mod['model']

        for k, v in list(model.theta_index.items()):
            if type(v) is tuple:
                model.theta_index[k] = slice(*v)
        powell_results = mod['powell']

    return model, powell_results

def import_module_from_string(source, name, add_to_sys_modules=True):
    """Well this seems dangerous.
    """
    import imp
    user_module = imp.new_module(name)
    exec(source, user_module.__dict__)
    if add_to_sys_modules:
        sys.modules[name] = user_module

    return user_module



def read_hdf5(filename, **extras):
    """Read an HDF5 file (with a specific format) into a dictionary of results.
    This HDF5 file is assumed to have the groups ``sampling`` and ``obs`` which
    respectively contain the sampling chain and the observational data used in
    the inference.
    All attributes of these groups as well as top-level attributes are loaded
    into the top-level of the dictionary using ``json.loads``, and therefore
    must have been written with ``json.dumps``.  This should probably use
    JSONDecoders, but who has time to learn that.
    :param filename:
        Name of the HDF5 file.
    """
    groups = {"sampling": {}, "obs": {},
              "bestfit": {}, "optimization": {}}
    res = {}
    with h5py.File(filename, "r") as hf:
        # loop over the groups
        for group, d in groups.items():
            # check the group exists
            if group not in hf:
                continue
            # read the arrays in that group into the dictionary for that group
            for k, v in hf[group].items():
                d[k] = np.array(v)
            # unserialize the attributes and put them in the dictionary
            for k, v in hf[group].attrs.items():
                try:
                    d[k] = json.loads(v)
                except:
                    try:
                        d[k] = unpick(v)
                    except:
                        d[k] = v
        # do top-level attributes.
        for k, v in hf.attrs.items():
            try:
                res[k] = json.loads(v)
            except:
                try:
                    res[k] = unpick(v)
                except:
                    res[k] = v
        res.update(groups['sampling'])
        res["bestfit"] = groups["bestfit"]
        res["optimization"] = groups["optimization"]
        res['obs'] = groups['obs']
        try:
            res['obs']['filters'] = load_filters([str(f) for f in res['obs']['filters']])
        except:
            pass
        try:
            res['rstate'] = unpick(res['rstate'])
        except:
            pass
        #try:
        #    mp = [names_to_functions(p.copy()) for p in res['model_params']]
        #    res['model_params'] = mp
        #except:
        #    pass

    return res


def sample_posterior(chain, weights=None, nsample=int(1e4),
                     start=0, thin=1, extra=None):
    """
    :param chain:
        ndarray of shape (niter, ndim) or (niter, nwalker, ndim)
    :param weights:
        weights for each sample, of shape (niter,)
    :param nsample: (optional, default: 10000)
        Number of samples to take
    :param start: (optional, default: 0.)
        Fraction of the beginning of the chain to throw away, expressed as a float in the range [0,1]
    :param thin: (optional, default: 1.)
        Thinning to apply to the chain before sampling (why would you do that?)
    :param extra: (optional, default: None)
        Array of extra values to sample along with the parameters of the chain.
        ndarray of shape (niter, ...)
    """
    start_index = np.floor(start * (chain.shape[-2] - 1)).astype(int)
    if chain.ndim > 2:
        flatchain = chain[:, start_index::thin, :]
        nwalker, niter, ndim = flatchain.shape
        flatchain = flatchain.reshape(niter * nwalker, ndim)
    elif chain.ndim == 2:
        flatchain = chain[start_index::thin, :]
        niter, ndim = flatchain.shape

    if weights is not None:
        p = weights[start_index::thin]
        p /= p.sum()
    else:
        p = None

    inds = np.random.choice(niter, size=nsample, p=p)
    if extra is None:
        return flatchain[inds, :]
    else:
        return flatchain[inds, :], extra[inds, ...]
