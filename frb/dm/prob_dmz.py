""" Code for calculations of P(DM|z) and P(z|DM)"""
import numpy as np
import os


from importlib.resources import files

from scipy.stats import norm, lognorm
from scipy.interpolate import interp1d, interp2d

from frb.dm import igm
from frb.dm import cosmic
from frb import defs

from IPython import embed

# Dict to resolve telescope/survey into the appropriate grid filename
telescope_dict = {
    'CHIME_repeaters': 'CHIME_pzdm_repeaters.npz',
    'CHIME': 'CHIME_pzdm.npz',
    'DSA': 'DSA_pzdm.npz',
    'Parkes': 'parkes_mb_class_I_and_II_pzdm.npz',
    'CRAFT': 'CRAFT_class_I_and_II_pzdm.npz',
    'CRAFT_ICS_1300': 'CRAFT_ICS_1300_pzdm.npz',
    'CRAFT_ICS_892': 'CRAFT_ICS_892_pzdm.npz',
    'CRAFT_ICS_1632': 'CRAFT_ICS_1632_pzdm.npz',
    'FAST': 'FAST_pzdm.npz',
    'perfect': 'PDM_z.npz'
}

class P_DM_z(object):
    pass


def prob_DMcosmic_FRB(frb, DM_min=0., DM_max=5000., step=1.,
                      ISMfrac=0.10, DM_MWhalo=50.):
    """
    Generate P(DM_cosmic) for an input FRP
    
    Args:
        frb (:class:`frb.frb.FRB`):
        DM_min (float, optional):
            Lowest DM for the calulation
        DM_max (float, optional):
            Highest DM for the calulation
        step (float, optional):
            Step size of DM array in units of pc/cm**3
        ISMfrac (float, optional):
            Fraction of DM_ISM to adopt as the 1-sigma error
        DM_MWhalo (float, optional):
            Fixed value to use for the MW halo

    Returns:
        tuple:  numpy.ndarray, numpy.ndarray
            DM_cosmic values (units of pc/cm**3), P(DM_cosmic) normalized to unity

    """
    # Init
    DMcosmics = np.arange(DM_min, DM_max+step, step)
    P_DM_cosmic = np.zeros_like(DMcosmics)

    # ISM
    scale = np.pi * ISMfrac * frb.DMISM.value
    p_ISM = norm(loc=frb.DMISM.value, scale=scale)

    # Pre calculate
    DM_ISMs = DMcosmics
    pdf_ISM = p_ISM.pdf(DM_ISMs)

    # Host 
    # TODO Should use the MCMC chains to do this right!
    #  And should fix Omega_b true
    exp_u = 68.2  # Median
    sigma_host = 0.88  # Median
    lognorm_floor=0.
    p_host = lognorm(s=sigma_host, loc=lognorm_floor, scale=exp_u)

    # Loop time
    for kk, DMcosmic in enumerate(DMcosmics):
        DM_host = frb.DM.value - DM_MWhalo - DM_ISMs - DMcosmic
        # Prob time
        Prob = pdf_ISM * p_host.pdf(DM_host*(1+frb.z))
        # Sum
        P_DM_cosmic[kk] = np.sum(Prob)

    # Normalize
    P_DM_cosmic = P_DM_cosmic / np.sum(P_DM_cosmic)

    # Return
    return DMcosmics, P_DM_cosmic


def grid_P_DMcosmic_z(beta=3., F=0.31, zvals=None, 
                      DM_cosmics=None,
                      cosmo=defs.frb_cosmo):
    """
    Generate a grid of P(DM_cosmic|z)

    Args:
        beta (float, optional):
            sigma_DM_cosmic parameter
        F (float, optional):
            Feedback parameter (higher F means weaker feedback)
        zvals (np.ndarray, optional):
            Redshifts for the grid
        DMcosmic (np.ndarray, optional):
            DMs for the grid
        cosmo (optional):
            Cosmology

    Returns:
        tuple: z, DM_cosmic, P(DM_cosmic|z)
    """
    # Check
    if not np.isclose(beta, 3.):
        raise IOError("Not prepared for this beta value (yet)")
    # Load
    # sigma_DM
    f_C0_3 = cosmic.grab_C0_spline()

    # Grid
    if zvals is None:
        zvals = np.linspace(0., 2., 200)
    if DM_cosmics is None:
        DM_cosmics = np.linspace(1., 5000., 1000)

    PDF_grid = np.zeros((DM_cosmics.size, zvals.size))

    # Loop
    for kk, zval in enumerate(zvals):
        # z=0
        if zval == 0:
            PDF_grid[0,0] = 1.
            continue
        avgDM = igm.average_DM(zval, cosmo=cosmo).value
        # Params
        sigma = F / np.sqrt(zval)
        C0 = f_C0_3(sigma)
        #  Delta
        Delta = DM_cosmics / avgDM
        # PDF time
        PDF = cosmic.DMcosmic_PDF(Delta, C0, sigma)
        # Normalize
        PDF_grid[:,kk] = PDF / np.sum(PDF)

    # Return
    return zvals, DM_cosmics, PDF_grid


def build_grid_for_repo(outfile:str):
    """
    Build a P(DM,z) grid for the Repository

    Args:
        outfile (str): Path+filename for output file
    """

    print("Generating a new PDM_z grid for the Repo")
    print("Please be patient (will take a few minutes)....")
    #
    zvals = np.linspace(0., 4., 200)
    z, DM, P_DM_z = grid_P_DMcosmic_z(zvals=zvals)
    # Write
    np.savez(outfile, z=z, DM=DM, PDM_z=P_DM_z)
    print(f"File written: {outfile}")
    print("This will be used going forth")



def grab_repo_grid(grid_name):
    """
    Grab the grid from the Repository based on the given grid name

    Args:
        grid_name (str): Name of the grid to grab

    Returns:
        dict: Numpy dict from the npz or npy save file
    """

    # File
    grid_file = files('frb.data.DM').joinpath(grid_name)
    
    # Build?
    if grid_name == 'PDM_z.npz':
        if not os.path.isfile(grid_file):
            build_grid_for_repo(grid_file)

    # Load
    print(f"Loading P(DM,z) grid from {grid_file}")
    sdict = np.load(grid_file)

    # Return
    return sdict



def get_DMcosmic_from_z(zarray, perc=0.5, redo_pdmz_grid=False, DMevals=np.linspace(1.,2000.,1000), beta=3., F=0.31, cosmo=defs.frb_cosmo):
    """
    Gives DMcosmic values of zarray, considering the percentile.
    Default is median (i.e. percentile=0.5)

    Parameters
    ----------
    zarray : float or np.array
        Redshift value or values
        If np.array, then it must start from 0
    perc : float
        Percentile of the PDF(DM_cosmic) for each z, from 0 to zmax. 
        Default = 0.5
    redo_pdmz_grid : bool
        Whether to recompute the DMcosmic-z grid
        Default is False
        If True, it will be computational expensive
    DMevals : np.array
        Array representing the DMcosmic evaluations of the PDF. Should start with 1.
        Default = np.linspace(1.,5000.,1000)
        Changing this only works when redo_pdmz_grid=True
    beta : float
        Slope of a galactic halo parameter (See Macquart+2020)
        Default = 3.0 
        Changing this only works when redo_pdmz_grid=True
    F : float
        Parameter that depends on galaxy feedback (See Macquart+2020)
        Default = 0.31
        Changing this only works when redo_pdmz_grid=True
    cosmo : astropy.cosmology object
        Cosmology


    Returns
    -------
    zarray : np.array
        z values where the DMcosmic was computed
    DMcosmic_array : np.array
        DMcosmic values corresponding to the z_array (in a given percentile)
    
    """
    
    #check z input
    if np.isscalar(zarray):
        z_new = np.array([zarray])
    else:
        z_new = np.array(zarray)
   
    # Get the DMcosmic-z grid (will read a default one if redo=False)
    filename = resource_filename('frb', os.path.join('data','DM', 'pdmz_default_grid.npz'))  
    if redo_pdmz_grid:
        # This is time consuming, could be done more efficiently
        print("Calculating the DMcosmic-z grid, this may take a while... [you can turn off this by setting `redo_pdmz_grid=False`]")
        zeval , DM_cosmics, PDF_grid = grid_P_DMcosmic_z(zvals=z_new, beta=beta, F=F, cosmo=cosmo, DM_cosmics=DMevals)
        grid = PDF_grid
    else:
        if not os.path.isfile(filename):
            # here we hardcode the default values for the grid to be those of grid_P_DMcosmic_z() defaults
            print("No {} found in the repo. Calculating it. This may take a while, but should only happen 1 time.")
            zeval , DM_cosmics, PDF_grid = grid_P_DMcosmic_z()
            np.savez_compressed(filename, zeval=zeval, DM_cosmics=DM_cosmics, PDF_grid=PDF_grid)
        else:
            data = np.load(filename)
            zeval = data['zeval']
            DM_cosmics = data['DM_cosmics']
            PDF_grid = data['PDF_grid']

        # make a interpolated version of the user's grid
        f = interp2d(zeval, DM_cosmics, PDF_grid, kind='cubic')
        grid = f(zarray, DM_cosmics)

   
    # Get the relation at a given percentile
    perc_macquart = []  
    for ii, column in enumerate(grid.T):
        # import pdb; pdb.set_trace()
        # this could be speed up using a finer sampling for the interpolated grid
        if redo_pdmz_grid is not True and np.min(zarray)<0.02:
            # TODO: solve a bug when values are too low for redshift... seems a problems with the 1d interpolation.
            raise NotImplementedError('Not yet fully implemented for this low-z regime given the coarse grid. Please use zmin>0.02')
        dm_pdf_interp = interp1d(np.cumsum(column.flatten()), DM_cosmics, fill_value=0., bounds_error=False)
        perc_macquart.append(dm_pdf_interp(perc))
    perc_macquart = np.array(perc_macquart).flatten()
    
    # reformat output
    return z_new, perc_macquart
