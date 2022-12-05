""" Code for calculations of P(DM|z) and P(z|DM)"""
import numpy as np
import os
from pkg_resources import resource_filename

from scipy.stats import norm, lognorm

from frb.dm import igm
from frb.dm import cosmic
from frb import defs

from IPython import embed

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

def grab_repo_grid():
    """
    Grab the grid from the Repository
    This may require the code to build it first!

    Returns:
        dict: Numpy dict from the npz save file
    """

    # File
    PDM_z_grid_file = os.path.join(
        resource_filename('frb', 'data'), 'DM',
        'PDM_z.npz')

    # Build?
    if not os.path.isfile(PDM_z_grid_file):
        build_grid_for_repo(PDM_z_grid_file)
            
    # Load
    print(f"Loading P(DM,z) grid from {PDM_z_grid_file}")
    sdict = np.load(PDM_z_grid_file)

    # Return
    return sdict