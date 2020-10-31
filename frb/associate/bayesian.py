"""Methods related to Bayesian association analysis"""

import numpy as np
import os
from pkg_resources import resource_filename

import copy

from astropy import units
from astropy.coordinates import SkyCoord

from frb.galaxies import hosts

from IPython import embed

# Globals -- to speed up calculations
r_dat, mag_uniq, _ = hosts.read_r_mags(
    resource_filename('frb', os.path.join('data', 'Galaxies', 'driver2016_t3data.fits')))
eb17_spl = hosts.interpolate.UnivariateSpline(x=mag_uniq,
                                   y=np.log10(r_dat),
                                   bbox=[-100, 100],
                                   k=3)
def n_gal(m_r):
    return 10 ** eb17_spl(m_r)


def add_contam(ncontam, frb_coord, cand_gal, cand_coord, fov,
               rstate=None, rmag=25.):
    if ncontam != 1:
        import pdb; pdb.set_trace()  # Not ready for more than 1
    if rstate is None:
        rstate = np.random.RandomState()
    # Generate a random contaminant -- Square inside a circle to keep it simple
    s = 2*fov / np.sqrt(2)
    x,y  = rstate.uniform(size=2*ncontam, low=0, high=s)
    r = np.sqrt(x**2 + y**2)
    pa = float(rstate.uniform(size=ncontam, low=0., high=360.))
    #
    c_coord = frb_coord.directional_offset_by(pa*units.deg, r*units.arcsec)
    # Add to table
    row = copy.deepcopy(cand_gal[0])
    row['DES_r'] = rmag
    row['ra'] = c_coord.ra.value
    row['dec'] = c_coord.dec.value
    cand_gal.add_row(row)
    # Coordinate
    cand_coord = SkyCoord([cand_coord, c_coord])
    # Return
    return cand_gal, cand_coord


def bloom_sigma(rmag):
    """
    Estimated incidence of galaxies per sq arcsec with r > rmag

    Args:
        rmag (float or np.ndarray):

    Returns:
        float or np.ndarray:  Galaxy density

    """
    # Sigma(m)
    sigma = 1. / (3600. ** 2 * 0.334 * np.log(10)) * 10 ** (0.334 * (rmag - 22.963) + 4.320)
    return sigma


def prior_uniform():
    pass


def pchance(rmag, sep, r_half, sigR, scale_rhalf=3., nsigma=3., ndens_eval='bloom'):
    """

    Args:
        rmag (np.ndarray):
        sep (np.ndarray):
        r_half (np.ndarray):
            Half light radius of the galaxy
        sigR (float):
            1 sigma error in FRB localization; assumed symmetric
        scale_rhalf:
        nsigma:

    Returns:

    """

    # Reff - More conservative than usual
    Rs = np.stack([scale_rhalf * r_half, np.ones_like(r_half)* nsigma * sigR,
                   np.sqrt(sep ** 2 + (scale_rhalf * r_half) ** 2)])
    reff = np.max(Rs, axis=0)

    # Number density
    if ndens_eval =='bloom':
        nden = bloom_sigma(rmag)
    elif ndens_eval =='eb17':
        embed(header='83 of pchance')
        # SPEED UP THE FOLLOWING!  SPLINE IT TOO
        ndens = hosts.quad(n_gal, 0, rmag)[0]
    else:
        raise IOError("Bad ndens evaluation")

    # Nbar
    Nbar = np.pi * reff ** 2 * nden

    # Return Pchance
    return 1. - np.exp(-Nbar)


def raw_prior_Mi(Pchance, method):
    """
    Raw prior for a given set of Pchance values

    Proper normalization requires P(S) so that is done below

    Args:
        Pchance (np.ndarray):
        method (str):

    Returns:
        np.ndarray:

    """
    if method == 'linear':
        return 1 - Pchance
    elif method == 'inverse':
        return 1./Pchance
    elif method == 'uniform':
        return np.ones_like(Pchance)
    else:
        raise IOError("Bad method for prior_Mi")


def prior_S_n(rlim, theta_max, rmax=30.):
    rval = np.linspace(rlim, rmax, 100)
    # Nbar
    sigma = bloom_sigma(rval)
    Nbar = np.pi * theta_max ** 2 * sigma
    # Average
    P_S = np.sum(sigma * np.exp(-Nbar)) / np.sum(sigma)
    # Return
    return P_S


def pw_Mi(r_w, r_half, theta_prior):
    """

    Args:
        r_w (np.ndarray):
            offset from galaxy center in arcsec
        r_half (float):
            Half-light radius of this galaxy in arcsec
        theta_prior (dict):
            Parameters for theta prior

    Returns:
        np.ndarray: Probability values; un-normalized

    """
    p = np.zeros_like(r_w)
    ok_w = r_w < theta_prior['max']*r_half
    if theta_prior['method'] == 'rcore':
        if np.any(ok_w):
            p[ok_w] = r_half / (r_w[ok_w] + r_half)
    elif theta_prior['method'] == 'uniform':
        if np.any(ok_w):
            p[ok_w] = 1.
    elif theta_prior['method'] == 'exp':
        if np.any(ok_w):
            p[ok_w] = (r_w[ok_w] / r_half) * np.exp(-r_w[ok_w]/r_half)
    else:
        raise IOError("Bad theta method")
    #
    return p


def px_Mi(box_radius, frb_coord, eellipse, cand_coords,
          theta_prior, nsamp=1000):
    """
    Calculate p(x|M_i)

    Args:
        box_radius (float):
            Maximum radius for analysis, in arcsec
        eellipse (dict):
            Error ellipse for the FRB
            a, b, theta
        frb_coord (SkyCoord):
        cand_coords (SkyCoord):
            Coordinates of the candidate hosts
        theta_prior (dict):
            Parameters for theta prior
        nsamp (int, optional):

    Returns:
        np.ndarray:

    """
    # Error ellipse
    pa_ee = eellipse['theta'] # deg
    dtheta = 90. - pa_ee  # Place a of ellipse along the x-axis
    # Set Equinox (for spherical offsets)
    frb_coord.equinox = cand_coords[0].equinox
    #
    x = np.linspace(-box_radius, box_radius, nsamp)
    xcoord, ycoord = np.meshgrid(x,x)

    # Build the grid around the FRB (orient a on our x axis)
    l_w = np.exp(-xcoord ** 2 / (2 * eellipse['a'] ** 2)) * np.exp(-ycoord ** 2 / (2 * eellipse['b'] ** 2))

    p_xMis = []
    for icand, cand_coord in enumerate(cand_coords):

        # Calculate observed FRB location
        #dra, ddec = cand_coord.spherical_offsets_to(frb_coord)
        #xFRB = -dra.to('arcsec').value
        #yFRB = ddec.to('arcsec').value

        # #####################
        # l(w) -- 2D Gaussian

        # Rotate the galaxy
        r = frb_coord.separation(cand_coord).to('arcsec')
        pa_gal = frb_coord.position_angle(cand_coord).to('deg')
        new_pa_gal = pa_gal + dtheta * units.deg

        #r_wsq = (xcoord-xFRB)**2 + (ycoord-yFRB)**2
        #l_w = np.exp(-r_wsq/(2*sigR**2)) / sigR / np.sqrt(2*np.pi)

        # p(w|M_i)
        # x, y gal
        x_gal = -r.value * np.sin(new_pa_gal).value
        y_gal = r.value * np.cos(new_pa_gal).value
        r_w = np.sqrt((xcoord-x_gal)**2 + (ycoord-y_gal)**2)
        p_wMi = pw_Mi(r_w, theta_prior['r_half'][icand], theta_prior)

        # Product
        grid_p = l_w * p_wMi

        # Average
        p_xMis.append(np.mean(grid_p))
    # Return
    return np.array(p_xMis)


def renorm_priors(raw_Mi, S):
    """
    Simple method to normalize the Priors

    Args:
        raw_Mi (np.ndarray):
        raw_S (float):

    Returns:
        tuple: Normalized priors

    """
    raw_sum = np.sum(raw_Mi)
    return (1.-S) * raw_Mi/raw_sum


def mock_run(sky, frbs, sigR, theta, fov, scale_rhalf=2., nsigma=2.,
             ncontam=0, seed=12345):
    rstate = np.random.RandomState(seed) # For contaminants
    # Coords
    obj_coord = SkyCoord(ra=sky['ra'], dec=sky['dec'], unit='deg')
    frb_coords = SkyCoord(ra=frbs['ra'], dec=frbs['dec'], unit='deg')
    #
    model_dict = {}
    for ifrb in range(len(frb_coords)):
        key = '{:04d}'.format(ifrb)
        if (ifrb % 10) == 0:
            print("Working on FRB: {}".format(key))
        model_dict[key] = {}
        # Grab candidates -- Assume always at least 1
        cands = frb_coords[ifrb].separation(obj_coord) < fov
        cand_gal = sky[cands]
        cand_coord = obj_coord[cands]
        # Contaminants?
        if ncontam != 0:
            cand_gal, cand_coord = add_contam(
                ncontam, frb_coords[ifrb], cand_gal, cand_coord,
                fov.to('arcsec').value, rstate=rstate)
        # Additional
        cand_gal['r_half'] = 1.
        cand_sep = frb_coords[ifrb].separation(cand_coord).to('arcsec')
        # Save
        model_dict[key]['ncand'] = len(cand_gal)
        model_dict[key]['cand_mag'] = cand_gal['DES_r'].data
        model_dict[key]['theta'] = cand_sep.to('arcsec').value
        # Raw P(Mi)
        raw_prior_Mi = prior_Mi_n(cand_gal['DES_r'].data, cand_sep.to('arcsec').value,
                                           cand_gal['r_half'].data, sigR.to('arcsec').value,
                                           scale_rhalf=scale_rhalf, nsigma=nsigma)
        # Raw P(S)
        raw_prior_S = 1. - np.max(raw_prior_Mi)
        # Normalize priors
        prior_Mi, prior_S = renorm_priors(raw_prior_Mi, raw_prior_S)
        # Calculate p(x|Mi)
        p_xMi = px_Mi(fov.to('arcsec').value, frb_coords[ifrb],
                               cand_coord, theta, sigR.to('arcsec').value)
        # Calcualte p(x|S)
        p_xS = px_Mi(fov.to('arcsec').value, frb_coords[ifrb],
                     frb_coords[ifrb:ifrb + 1],  # Needs to be an array
                     theta, sigR.to('arcsec').value)[0]  # Converts array to float
        # Calculate p(x)
        p_x = prior_S * p_xS + np.sum(prior_Mi * p_xMi)
        # Evaluate p(Mi|x)
        P_Mix = prior_Mi * p_xMi / p_x
        model_dict[key]['P_Mix'] = P_Mix
        # Evaluate p(S|x)
        P_Sx = prior_S * p_xS / p_x
        model_dict[key]['P_Sx'] = P_Sx
    # Return
    return model_dict


def parse_model(model_dict, scale=1000, mag_lim=25.5):
    """

    Args:
        model_dict:
        scale:
        mag_lim:

    Returns:

    """
    model_mags = []
    model_theta = []
    max_PMix = []
    #
    for key in model_dict.keys():
        # Observed candidates
        for kk, p in enumerate(model_dict[key]['P_Mix']):
            N = int(np.round(p * scale))
            model_mags += [model_dict[key]['cand_mag'][kk]] * N
            model_theta += [model_dict[key]['theta'][kk]] * N
        # Unseen
        N = int(np.round(model_dict[key]['P_Sx'] * scale))
        model_mags += [mag_lim] * N
        # Misc
        max_PMix.append(np.max(model_dict[key]['P_Mix']))
        # import pdb; pdb.set_trace()
    # Arrays
    model_mags = np.array(model_mags)
    model_theta = np.array(model_theta)
    max_PMix = np.array(max_PMix)
    # Return
    return model_mags, model_theta, max_PMix
