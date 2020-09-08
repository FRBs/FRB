"""Methods related to Bayesian association analysis"""

import numpy as np

from astropy.table import Table
from astropy import units
from astropy.coordinates import SkyCoord

def bloom_sigma(rmag):
    # Sigma(m)
    sigma = 1. / (3600. ** 2 * 0.334 * np.log(10)) * 10 ** (0.334 * (rmag - 22.963) + 4.320)
    return sigma


def prior_uniform():
    pass

def prior_Mi_n(rmag, sep, r_half, sigR, scale_rhalf=3., nsigma=3.):
    # Reff - More conservative than usual
    Rs = np.stack([scale_rhalf * r_half, np.ones_like(r_half)* nsigma * sigR,
                   np.sqrt(sep ** 2 + (scale_rhalf * r_half) ** 2)])
    reff = np.max(Rs, axis=0)
    # Nbar
    Nbar = np.pi * reff ** 2 * bloom_sigma(rmag)
    # Return
    return np.exp(-Nbar)


def prior_S_n(rlim, theta_max, rmax=30.):
    rval = np.linspace(rlim, rmax, 100)
    # Nbar
    sigma = bloom_sigma(rval)
    Nbar = np.pi * theta_max ** 2 * sigma
    # Average
    P_S = np.sum(sigma * np.exp(-Nbar)) / np.sum(sigma)
    # Return
    return P_S


def px_Mi(box_radius, frb_coord, cand_coords, theta, sigR, nsamp=1000):
    x = np.linspace(-box_radius, box_radius, nsamp)
    xcoord, ycoord = np.meshgrid(x,x)
    r_w = np.sqrt(xcoord**2 + ycoord**2)
    p_xMis = []
    for cand_coord in cand_coords:
        # Center on the galaxy
        # Calculate observed FRB location
        dra, ddec = cand_coord.spherical_offsets_to(frb_coord)
        xFRB = -dra.to('arcsec').value
        yFRB = ddec.to('arcsec').value
        # l(w) -- Gaussian
        r_wsq = (xcoord-xFRB)**2 + (ycoord-yFRB)**2
        l_w = np.exp(-r_wsq/(2*sigR**2)) / sigR / np.sqrt(2*np.pi)
        # p(w|M_i)
        p_wMi = np.zeros_like(xcoord)
        if theta['method'] == 'rcore':
            ok_w = r_w < theta['max']
            p_wMi[ok_w] = 1./(r_w + theta['core'])
        elif theta['method'] == 'uniform':
            ok_w = r_w < theta['max']
            p_wMi[ok_w] = 1.
        else:
            raise IOError("Bad theta method")
        # Product
        grid_p = l_w * p_wMi
        # Average
        p_xMis.append(np.mean(grid_p))
    # Return
    return np.array(p_xMis)


def renorm_priors(raw_Mi, raw_S):
    raw_sum = np.sum(raw_Mi) + raw_S
    return raw_Mi/raw_sum, raw_S/raw_sum


def mock_run(sky, frbs, sigR, theta, fov, scale_rhalf=2., nsigma=2.):
    # Coords
    obj_coord = SkyCoord(ra=sky['ra'], dec=sky['dec'], unit='deg')
    frb_coord = SkyCoord(ra=frbs['ra'], dec=frbs['dec'], unit='deg')
    #
    model_dict = {}
    for ifrb in range(len(frb_coord)):
        key = '{:04d}'.format(ifrb)
        if (ifrb % 10) == 0:
            print("Working on FRB: {}".format(key))
        model_dict[key] = {}
        # Grab candidates -- Assume always at least 1
        cands = frb_coord[ifrb].separation(obj_coord) < fov
        cand_gal = sky[cands]
        cand_gal['r_half'] = 1.
        cand_coord = obj_coord[cands]
        cand_sep = frb_coord[ifrb].separation(cand_coord).to('arcsec')
        # Save
        model_dict[key]['ncand'] = len(cands)
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
        p_xMi = px_Mi(fov.to('arcsec').value, frb_coord[ifrb],
                               cand_coord, theta, sigR.to('arcsec').value)
        # Calcualte p(x|S)
        p_xS = px_Mi(fov.to('arcsec').value, frb_coord[ifrb], frb_coord[ifrb:ifrb + 1],
                              theta, sigR.to('arcsec').value)[0]
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
