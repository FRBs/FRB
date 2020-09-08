"""Methods related to Bayesian association analysis"""

import numpy as np


def bloom_sigma(rmag):
    # Sigma(m)
    sigma = 1. / (3600. ** 2 * 0.334 * np.log(10)) * 10 ** (0.334 * (rmag - 22.963) + 4.320)
    return sigma


def prior_uniform():
    pass

def prior_Mi_n(rmag, sep, r_half, sigR, scale_rhalf=3., nsigma=3.):
    # Reff - More conservative than usual
    Rs = np.stack([3 * r_half, np.ones_like(r_half)* nsigma * sigR,
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
