""" Codes realted to developing what has become PATH (aka astropath)"""

import pdb; pdb.set_trace()  # This module is defunct

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
