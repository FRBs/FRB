
from astropy import units

from frb.associate import frbassociate

import pytest

def test_step_by_step():
    from frb import frb

    # Testing

    # #########################
    # Prep
    frb180924 = frb.FRB.by_name('FRB180924')

    # Instantiate
    max_radius = 10.  # in arcseconds
    frbA_180924 = frbassociate.FRBAssociate(frb180924, image_file='dev/FRB180924_DESr.fits',
                               max_radius=max_radius)

    # Threshold
    frbA_180924.threshold()

    # Segment
    frbA_180924.segment()#show=True, outfile='dev/FRB180924_segm.png')

    # Photometry
    frbA_180924.photometry('MAGZERO', 'r')#, show=True, outfile='dev/FRB180924_aper.png')

    # Candidates
    plate_scale = frbA_180924.header['CD2_2'] * 3600. * units.arcsec  # arcsec
    frbA_180924.cut_candidates(plate_scale, bright_cut=18.1)


    # #########################
    # Bayesian time

    # Pchance
    frbA_180924.calc_pchance()

    # Priors
    frbA_180924.calc_priors(0.01)

    # Theta
    theta_max = 10.  # in half-light units
    theta_u = dict(method='uniform',
                   r_half=frbA_180924.candidates['half_light'].data,
                   max=theta_max)
    frbA_180924.set_theta_prior(theta_u)

    # Calcuate p(M_i|x)
    frbA_180924.calc_PMx()

    final_cands = frbA_180924.P_Mix > 0.01
    print(frbA_180924.candidates[['id', 'r', 'half_light',
                                  'separation', 'P_M', 'P_Mx']][final_cands])
    frbA_180924.candidates['P_Mx_u'] = frbA_180924.P_Mix.copy()

    # Now Linear prior
    # Theta
    theta_max = 10.  # in half-light units
    theta_c = dict(method='rcore', max=theta_max)
    theta_c['r_half'] = frbA_180924.candidates['half_light'].value
    frbA_180924.set_theta_prior(theta_c)

    # Calcuate p(M_i|x)
    frbA_180924.calc_PMx()

    final_cands = frbA_180924.P_Mix > 0.01
    print(frbA_180924.candidates[['id', 'r', 'half_light',
                                  'separation', 'P_M', 'P_Mx']][final_cands])
    frbA_180924.candidates['P_Mx_c'] = frbA_180924.P_Mix.copy()


