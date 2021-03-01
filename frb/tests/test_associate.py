import os
import numpy as np
from pkg_resources import resource_filename

import pandas

from astropy import units

from frb.associate import frbassociate

import pytest

remote_data = pytest.mark.skipif(os.getenv('FRB_GDB') is None,
                                 reason='test requires dev suite')

@remote_data
def test_individual():
    from frb.associate import frbs
    config = getattr(frbs, 'FRB180924'.lower())
    frbA = frbassociate.run_individual(config)
    # Test
    assert isinstance(frbA.candidates, pandas.DataFrame)


@remote_data
def test_step_by_step():
    from frb import frb

    # Testing

    # #########################
    # Prep
    frb180924 = frb.FRB.by_name('FRB180924')

    # Instantiate
    max_radius = 10.  # in arcseconds
    frbA_180924 = frbassociate.FRBAssociate(frb180924,
                                            image_file=resource_filename('frb', 'associate/dev/FRB180924_DESr.fits'),
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
                   ang_size=frbA_180924.candidates['half_light'].values,
                   max=theta_max)
    frbA_180924.set_theta_prior(theta_u)

    # Calcuate p(O_i|x)
    frbA_180924.calc_POx()

    final_cands = frbA_180924.P_Oix > 0.01
    print(frbA_180924.candidates[['id', 'r', 'half_light',
                                  'separation', 'P_O', 'P_Ox']][final_cands])

    # Now Linear prior
    # Theta
    theta_max = 10.  # in half-light units
    theta_c = dict(method='core', max=theta_max)
    theta_c['ang_size'] = frbA_180924.candidates['half_light'].values
    frbA_180924.set_theta_prior(theta_c)

    # Calcuate p(O_i|x)
    frbA_180924.calc_POx()

    final_cands = frbA_180924.P_Oix > 0.01
    print(frbA_180924.candidates[['id', 'r', 'half_light',
                                  'separation', 'P_O', 'P_Ox']][final_cands])
    frbA_180924.candidates['P_Ox_core'] = frbA_180924.P_Oix.copy()


    # Now Exponential prior
    # Theta
    theta_max = 10.  # in half-light units
    theta_e = dict(method='exp', max=theta_max)
    theta_e['ang_size'] = frbA_180924.candidates['half_light'].values
    frbA_180924.set_theta_prior(theta_e)

    # Calcuate p(O_i|x)
    frbA_180924.calc_POx()

    final_cands = frbA_180924.P_Oix > 0.01
    print(frbA_180924.candidates[['id', 'r', 'half_light',
                                  'separation', 'P_O', 'P_Ox']][final_cands])

    assert np.isclose(np.max(frbA_180924.candidates.P_Ox), 0.9990933441057845)



