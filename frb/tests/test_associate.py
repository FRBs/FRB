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


    '''  This works on the 180301 branch
    # We skirt the usual candidate init
    frbA.candidates['mag'] = frbA.candidates[frbA.filter]
    frbA.init_cand_coords()
    # Set priors
    frbA.init_cand_prior('inverse', P_U=0.)
    frbA.init_theta_prior('inverse', 6.)

    # Localization
    frbA.init_localization('eellipse', 
                            center_coord=frbA.frb.coord,
                            eellipse=frbA.frb_eellipse)
    
    # Calculate priors
    frbA.calc_priors()                            

    # Calculate p(O_i|x)
    frbA.calc_posteriors('fixed', box_hwidth=frbA.max_radius)
    '''