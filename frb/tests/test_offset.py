from pkg_resources import resource_filename

import numpy as np


from frb.frb import FRB, build_table_of_frbs
from frb.galaxies import offsets as frb_offsets



def test_offset():
    ifrb = FRB.by_name('FRB20121102A')
    host = ifrb.grab_host()

    ra_sig_source = host.positional_error['ra_source']
    ra_sig_astro = host.positional_error['ra_astrometric']
    dec_sig_source = host.positional_error['dec_source']
    dec_sig_astro = host.positional_error['dec_astrometric']

    host_ra_sig = np.sqrt(ra_sig_astro ** 2 + ra_sig_source ** 2)
    host_dec_sig = np.sqrt(dec_sig_astro ** 2 + dec_sig_source ** 2)

    ang_avg, avg_err, ang_best, best_err = frb_offsets.angular_offset(
        ifrb, host, gal_sig=(host_ra_sig, host_dec_sig))

    assert np.isclose(ang_best, 0.22687, rtol=1e-5)