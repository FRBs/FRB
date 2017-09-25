""" Module for assessing impact of intervening galaxies
   (DLAs) on FRB measurements
   Based on calclations presented in Prochaska & Neeleman 2017
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

def avgDM(z, use_boot=False):
    """ Calculate the average DM from intervening galaxies

    Parameters
    ----------
    z : float or ndarray
    use_boot

    Returns
    -------
    avgDM : float or ndarray (depending on type of input z)

    """
    # Init
    if isinstance(z, float):
        flg_float = True
        z = np.array([z])
    else:
        flg_float = False
    #
    if use_boot:
        pdb.set_trace()
    else:
        # Load model
        atan_lz = fda.model_lz('arctan', path=analy_pth)

        dz = np.median(z - np.roll(z, 1))
        lz = atan_lz['eval'](z, param=atan_lz['param'])

        # Integral h(N) dN
        # zero_h = fda.int_dbl_pow()

        # Average NHI -- This includes the zero_h term, so we are fine
        avgNHI = fda.avgN_dbl_pow()

        # Take ne/nH
        nenH_p = fda.nenH_param
        nenH = nenH_p['bp'] + nenH_p['m'] * (avgNHI - 20.3)

        # Integrate lz for n(z)
        cumul = np.cumsum(lz * dz)

        # Average <z>
        avgz = np.cumsum(z * lz * dz) / cumul

        # <DM> for a single DLA (rest-frame)
        DM_DLA = 10. ** (avgNHI + nenH) / u.cm ** 2
        print("DM for an average DLA = {} (rest-frame)".format(DM_DLA.to('pc/cm**3')))

        # Altogether now
        avgDM = 10 ** avgNHI * 10 ** nenH * cumul / (1 + avgz) / u.cm ** 2
        # z=1
        iminz = np.argmin(np.abs(z - 1))
        print("<DM> at z=1 is {}".format(avgDM[iminz].to('pc/cm**3')))
