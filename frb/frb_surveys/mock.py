""" Generate a set of fake FRBs for testing """

import numpy as np

from scipy.interpolate import interp1d
import pandas

from zdm.chime import grids

from frb.frb_surveys import chime
from frb.galaxies import hosts
from frb.defs import frb_cosmo

from IPython import embed

def for_chime(NFRB:int=10000, cut_fluence:float=5.):

    # Load CHIME Dr1
    df_dr1 = chime.load_catalog()
    flu_cut = df_dr1['fluence'] >= cut_fluence
    df_dr1 = df_dr1[flu_cut].copy()

    # Load host galaxy M_r
    xvals, prob1 = hosts.load_Mr_pdf()

    # Load p(z|DM)
    dmvals, zvals, all_rates, all_singles, all_reps =\
        grids.load()
    
    # Cumulative
    cum_all = np.cumsum(all_rates, axis=0)
    norm = np.outer(np.ones(zvals.size), cum_all[-1,:]) 
    cum_all /= norm
    cum_all[0,:] = 0.

    # Interpolators
    print("Building interpolators")
    fs = [interp1d(cum_all[:,ii], zvals) for ii in range(dmvals.size)]
    DM_hosts = np.linspace(1e-3,2000.,2000)
    DMh_pdf = lognorm_pdf(DM_hosts)
    cum_DMh = np.cumsum(DMh_pdf)
    cum_DMh[0] = 0.
    fh = interp1d(cum_DMh/cum_DMh[-1], DM_hosts)


    # Random numbers
    rand = np.random.uniform(size=NFRB)
    rand_DMex = np.random.choice(
        df_dr1['DMex'], size=NFRB)

    # DM_cosmic 
    DM_cosmic = np.zeros(NFRB)
    DM_host = np.zeros(NFRB)
    while np.any(DM_cosmic == 0.):
        rand_DMh = np.random.uniform(size=NFRB)
        DMh = fh(rand_DMh)
        tmp = rand_DMex - DMh
        # Keep
        gd = (tmp > 0.) & (DM_cosmic == 0.)
        DM_cosmic[gd] = tmp[gd]
        DM_host[gd] = DMh[gd]

    # Ugly for loops
    zs = []
    for kk,DMc in enumerate(DM_cosmic):
        imin = np.argmin(np.abs(dmvals-DMc))
        z = fs[imin](rand[kk])
        zs.append(float(z))
    zs = np.array(zs)

    # Now apparent magnitude
    M_r = np.random.choice(xvals, NFRB, p=prob1)
    #Index_ = np.random.shuffle(np.arange(len(Ar)))
    #mw_extinction = np.squeeze(Ar[Index_])
    dist_mod = frb_cosmo.distmod(zs).value
    host_m_r = dist_mod + M_r #+ mw_extinction

    # Build FRB table
    df = pandas.DataFrame()
    df['DMex'] = rand_DMex
    df['DM_cosmic'] = DM_cosmic
    df['DM_host'] = DM_host
    df['z'] = zs
    df['M_r'] = M_r
    df['m_r'] = host_m_r

    #embed(header='87')

    # Return
    return df


def lognorm_pdf(DM, logs:tuple=None):
    if logs is None:
        logmean = 2.16 / 0.4342944619
        logsigma = 0.51 / 0.4342944619
    else:
        logmean, logsigma = logs
    # 
    logDM = np.log(DM)
    
    norm = (2.0 * np.pi) ** -0.5 / DM / logsigma
    return norm * np.exp(-0.5 * ((logDM - logmean) / logsigma) ** 2)


if __name__ == '__main__':
    for_chime()