import sys, os
import numpy as np
from scipy.stats import norm

from astropy.table import Table

sys.path.append(os.path.abspath("../Analysis/py"))
import alopeke_defs

def calc_zeropoint(data_dict, camera:str):

    # Instrumental magnitude for star 1 and star 2
    m_star1 = -2.5 * np.log10(data_dict['C_star1_full']-data_dict['C_bkg_full'])
    m_star2 = -2.5 * np.log10(data_dict['C_star2_full']-data_dict['C_bkg_full'])

    # Uncertainties on instrumental magnitudes - assume Poisson
    merr_star1 = 1.086 * np.sqrt(data_dict['C_star1_full'])/data_dict['C_star1_full']
    merr_star2 = 1.086 * np.sqrt(data_dict['C_star1_full'])/data_dict['C_star1_full']

    # Now get zero points
    zpt_star1 = 0.
    zpt_star2 = 0.
    if camera=='red':
        zpt_star1 = alopeke_defs.i_1 - m_star1
        zpt_star2 = alopeke_defs.i_2 - m_star2
    elif camera=='blue':
        zpt_star1 = alopeke_defs.r_1 - m_star1
        zpt_star2 = alopeke_defs.r_2 - m_star2
    else:
        raise Exception(f'ERROR: Unrecognized camera {camera}')

    zpt = (merr_star1 * zpt_star1 + merr_star2 * zpt_star2)/(merr_star1 + merr_star2)
    zpt_err = 1./np.sqrt(2) * np.sqrt(merr_star1**2 + merr_star2**2)

    return(zpt, zpt_err)

def get_star_flux(date:str, data_dict):
    mask = np.abs(data_dict['MJD_star1_full']-alopeke_defs.FRB_time[date].mjd)*24*3600<0.163

    data_dict['mean_bkg'] = np.mean(data_dict['C_bkg_full'][mask])
    data_dict['mean_star1'] = np.mean(data_dict['C_star1_full'][mask])
    data_dict['mean_star2'] = np.mean(data_dict['C_star2_full'][mask])
    data_dict['mean_star3'] = np.mean(data_dict['C_star3_full'][mask])

    return(data_dict)

def load_camera(camera:str, date:str, cut=True):

    if camera == 'red':
        data = Table.read(f'../../Data/master_table_{date}_r.fits')
    elif camera == 'blue':
        data = Table.read(f'../../Data/master_table_{date}_b.fits')
    else:
        raise IOError("Bad camera")

    # Cut down to Galaxy
    i_FRB = (data['MJD']>=alopeke_defs.mjd_low[date]) & (
        data['MJD']<=alopeke_defs.mjd_high[date])
    i_gal = np.invert(i_FRB)

    # Expunge first reads
    if cut:
        gd_data = data['frame'] > 1
    else:
        gd_data = np.ones(len(data), dtype='bool')

    # ADUs
    C_gal = data['flux_2FWHM'][i_gal & gd_data]
    C_gal_full = data['flux_2FWHM'][gd_data]
    C_frb = data['flux_2FWHM'][i_FRB & gd_data]

    C_star1 = data['flux_star1_2FWHM'][i_gal & gd_data]
    C_star1_full = data['flux_star1_2FWHM'][gd_data]

    C_star2 = data['flux_star2_2FWHM'][i_gal & gd_data]
    C_star2_full = data['flux_star2_2FWHM'][gd_data]

    C_star3 = data['flux_star3_2FWHM'][i_gal & gd_data]
    C_star3_full = data['flux_star3_2FWHM'][gd_data]

    C_bkg = data['flux_bkg_2FWHM'][i_gal & gd_data]
    C_bkg_full = data['flux_bkg_2FWHM'][gd_data]
    

    out_dict = {}

    # Time
    out_dict['MJD_gal'] = data['MJD'][i_gal & gd_data]
    out_dict['MJD_FRB'] = data['MJD'][i_FRB & gd_data]
    out_dict['MJD_star1'] = data['MJD'][i_gal & gd_data]
    out_dict['MJD_star1_full'] = data['MJD'][gd_data]
    out_dict['MJD_star2'] = data['MJD'][i_gal & gd_data]
    out_dict['MJD_star2_full'] = data['MJD'][gd_data]
    out_dict['MJD_star3'] = data['MJD'][i_gal & gd_data]
    out_dict['MJD_star3_full'] = data['MJD'][gd_data]
    out_dict['dt'] = (data['MJD'][2]-data['MJD'][1])*24*3600 # seconds

    # Counts
    gain = alopeke_defs.gain_red if camera == 'red' else alopeke_defs.gain_blue
    out_dict['C_FRB'] = C_frb * gain / alopeke_defs.EM
    out_dict['C_gal'] = C_gal* gain / alopeke_defs.EM
    out_dict['C_gal_full'] = C_gal_full * gain / alopeke_defs.EM

    out_dict['C_star1'] = C_star1 * gain / alopeke_defs.EM
    out_dict['C_star1_full'] = C_star1_full * gain / alopeke_defs.EM

    out_dict['C_star2'] = C_star2 * gain / alopeke_defs.EM
    out_dict['C_star2_full'] = C_star2_full * gain / alopeke_defs.EM

    out_dict['C_star3'] = C_star3 * gain / alopeke_defs.EM
    out_dict['C_star3_full'] = C_star3_full * gain / alopeke_defs.EM

    out_dict['C_bkg'] = C_bkg * gain / alopeke_defs.EM
    out_dict['C_bkg_full'] = C_bkg_full * gain / alopeke_defs.EM
    

    # Gauss
    out_dict['Gauss_gal'] = norm.fit(out_dict['C_gal'])
    out_dict['Gauss_star1'] = norm.fit(out_dict['C_star1'])
    out_dict['Gauss_star2'] = norm.fit(out_dict['C_star2'])
    out_dict['Gauss_star3'] = norm.fit(out_dict['C_star3'])

    out_dict = get_star_flux(date, out_dict)

    zpt, zpt_err = calc_zeropoint(out_dict, camera)
    out_dict['zeropoint'] = zpt
    out_dict['zeropoint_error'] = zpt_err

    # Return
    return out_dict
