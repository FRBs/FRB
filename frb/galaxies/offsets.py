""" Methods related to offsets"""

import numpy as np

import pandas

from astropy import units
from astropy.coordinates import SkyCoord

from IPython import embed

def angular_offset(frb, galaxy, nsigma=5., nsamp=2000,
                   gal_sig=None):
    """

    Warning: All calculations in arcsec -- Do not use for *large* localization error

    Args:
        frb (frb.frb.FRB):
        galaxy (frb.galaxies.FRBGalaxy):
        nisgma: Number of sigma around the FRB uncertainty ellipse from which the grid is build
        nsamp: Grid sample points for each 1D Gaussian (one in RA and Dec)
        gal_sig (tuple): RA, DEC errors in arcsec as floats

    Returns:
        tuple: float, float, float, float
            avg_off, sig_off, best_off, sig_best

    """
    # Error ellipse
    sig_a = frb.eellipse['a']  # arcsec
    sig_b = frb.eellipse['b']  # arcsec
    pa_ee = frb.eellipse['theta'] # deg

    # If systematic exists, add in quadrature
    if 'a_sys' in frb.eellipse.keys():
        sig_a = np.sqrt(frb.eellipse['a_sys']**2 + sig_a**2)
        sig_b = np.sqrt(frb.eellipse['b_sys']**2 + sig_b**2)

    # Add in Galaxy error : THIS TREATMENT IS APPROXIMATE
    if gal_sig is not None:
        # Project
        sig2_gal_a = gal_sig[1]**2 * np.cos(pa_ee)**2 + gal_sig[0]**2 * np.sin(pa_ee)**2
        sig2_gal_b = gal_sig[0]**2 * np.cos(pa_ee)**2 + gal_sig[1]**2 * np.sin(pa_ee)**2
        # Add em in
        sig_a = np.sqrt(sig_a**2 + sig2_gal_a)
        sig_b = np.sqrt(sig_b**2 + sig2_gal_b)

    # FRB is the reference frame
    dtheta = 90. - pa_ee  # Place a of ellipse along the x-axis

    # Build the grid around the FRB (orient a on our x axis)
    x = np.linspace(-nsigma*sig_a, nsigma*sig_a, nsamp)
    y = np.linspace(-nsigma*sig_b, nsigma*sig_b, nsamp)
    xx, yy = np.meshgrid(x, y)

    # Rotate the galaxy
    r = frb.coord.separation(galaxy.coord).to('arcsec')
    pa_gal = frb.coord.position_angle(galaxy.coord).to('deg')
    new_pa_gal = pa_gal + dtheta*units.deg

    # x, y gal
    x_gal = -r.value * np.sin(new_pa_gal).value
    y_gal = r.value * np.cos(new_pa_gal).value

    # Offset
    ang_off = np.sqrt((xx-x_gal)**2 + (yy-y_gal)**2)
    p_xy = np.exp(-xx**2 / (2*sig_a**2)) * np.exp(-yy**2 / (2*sig_b**2))

    # Best offset
    best_off = r.value
    var_best = np.sum((ang_off-best_off)**2 * p_xy) / np.sum(p_xy)
    sig_best = np.sqrt(var_best)

    # Average over the grid
    avg_off = np.sum(ang_off * p_xy) / np.sum(p_xy)
    var_off = np.sum((ang_off-avg_off)**2 * p_xy) / np.sum(p_xy)
    sig_off = np.sqrt(var_off)

    # Return
    return avg_off, sig_off, best_off, sig_best


def incorporate_hst(hst_astrom:pandas.DataFrame, host):
    """Updates coordinates and offsets for galaxies
    observed with HST

    Currently for Mannings+2021 only

    Args:
        hst_astrom (pandas.DataFrame): [description]
        host (frb.galaxies.frbgalaxy.FRBHost): [description]
    """

    frb_index = host.frb.frb_name[3:]
    hst_idx = np.where(hst_astrom.FRB.values.astype('str') == frb_index)[0]  # get from excel csv

    hst_row = hst_astrom.iloc[hst_idx[0]]

    # Set coordinate
    host.coord = SkyCoord(
        hst_row.RA + ' ' + hst_row.Dec, 
        unit=(units.hourangle, units.deg), frame='icrs')

    # Positional errors
    host.positional_error['ra_astrometric'] = hst_row['RA_astrom_sig']
    host.positional_error['dec_astrometric'] = hst_row['DEC_astrom_sig']
    host.positional_error['ra_source'] = hst_row['RA_source_sig']
    host.positional_error['dec_source'] = hst_row['DEC_source_sig']

    # combine errors for offset estimation
    host_ra_sig = np.sqrt(host.positional_error['ra_astrometric'] ** 2 + host.positional_error['ra_source'] ** 2)
    host_dec_sig = np.sqrt( host.positional_error['dec_astrometric'] ** 2 + host.positional_error['dec_source'] ** 2)

    # Angular offset
    host.offsets['ang_avg'], host.offsets['ang_avg_err'], \
    host.offsets['ang_best'], host.offsets['ang_best_err'] \
        = angular_offset(host.frb, host, gal_sig=(host_ra_sig, host_dec_sig))


if __name__ == '__main__':

    from frb.frb import FRB
    ifrb = FRB.by_name('FRB200430')
    host = ifrb.grab_host()
    angular_offset(ifrb, host)
