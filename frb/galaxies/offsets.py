""" Methods related to offsets"""

import numpy as np

from astropy import units

from IPython import embed

def angular_offset(frb, galaxy, nsigma=5., nsamp=2000):
    """

    Warning: All calculations in arcsec -- Do not use for *large* localization error

    Args:
        frb (frb.frb.FRB):
        galaxy (frb.galaxies.FRBGalaxy):

    Returns:
        tuple: float, float, float, float
            avg_off, sig_off, best_off, sig_best

    """
    # Error ellipse
    sig_a = frb.eellipse['a']  # arcsec
    sig_b = frb.eellipse['b']  # arcsec

    # If systematic exists, add in quadrature
    if 'a_sys' in frb.eellipse.keys():
        sig_a = np.sqrt(frb.eellipse['a_sys']**2 + sig_a**2)
        sig_b = np.sqrt(frb.eellipse['b_sys']**2 + sig_b**2)

    # FRB is the reference frame
    pa_ee = frb.eellipse['theta'] # deg
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

if __name__ == '__main__':

    from frb.frb import FRB
    ifrb = FRB.by_name('FRB200430')
    host = ifrb.grab_host()
    angular_offset(ifrb, host)
