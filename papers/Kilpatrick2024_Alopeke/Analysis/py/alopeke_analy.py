""" Analysis methods """
import sys, os
import numpy as np
import requests

from scipy.stats import poisson
from scipy.interpolate import interp1d

from astropy import units
from astropy.time import Time

from frb.galaxies import photom, nebular

sys.path.append(os.path.abspath("../../Analysis/py"))
import alopeke_defs
import alopeke_utils2

from IPython import embed

import warnings
warnings.filterwarnings('ignore')

def dust_extinction(camera:str):
    """Calculate the dust extinction 

    Args:
        camera (str): [description]

    Raises:
        IOError: [description]

    Returns:
        [type]: [description]
    """
    try:
        EBV = float(nebular.get_ebv(alopeke_defs.frb180916.coord,
                              definition='SandF')['meanValue'])
    except requests.exceptions.ConnectionError:
        raise IOError("No internet connection!")

    if camera == 'red':
        return EBV, photom.extinction_correction('GMOS_N_i', EBV)
    elif camera == 'blue':
        return EBV, photom.extinction_correction('GMOS_N_r', EBV)
    else:
        raise IOError("Bad camera!")

def time_obs(date:str, camera:str):
    # Load
    data_dict = alopeke_utils2.load_camera(camera, date)

    # t0
    FRB_MJD = alopeke_defs.chime_time_arrival[date].mjd
    t_start = (FRB_MJD - data_dict['MJD_gal'][0])*24*3600 / 60. # min
    t_end = (data_dict['MJD_gal'][-1] - FRB_MJD)*24*3600 / 60. # min

    time_of_start = Time(data_dict['MJD_gal'][0], format='mjd')
    sec_frac =  float(time_of_start.datetime.strftime(".%f"))
    sec_frac = '%.3f'%sec_frac
    sec_frac = sec_frac[1:]
    time_of_start_str = time_of_start.datetime.strftime("%Y-%m-%dT%H:%M:%S")+sec_frac
    
    # Report
    print(f"This camera = {camera}")
    print(f"Time of start (UTC) {time_of_start_str}")
    print(f"We started {t_start} minutes before")
    print(f"And ended {t_end} minutes after")

def time_evolution(date:str, camera:str, t_rel=59145):
    # Load
    data_dict = alopeke_utils2.load_camera(camera, date)

    time_array = (data_dict['MJD_gal'].data-alopeke_defs.FRB_time[date].mjd)*24*3600
    mask = np.abs(time_array) < 1.0

    # Fit
    p, V = np.polyfit(
        time_array[mask], data_dict['C_star1'][mask], 1,
        cov=True)
    
    # Report
    print("Intercept is: {0} electrons/s".format(p[1]))
    print("Slope is: {} electrons/s".format(p[0]))
    return p, V

def prob_of_max(date:str, camera:str, nsamp=10000):
    data_dict = alopeke_utils2.load_camera(camera, date, cut=True)

    n_event = len(data_dict['C_FRB'])
    max_CFRB = data_dict['C_FRB'].data.max()

    # Random draws from C_gal
    n_gal = len(data_dict['C_gal'])
    rand_num = np.random.choice(data_dict['C_gal'], size=(n_event, nsamp), replace=True)

    exceed = np.sum(rand_num > max_CFRB, axis=0)
    f_exceed = np.sum(exceed > 0)/exceed.size

    # Return
    print("Fraction exceeding: {}".format(f_exceed))
    return f_exceed

def values_around_frb(date:str, camera:str):

    # Load
    data_dict = alopeke_utils2.load_camera(camera, date)

    n_values = len(data_dict['C_FRB'])
    max_flux = np.max(data_dict['C_FRB'])
    std_flux = np.std(data_dict['C_FRB'])
    mean_flux = np.mean(data_dict['C_FRB'])

    ata = alopeke_defs.ata

    print(f'N values (+/-{ata}): {n_values}')
    print(f'Max flux: {max_flux}')
    print(f'Mean flux: {mean_flux}')
    print(f'Standard deviation: {std_flux}')

def sensitivity_function(date:str, camera:str):

    # Load
    data_dict = alopeke_utils2.load_camera(camera, date)

    zpt = np.nanmedian(data_dict['zeropoint'])
    zpt_err = np.nanmedian(data_dict['zeropoint_error'])

    print(f'Zeropoint: {zpt} AB mag')
    print(f'Zeropoint error: {zpt_err} mag')
    
def total_time(camera, date):
    # Load
    data_dict = alopeke_utils2.load_camera(camera, date)
    Dt = data_dict['MJD_gal'][-1] - data_dict['MJD_gal'][0]
    Dt *= 24*3600
    print("Total observing = {}s".format(Dt))

def upper_limit(date:str, camera:str, nsamp:int=100, step:int=2, cl:float=2.7e-3):
    """ Calculate the upper limit to the counts

    Args:
        camera (str): [description]
        nsamp (int, optional): [description]. Defaults to 100.
        step (int, optional): [description]. Defaults to 2.
        cl (float, optional): [description]. Defaults to 1e-3.

    Returns:
        [type]: [description]
    """
                
    if camera == 'red':
        assert alopeke_defs.EM == 1000
        # Range of analysis
        min_FRB=20
        max_FRB=160
    elif camera == 'blue':
        assert alopeke_defs.EM == 1000
        # Range of analysis
        min_FRB=20
        max_FRB=160
    else:
        raise IOError(f"Bad camera: {camera}")

    data_dict = alopeke_utils2.load_camera(camera, date)
    C_gal = data_dict['C_gal']
    C_grid = np.outer(C_gal, np.ones(nsamp))

    # Max observed
    max_C_obs_FRB = np.max(data_dict['C_FRB'])

    C_FRB_val = np.arange(min_FRB, max_FRB+step, step=step)
    p_exceed = []

    for C_FRB in C_FRB_val:
        rand_C_FRB = poisson.rvs(C_FRB, size=(C_gal.size, nsamp))
        rand_C_obs = rand_C_FRB + C_grid
        # Stats
        p_exceed.append(np.sum(rand_C_obs > max_C_obs_FRB)/rand_C_obs.size)
    p_exceed = np.array(p_exceed)

    # Interpolate
    f_int = interp1d(p_exceed, C_FRB_val)
    # import pdb;pdb.set_trace()
    upper_FRB = float(f_int(1-cl))

    print(f"The upper limit is {upper_FRB}")

    # Return
    return C_FRB_val, p_exceed, upper_FRB

def calc_fluence(camera:str, date:str):
    # Grab the values
    data_dict = alopeke_utils2.load_camera(camera, date)
    C_FRB_val, p_exceed, upper_FRB = upper_limit(date, camera)
    EBV, A = dust_extinction(camera)

    # Magnitude
    Amag = 2.5*np.log10(A)

    # Star-1
    mean_1_tot = data_dict['mean_star1']
    bkg_median = np.median(data_dict['mean_bkg']) # Background
    mean_1 = mean_1_tot - bkg_median #(alopeke_defs.C_red_bg if camera == 'red' else alopeke_defs.C_blue_bg)

    # Reference Zero point
    ZP = 2.5*np.log10(mean_1) + (alopeke_defs.i_1 if camera == 'red' else alopeke_defs.r_1)

    # Convert FRB limit to apparent magnitude and apply extinction
    mag_FRB = -2.5*np.log10(upper_FRB) + ZP - Amag
    mag_FRB_uncorr = -2.5*np.log10(upper_FRB) + ZP

    # AB -- Include a color term??
    fnu = 10**(-(48.6+mag_FRB)/2.5) * units.erg/units.s/units.Hz/units.cm**2
    fnu_uncorr = 10**(-(48.6+mag_FRB_uncorr)/2.5) * units.erg/units.s/units.Hz/units.cm**2

    # Fluence
    fluence = (fnu * alopeke_defs.dt_alopeke).to('uJy s')
    fluence_uncorr = (fnu_uncorr * alopeke_defs.dt_alopeke).to('uJy s')

    freq = None
    if camera=='blue':
        freq = alopeke_defs.r_nu
    elif camera=='red':
        freq = alopeke_defs.i_nu
    else:
        raise IOError("Bad camera!")

    energy = fnu * alopeke_defs.dt_alopeke * freq * 4 * np.pi * alopeke_defs.distance**2
    energy = energy.to('erg')

    # Get the radio fluence
    radio_fluence = alopeke_defs.radio_data[date]["fit_statistics"]["bestfit_parameters"]["fluence"][0]
    # Radio fluence is in Jy ms
    radio_fluence = radio_fluence * 1.0e-23 * 1.0e-3 * units.erg/units.Hz/units.cm**2

    fluence_ratio = fnu * alopeke_defs.dt_alopeke / radio_fluence

    print(f"Extinction (A_filter): {Amag}")
    print(f"AB mag (no dust correction): {mag_FRB_uncorr}")
    print(f"AB mag (MW dust correction): {mag_FRB}")
    print(f"fluence (no dust correction): {fluence_uncorr}")
    print(f"fluence (MW dust correction): {fluence}")
    print(f"fluence ratio (opt/radio): {fluence_ratio}")
    print(f"Equivalent isotropic energy: {energy}")
    return fluence, mag_FRB

# See equation 63 from Metzger paper, defined piecewise in time-evolving
# cooling and synchrotron frequency.
# Assume t in seconds, nu in Hz
def luminosity(nu, t, E43=100.0, sigma=0.3, beta=0.5, M21=1.0, T=1.4e6, dt=1e-4, nH=None):
    # nu_syn is defined piecewise.  See equations 56-57 in Metzger
    nu_syn_0 = 1.38e22 * (10 * sigma)**0.5 * E43**0.5 * (1e3 * dt)**-1.5
    # Note that t_dec in eqn 56-57 is approximately dt (eqn 13)
    if t < dt:
        nu_syn = nu_syn_0 * (t/dt)**-1.0
    else:
        nu_syn = nu_syn_0 * (t/dt)**-1.5

    # Can use circumburst density instead of T
    if nH:
        T = 1e5 * np.sqrt((4e3/nH) * M21 * (beta/0.5)**-3)

    nu_c = 2.17e18 * (10 * sigma)**-1.5 * (2 * beta)**3 * M21**-1 *\
        (1.0e3 * t)**-0.5 * (1.0e-5 * T)**2 # in Hz

    L_pk = 1e45 * E43 * (1.0e3 * t)**-1 # in erg/s
    tc = 6.4 * (10 * sigma)**2 * E43**0.5 * (2 * beta)**-3 * M21 * (1.0e-5*T)**-2

    if nu > nu_syn:
        L_pk = L_pk * np.exp(-(nu/nu_syn-1))

    if nu < nu_c:
        return(L_pk*(nu/nu_c)**(1.333333333)*(nu_c/nu_syn)**0.5)
    else:
        return(L_pk*(nu/nu_syn)**0.5)

def afterglow_limits(date:str, Erange=[-1.8, 4.414], nrange=[1,6.2]):

    Evals = np.linspace(*Erange, 400)
    nvals = np.linspace(*nrange, 300)

    grid = np.zeros((len(nvals),len(Evals)))

    print(f'Getting fluence limits for {date}')
    #rfluence, rmag = calc_fluence('red', date)
    #bfluence, bmag = calc_fluence('blue', date)

    rmag=16.63982358813184
    bmag=16.38716776637639

    rfreq = alopeke_defs.i_nu
    bfreq = alopeke_defs.r_nu

    # Convert magnitude to a luminosity, i.e., nu * Lnu
    r_flux = 3631.0 * 1.0e-23 * 10**(-0.4 * rmag)
    b_flux = 3631.0 * 1.0e-23 * 10**(-0.4 * bmag)

    r_flux = r_flux * units.erg / units.s / units.Hz / (units.cm)**2
    b_flux = b_flux * units.erg / units.s / units.Hz / (units.cm)**2

    r_lum = r_flux * rfreq * 4 * np.pi * alopeke_defs.distance**2
    b_lum = b_flux * bfreq * 4 * np.pi * alopeke_defs.distance**2

    r_lum = r_lum.to(units.erg/units.second).value
    b_lum = b_lum.to(units.erg/units.second).value

    for i,modE in enumerate(Evals):
        for j,modn in enumerate(nvals):
            nval = 10**modn
            Eval = 10**modE

            model_lum_r = luminosity(rfreq.to('Hz').value,0.01,E43=Eval,nH=nval)
            model_lum_b = luminosity(bfreq.to('Hz').value,0.01,E43=Eval,nH=nval)

            if model_lum_r > r_lum or model_lum_b > b_lum:
                grid[j,i] = 1

    return(Evals, nvals, grid)

def summary_dict(date, calc_upper=False, t_rel=59145):

    summary = {}

    for camera in ['blue', 'red']: 
        summary[camera] = {}
        data_dict = alopeke_utils2.load_camera(camera, date)

        # Extinction
        EBV, A = dust_extinction(camera)
        summary[camera]['A'] = {}
        summary[camera]['A']['variable'] = 'A'
        summary[camera]['A']['value'] = 2.5*np.log10(A)  # Magnitudes! 
        summary[camera]['A']['vformat'] = '{:0.1f}'
        summary[camera]['A']['error'] = ''
        summary[camera]['A']['eformat'] = '{}'
        summary[camera]['A']['units'] = 'mag'
        summary[camera]['A']['desc'] = 'Galactic extinction'

        # Star 1
        summary[camera]['C_1'] = {}
        summary[camera]['C_1']['variable'] = '\\cstar'
        summary[camera]['C_1']['value'] = data_dict['Gauss_star1'][0]
        summary[camera]['C_1']['vformat'] = '{:0.1f}'
        summary[camera]['C_1']['error'] = data_dict['Gauss_star1'][1] / np.sqrt(len(data_dict['C_star1']))
        summary[camera]['C_1']['eformat'] = '{:0.3f}'
        summary[camera]['C_1']['units'] = '\\cunit'
        summary[camera]['C_1']['desc'] = 'Mean count rate of \\refstar'

        # counts at galaxy, away from FRB event
        summary[camera]['C_gal'] = {}
        summary[camera]['C_gal']['variable'] = '\\ctfrb'
        summary[camera]['C_gal']['value'] = data_dict['Gauss_gal'][0]
        summary[camera]['C_gal']['vformat'] = '{:0.1f}'
        summary[camera]['C_gal']['error'] = '' # We report the RMS not an error
        summary[camera]['C_gal']['eformat'] = '{}'
        summary[camera]['C_gal']['units'] = '\\cunit'
        summary[camera]['C_gal']['desc'] = 'Mean count rate of galaxy at FRB location'

        summary[camera]['sC_gal'] = {}
        summary[camera]['sC_gal']['variable'] = '\\sgal'
        summary[camera]['sC_gal']['value'] = data_dict['Gauss_gal'][1]
        summary[camera]['sC_gal']['vformat'] = '{:0.1f}'
        summary[camera]['sC_gal']['error'] = '' # We report the RMS not an error
        summary[camera]['sC_gal']['eformat'] = '{}'
        summary[camera]['sC_gal']['units'] = '\\cunit'
        summary[camera]['sC_gal']['desc'] = 'RMS in count rate of galaxy at FRB location'

        # Drift in C_gal
        p, V = time_evolution(camera, t_rel=t_rel)
        summary[camera]['dCdt'] = {}
        summary[camera]['dCdt']['variable'] = '$d\\ctfrb/dt$'
        summary[camera]['dCdt']['value'] = p[0]*1e3
        summary[camera]['dCdt']['vformat'] = '{:0.1f}'
        summary[camera]['dCdt']['error'] = np.sqrt(V[0,0])*1e3
        summary[camera]['dCdt']['eformat'] = '{:0.2f}'
        summary[camera]['dCdt']['units'] = '$10^{-3}$ counts s$^{-1}$'
        summary[camera]['dCdt']['desc'] = 'Drift in count rate during the observations'

        # Maximum C during FRB event, including galaxy
        summary[camera]['mxFRB'] = {}
        summary[camera]['mxFRB']['variable'] = '\\mxfrb'
        summary[camera]['mxFRB']['value'] = np.max(data_dict['C_FRB'])
        summary[camera]['mxFRB']['vformat'] = '{:0.1f}'
        summary[camera]['mxFRB']['error'] = '' # We report the RMS not an error
        summary[camera]['mxFRB']['eformat'] = '{}'
        summary[camera]['mxFRB']['units'] = '\\cunit'
        summary[camera]['mxFRB']['desc'] = 'Maximum count rate at FRB location during event'

        # Upper limit
        if calc_upper:
            _, _, upper_FRB = upper_limit(date, camera)
            summary[camera]['uppFRB'] = {}
            summary[camera]['uppFRB']['variable'] = '\\umufrb'
            summary[camera]['uppFRB']['value'] = upper_FRB
            summary[camera]['uppFRB']['vformat'] = '{:0.1f}'
            summary[camera]['uppFRB']['error'] = '' # We report the RMS not an error
            summary[camera]['uppFRB']['eformat'] = '{}'
            summary[camera]['uppFRB']['units'] = '\\cunit'
            summary[camera]['uppFRB']['desc'] = 'Upper limit (99.9\%) on count rate of FRB emission'

        # Fluences
        fluence, mag_frb = calc_fluence(camera)
        summary[camera]['uppFlu'] = {}
        summary[camera]['mag_lim'] = {}
        if camera == 'red':
            summary[camera]['uppFlu']['variable'] = '\\iflufrb'
        else:
            summary[camera]['uppFlu']['variable'] = '\\rflufrb'
        summary[camera]['uppFlu']['value'] = fluence.to('uJy s').value
        summary[camera]['uppFlu']['vformat'] = '{:0.3f}'
        summary[camera]['uppFlu']['error'] = '' # We report upper limit
        summary[camera]['uppFlu']['eformat'] = '{}'
        summary[camera]['uppFlu']['units'] = '\\fluunit'
        summary[camera]['uppFlu']['desc'] = 'Upper limit (99.9\%) to the fluence'
        summary[camera]['mag_lim']['value'] = mag_frb
    # Return
    return summary


if __name__ == '__main__':
    camera = ['red', 'blue']

    date = sys.argv[1]
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]

    t = Time(year+'-'+month+'-'+day)
    t_rel=t.mjd
    for cam in camera:
        time_obs(date, cam)
        values_around_frb(date, cam)
        sensitivity_function(date, cam)
        calc_fluence(cam, date)
        time_evolution(date, cam, t_rel=t_rel)
        upper_limit(date, cam)
        prob_of_max(date, cam)
        total_time(cam, date)
        dust_extinction(cam)
