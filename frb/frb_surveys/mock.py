""" Generate a set of fake FRBs for testing """

import numpy as np

from scipy.interpolate import interp1d
import pandas

from astropy import units
from astropy.coordinates import SkyCoord, match_coordinates_sky

from zdm.chime import grids

from frb.frb_surveys import chime
from frb.galaxies import hosts
from frb.defs import frb_cosmo

from IPython import embed

def frbs_for_chime(NFRB:int=10000, cut_fluence:float=5.):
    """
    Generate a mock catalog of Fast Radio Bursts (FRBs) for the CHIME survey.

    Parameters:
    - NFRB (int): Number of FRBs to generate (default: 10000)
    - cut_fluence (float): Minimum fluence threshold for FRB selection (default: 5.0)
        in mJy

    Returns:
    - df (pandas.DataFrame): DataFrame containing the mock FRB catalog with the following columns:
        - 'DMex': Extragalactic Dispersion Measure (DM) values
        - 'DM_cosmic': Cosmic DM values
        - 'DM_host': Host galaxy DM values
        - 'z': Redshift values
        - 'M_r': Absolute magnitude values
        - 'm_r': Apparent magnitude values
    """
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

    # Return
    return df


def lognorm_pdf(DM:np.ndarray, logs:tuple=None):
    """
    Calculate the probability density function (PDF) of a log-normal distribution.

    Args:
        DM (np.ndarray): Array of values representing the dispersion measure.
        logs (tuple, optional): Tuple containing the log-mean and log-standard deviation of the log-normal distribution. If not provided, default values are used.

    Returns:
        np.ndarray: Array of values representing the PDF of the log-normal distribution.

    """
    if logs is None:
        logmean = 2.16 / 0.4342944619
        logsigma = 0.51 / 0.4342944619
    else:
        logmean, logsigma = logs
    
    logDM = np.log(DM)
    
    norm = (2.0 * np.pi) ** -0.5 / DM / logsigma
    return norm * np.exp(-0.5 * ((logDM - logmean) / logsigma) ** 2)

def frbs_in_hosts(frb_tbl:pandas.DataFrame,
                  galaxy_df:pandas.DataFrame,
                  localization:tuple,
                  scale:float=2., # half-light
                  trim_catalog:units.Quantity=1*units.arcmin,
                  debug:bool=False,
                  plots:bool=False):
    """ Genereate a set of FRBs associated to galaxies
    in the input galaxy catalog
    
    With random placement in the galaxy
    and a random location error

    A pandas table is generated and saved to disk with columns:
        - ra (deg) FRB RA as observed (i.e. allows for localization error)
        - dec (deg) FRB Dec as observed (i.e. allows for localization error)
        - true_ra (deg) True FRB RA in the galaxy
        - true_dec (deg) True FRB Dec in the galaxy
        - gal_ID (int) ID of galaxy in input catalog
        - gal_off (arcsec) Offset of true FRB coord from galaxy center (true)
        - mag (float) Galaxy magnitude (m_r)
        - half_light (arcsec)
        - loc_off (arcsec) Offset of observed FRB coord from its true position

    Args:
        frb_tbl (pandas.DataFrame): Table of FRB information
            Required columns:
                m_r
                ID (int)
        galaxy_df (pandas.DataFrame):
            Galaxy catalog.  Required columns:
                ra
                dec
                mag_best
                half_light
                ID (int)
        localization (tuple): (sig_a, sig_b, PA)
        scale (float, optional): Defaults to 2..
            Scale factor for the galaxy half-light radius
        trim_catalog (units.Quantity):
            Trim edges off the catalog so that we can maintain
            PATH analysis.  Default = 1*units.arcmin
        debug (bool, optional): Defaults to False.
            Debugging mode
        plots (bool, optional): Defaults to False.
            Generate plots

    Returns:
        pandas.DataFrame: Table of FRB/Host information
    """
    # Prep
    nsample = len(frb_tbl['m_r'])

    ra_deg = galaxy_df.ra.values * units.deg
    dec_deg = galaxy_df.dec.values * units.deg
    ra_min, ra_max = ra_deg.min(), ra_deg.max()
    dec_min, dec_max = dec_deg.min(), dec_deg.max()

    # Trim
    cut_ra = (ra_deg > (ra_min + trim_catalog)) & (
        ra_deg < (ra_max - trim_catalog))
    cut_dec = (dec_deg > (dec_min + trim_catalog)) & (
        dec_deg < (dec_max - trim_catalog))

    # Cut me
    cuts = cut_dec & cut_ra
    galaxy_cut = galaxy_df[cuts]

    # Choose the galaxies
    fake_coords = SkyCoord(ra=np.ones(nsample),
                           dec=frb_tbl['m_r'], unit='deg')
    fake_galaxy = SkyCoord(ra=np.ones(len(galaxy_cut)),
                           dec=galaxy_cut.mag_best.values,
                           unit='deg')

    # Prep for matching
    galaxy_flag = np.ones(len(galaxy_cut), dtype=bool)
    galaxy_flag_idx = np.arange(len(galaxy_cut))
    galaxy_idx = galaxy_cut.index.values.copy()
    frb_idx = -1*np.ones(len(fake_coords), dtype=int)

    # Match
    while(np.any(frb_idx < 0)):

        print(f"Remaining: {np.sum(frb_idx < 0)}")
        print(f"Brightest: {np.min(fake_coords.dec)}")

        # Sub me
        sub_fake_coords = fake_coords[frb_idx < 0]
        sub_frb_idx = np.where(frb_idx < 0)[0]

        sub_fake_galaxy = fake_galaxy[galaxy_flag]
        sub_galaxy_idx = galaxy_idx[galaxy_flag] # Index in the full galaxy table
        sub_galaxy_flag_idx = galaxy_flag_idx[galaxy_flag] # Index for the flagging

        # Ran out of bright ones?
        if np.max(sub_fake_coords.dec.deg) < np.min(sub_fake_galaxy.dec.deg):
            srt_galaxy = np.argsort(sub_fake_galaxy.dec.deg)
            srt_frb = np.argsort(sub_fake_coords.dec.deg)
            # Set
            frb_idx[sub_frb_idx[srt_frb]] = sub_galaxy_idx[srt_galaxy[:len(srt_frb)]]
            assert np.all(frb_idx >= 0)
            mag_bright_cut = sub_fake_galaxy.dec.deg[srt_galaxy][len(sub_fake_coords)]
            print(f'Ran out of bright ones at {mag_bright_cut}')
            #embed(header='monte_carlo.py: 153')
            break


        print(f"Min: {sub_fake_coords.dec.min()}")
        print(f"Max: {sub_fake_coords.dec.max()}")

        # Match
        idx, d2d, _ = match_coordinates_sky(
            sub_fake_coords, sub_fake_galaxy,
            nthneighbor=1)

        # Worst case
        imx = np.argmax(d2d)
        #print(f'Max: {sub_fake_coords[imx]}')
        print(f'sep = {d2d[imx]}')

        # Take a cosmo galaxy only once
        uni, uni_idx = np.unique(idx, return_index=True)

        frb_idx[sub_frb_idx[uni_idx]] = sub_galaxy_idx[uni]
        galaxy_flag[sub_galaxy_flag_idx[uni]] = False

        #if debug:
        #    imn = np.argmin(fake_coords.dec)
        #    embed(header='monte_carlo.py: 116')

    # Generating the FRB coordinates
    galaxy_sample = galaxy_cut.loc[frb_idx]
    galaxy_coords = SkyCoord(ra=galaxy_sample.ra.values,
                             dec=galaxy_sample.dec.values,
                             unit='deg')

    # Offset the FRB in the galaxy
    theta_max = galaxy_sample.half_light.values / scale
    randn = np.random.normal(size=10*nsample)
    gd = np.abs(randn) < (6.*scale)
    randn = randn[gd][0:nsample]
    galaxy_offset = randn * theta_max * units.arcsec
    gal_pa = np.random.uniform(size=nsample, low=0., high=360.)

    print("Offsetting FRBs in galaxy...")
    frb_coord = [coord.directional_offset_by(
        gal_pa[kk]*units.deg, galaxy_offset[kk]) 
                 for kk, coord in enumerate(galaxy_coords)]

    true_frb_coord = frb_coord.copy()

    # Offset by Localization
    randn = np.random.normal(size=10*nsample)
    gd = np.abs(randn) < 3. # limit to 3 sigma
    randn = randn[gd]

    # TODO -- Make sure this is right
    a_off = randn[0:nsample] * localization[0] * units.arcsec
    b_off = randn[nsample:2*nsample] * localization[1] * units.arcsec
    #pa = np.arctan2(decoff, raoff) * 180./np.pi - 90.
    local_offset = np.sqrt(a_off**2 + b_off**2) 

    print("Offsetting FRB by localization...")
    frb_coord = [coord.directional_offset_by(
        localization[2]*units.deg, a_off[kk])
                 for kk, coord in enumerate(frb_coord)]
    frb_coord = [coord.directional_offset_by(
        (localization[2]+90)*units.deg, b_off[kk])
                 for kk, coord in enumerate(frb_coord)]

    # Write to disk
    df = pandas.DataFrame()
    df['ra'] = [coord.ra.deg for coord in frb_coord]
    df['dec'] = [coord.dec.deg for coord in frb_coord]
    df['true_ra'] = [coord.ra.deg for coord in true_frb_coord]
    df['true_dec'] = [coord.dec.deg for coord in true_frb_coord]
    df['gal_ID'] = galaxy_sample.ID.values
    df['gal_off'] = galaxy_offset.value # arcsec
    df['mag'] = galaxy_sample.mag_best.values
    df['half_light'] = galaxy_sample.half_light.values
    df['loc_off'] = local_offset.value # arcsec

    # Return
    return df
