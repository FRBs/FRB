""" Module to correct pulsar and FRB DMs for the MW ISM """

import sys
from ne2001 import ne_io, density #ne2001 ism model
import pygedm #ymw ism model
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)
ne = density.ElectronDensity()

def find_delta_dm(transient_type,transient_data,ism_model,b_val,mc_deg,save_df=True):
    """
    Find pulsar/FRB DMs corrected for by the MW ISM DM and remove observations in complex DM regions.
    Returns array of DMs
    To use this code, please download the FRB catalogue from http://frbcat.org/ [Petroff et al. 2017]
    and the ATNF pulsar catalogue [Manchester et al. 2005].
    Check columns match the template used here.
    
    Arguments:
        transient_type (str):
            Accepts 'frb' or 'pulsar'.
        transient_data (str):
            Path to data (in .csv format).
        ism_model (str):
            Model used to calculated the MW halo DM.
            Accepts 'ymw16'  [Yao et al. 2017] or 'ne2001' [Cordes & Lazio 2003].
        b_val (int):
            Galactic latitude considered (b>b_val, b<-b_val).
        mc_deg (int):
            Number of degrees from Magellanic clouds within which transients are removed.
        save_df (str, optional):
            Save transient DMs and coords to csv.
    
    Outputs:
    """
    # Sort data and get coords
    if transient_type=='frb':
        transcat_df = pd.read_csv(transient_data, skiprows=1, usecols= [0,5,6,7], names=['Name','l','b','dm']) 
        transcat_df['dm'] = transcat_df['dm'].str.split('&').str[0].astype(float).values
        coords = SkyCoord(l=transcat_df['l'], b=transcat_df['b'], unit=(u.degree),frame=Galactic)
    elif transient_type=='pulsar':
        transcat_df = pd.read_csv(transient_data, skiprows=2, usecols = [1,2,3,9,10], names=['Name','Pref','dm','RAJD','DECJD'])
        transcat_df = transcat_df[~transcat_df['dm'].str.contains('*', regex=False)].reset_index(drop=True)
        transcat_df['dm'] = transcat_df['dm'].astype(float)
        c_icrs = SkyCoord(ra=transcat_df['RAJD'], dec=transcat_df['DECJD'], unit=(u.degree), frame='icrs')
        transcat_df['l'] = pd.DataFrame(c_icrs.galactic.l.value)
        transcat_df['b'] = pd.DataFrame(c_icrs.galactic.b.value)
        coords = SkyCoord(l=transcat_df['l'], b=transcat_df['b'], unit=(u.degree),frame=Galactic)

    # Find transients in line of sight of MCs
    logging.info('Removing transients near Magellanic clouds...')
    # LMC
    lmc_distance = 50*u.kpc
    lmc_coord = SkyCoord('J052334.6-694522',unit=(u.hourangle, u.deg),distance=lmc_distance)
    close_to_lmc = lmc_coord.separation(coords) < mc_deg*u.deg
    lmc_trans = list(transcat_df[close_to_lmc]['Name'])
    # SMC
    smc_distance = 61*u.kpc
    smc_coord = SkyCoord('J005238.0-724801',unit=(u.hourangle, u.deg),distance=smc_distance)
    close_to_smc = smc_coord.separation(coords) < mc_deg*u.deg
    smc_trans = list(transcat_df[close_to_smc]['Name'])
    
    transcat_df = transcat_df[~transcat_df['Name'].isin(lmc_trans)].reset_index(drop=True)
    transcat_df = transcat_df[~transcat_df['Name'].isin(smc_trans)].reset_index(drop=True)
    if transient_type=='pulsar':
        transcat_df = transcat_df[~transcat_df['Pref'].str.contains('mfl+06', regex=False)].reset_index(drop=True)
    elif transient_type=='frb':
        pass
    
    # Remove transients with low Galactic lattitudes
    logging.info('Removing transients with low Galactic lattitudes...')
    transcat_df = pd.concat([transcat_df[transcat_df.b > b_val], transcat_df[transcat_df.b < -b_val]], ignore_index=True)

    # ISM model
    logging.info('Correcting transient DMs for ISM...')
    trans_ism = []
    if ism_model=='ymw16':
        for i in range(len(transcat_df['dm'])):
            trans_ism_ = pygedm.dist_to_dm(transcat_df['l'].iloc[i], transcat_df['b'].iloc[i], 100000)[0].value
            trans_ism = np.append(trans_ism,trans_ism_)
    elif ism_model=='ne2001':
        for i in range(len(transcat_df['dm'])):
            trans_ism_ = ne.DM(transcat_df['l'].iloc[i], transcat_df['b'].iloc[i], 100.).value
            trans_ism = np.append(trans_ism,trans_ism_)

    transcat_df['trans_ism'] = pd.DataFrame(trans_ism)
    transcat_df['deltaDM'] = pd.DataFrame(transcat_df['dm']-transcat_df['trans_ism'])

    if save_df==True:
        transcat_df.to_csv('transient_data/'+transient_type+'cat_df_'+ism_model+'_'+str(int(b_val))+'.csv')
        logging.info('Transient data saved to csv.')
    else:
        pass

    return np.array(transcat_df['deltaDM'])

