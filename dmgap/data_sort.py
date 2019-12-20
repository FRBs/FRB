import sys
sys.path.insert(1, '/eplatts_UCSC_Server/dm_gap/ne2001-master/src')
import numpy as np
import pandas as pd
from ne2001 import ne_io
from ne2001 import density
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic

b_val = 20. #latitude considered (b>b_val, b<-b_val)

ne = density.ElectronDensity()

psrcat_df = pd.read_csv('transient_data/psrcat.csv', skiprows=2, usecols = [1,2,3,9,10], names=['Name','Pref','dm','RAJD','DECJD']) 

# FRB DATA
frbcat_df = pd.read_csv('transient_data/frbcat_20191111.csv', skiprows=1, usecols= [0,5,6,7], names=['Name','l','b','dm']) 
frbcat_df['dm'] = frbcat_df['dm'].str.split('&').str[0].astype(float).values

# Find FRBs in line of sight of MCs
coords_frb = SkyCoord(l=frbcat_df['l'], b=frbcat_df['b'], unit=(u.degree),frame=Galactic)
mfl = psrcat_df['Pref'] == 'mfl+06'
lmc_distance = 50*u.kpc
lmc_coord = SkyCoord('J052334.6-694522',unit=(u.hourangle, u.deg),distance=lmc_distance)
close_to_lmc = lmc_coord.separation(coords_frb) < 3*u.deg
lmc_frb = list(frbcat_df[close_to_lmc]['Name'])
# SMC
smc_distance = 61*u.kpc
smc_coord = SkyCoord('J005238.0-724801',unit=(u.hourangle, u.deg),distance=smc_distance)
# smc_coord_ = smc_coord.separation(coords_frb[mfl]).to('deg').value
close_to_smc = smc_coord.separation(coords_frb) < 3*u.deg
smc_frb = list(frbcat_df[close_to_smc]['Name'])
frbcat_df = frbcat_df[~frbcat_df['Name'].isin(lmc_frb)].reset_index(drop=True)
frbcat_df = frbcat_df[~frbcat_df['Name'].isin(smc_frb)].reset_index(drop=True)

frbcat_df = pd.concat([frbcat_df[frbcat_df.b > b_val], frbcat_df[frbcat_df.b < -b_val]], ignore_index=True)
print('FRB datasize is:', len(frbcat_df))
# FRB ne2001
frb_dmmax = []
for i in range(len(frbcat_df['dm'])):
    frb_dmmax_ = ne.DM(frbcat_df['l'].iloc[i], frbcat_df['b'].iloc[i], 100.).value
    frb_dmmax = np.append(frb_dmmax,frb_dmmax_)

frbcat_df['dmmax'] = pd.DataFrame(frb_dmmax)
frbcat_df['dmdiff'] = pd.DataFrame(frbcat_df['dm']-frbcat_df['dmmax'])
frbcat_df.to_csv('transient_data/frbcat_df.csv')
print('FRB data saved')

# PULSAR DATA
psrcat_df = pd.read_csv('transient_data/psrcat.csv', skiprows=2, usecols = [1,2,3,9,10], names=['Name','Pref','dm','RAJD','DECJD'])
psrcat_df = psrcat_df[~psrcat_df['dm'].str.contains('*', regex=False)].reset_index(drop=True)
psrcat_df['dm'] = psrcat_df['dm'].astype(float)
print(len(psrcat_df))
coords = SkyCoord(ra=psrcat_df['RAJD'], dec=psrcat_df['DECJD'], unit=(u.degree))

# Find pulsars within Magellanic clouds
# LMC
mfl = psrcat_df['Pref'] == 'mfl+06'
lmc_distance = 50*u.kpc
lmc_coord = SkyCoord('J052334.6-694522',unit=(u.hourangle, u.deg),distance=lmc_distance)
lmc_coord_ = lmc_coord.separation(coords[mfl]).to('deg').value
close_to_lmc = lmc_coord.separation(coords) < 3*u.deg
lmc_pulsars = list(psrcat_df[close_to_lmc]['Name'])
# SMC
smc_distance = 61*u.kpc
smc_coord = SkyCoord('J005238.0-724801',unit=(u.hourangle, u.deg),distance=smc_distance)
smc_coord_ = smc_coord.separation(coords[mfl]).to('deg').value
close_to_smc = smc_coord.separation(coords) < 3*u.deg
smc_pulsars = list(psrcat_df[close_to_smc]['Name'])

# Remove pulsars in/near MCs
psrcat_df = psrcat_df[~psrcat_df['Name'].isin(lmc_pulsars)].reset_index(drop=True)
psrcat_df = psrcat_df[~psrcat_df['Name'].isin(smc_pulsars)].reset_index(drop=True)
psrcat_df = psrcat_df[~psrcat_df['Pref'].str.contains('mfl+06', regex=False)].reset_index(drop=True)
print(len(psrcat_df))
c_icrs = SkyCoord(ra=psrcat_df['RAJD'], dec=psrcat_df['DECJD'], unit=(u.degree), frame='icrs')
psrcat_df['l'] = pd.DataFrame(c_icrs.galactic.l.value)
psrcat_df['b'] = pd.DataFrame(c_icrs.galactic.b.value)
psrcat_df = pd.concat([psrcat_df[psrcat_df.b > b_val], psrcat_df[psrcat_df.b < -b_val]], ignore_index=True)
print(len(psrcat_df))
# Pulsar ne2001
psr_dmmax = []
for i in range(len(psrcat_df['dm'])):
    psr_dmmax_ = ne.DM(psrcat_df['l'].iloc[i], psrcat_df['b'].iloc[i], 100.).value
    psr_dmmax = np.append(psr_dmmax,psr_dmmax_)

psrcat_df['dmmax'] = pd.DataFrame(psr_dmmax)
psrcat_df['dmdiff'] = pd.DataFrame(psrcat_df['dm']-psrcat_df['dmmax'])
psrcat_df = psrcat_df.drop('RAJD', axis=1)
psrcat_df = psrcat_df.drop('DECJD', axis=1)
psrcat_df.to_csv('transient_data/psrcat_df.csv')
print('Pulsar data saved')


