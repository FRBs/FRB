import numpy as np
from astropy.time import Time
from astropy import units
import json

from frb.frb import FRB
frb180916 = FRB.by_name('FRB20180916B')

# CHIME
# Tendulkar (Feb 11, 2021): 
# This pulse arrival time is topocentric at 400 MHz.
# In the database I see 2020-10-23 07:48:30.777667 UTC+00:00 
# The uncertainty is 1 ms. 

chime_time_arrival = {
    '20201023':Time('2020-10-23T07:48:30.777667'),
    '20220908':Time('2022-09-08T10:53:26.888596'),
}
ata_chime = 1*units.s/1000.  # see comment above

#  Alopeke absolute time accuracy (from email Feb 8, 2021 Nic Scott)
#  Nic: The major contributor in this uncertainty is thought to be the variable lag
#  between the computer receipt from the NTP server and the triggering of the
#  cameras.
ata_alopeke = 163*units.s/1000.

# final absolute time accuracy (1 sigma)
ata = np.sqrt(ata_alopeke**2 + ata_chime**2)  # quadrature of absolute uncertainties

# According to the calculation of Kilpatrick, from 400 MHz to optical
# frequencies, a delay of 9.1s should be applied to the radio.
# Also account for light-travel time for CHIME to Alopeke (14.8ms)
chime_time_arrival_optical = {}
FRB_time = {}
mjd_low = {}
mjd_high = {}
for key in chime_time_arrival.keys():
    chime_time_arrival_optical[key]=chime_time_arrival[key] - 9.08*units.s - 0.0148*units.s
    FRB_time[key]=chime_time_arrival_optical[key]
    mjd_low[key] = (FRB_time[key]-ata).mjd
    mjd_high[key] = (FRB_time[key]+ata).mjd

# Gains 
gain_red = 5.61
gain_blue = 5.54
EM = 1000

# Center of band
XSDSS_r = 620. # nm
XSDSS_i = 765. # nm

# Exposure time
#dt_alopeke = 11.6 * (1e-3 * units.s)
dt_alopeke = 10.419 * (1e-3 * units.s)  # on-sky exposure time

# Photometry Star-1
# From Charlie Kilpatrick #frb180916-speckle channel (5 deb 2021)
# r=15.7481+/-0.00017 (PS1)
# i=15.1387+/-0.00237 (PS1)
r_1 = 15.7481  # Panstarrs-1 magnitude
i_1 = 15.1387  # Panstarrs-1 magnitude

r_2 = 17.1016  # Panstarrs-1 magnitude
i_2 = 16.6247  # Panstarrs-1 magnitude

redshift = 0.0337

r_nu = (2.998e18 * units.angstrom/units.second) / (6231 * units.angstrom)
i_nu = (2.998e18 * units.angstrom/units.second) / (7625 * units.angstrom)

distance = 150.0e6 * units.parsec

# Radio data from CHIME (email from Emmanuel Fonseca, 2023-04-21)
burst1_file = '../../Data/results_R3/results_fitburst_139459007.json'
burst2_file = '../../Data/results_R3/results_fitburst_244202260.json'

burst1_data = {}
burst2_data = {}

with open(burst1_file, 'r') as f:
    burst1_data = json.load(f)
with open(burst2_file, 'r') as f:
    burst2_data = json.load(f)

radio_data = {
    '20201023':burst1_data,
    '20220908':burst2_data,
}
