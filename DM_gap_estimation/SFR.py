"This calculates the SFR at a given redshift, which we use as an approximation of the FRB rate"

import numpy as np
import sklearn
from astropy.cosmology import WMAP9 as cosmo

"Redshift and conversion to Mpc"
z_stepsize = 0.0001
z_max = 1.2
z = np.arange(0,z_max,z_stepsize)
Mpc = cosmo.comoving_distance(z).value 
Mpc_stepsize = cosmo.comoving_distance(z_stepsize).value

psi = 0.015*(1+z)**2.7/(1+((1+z)/2.9)**5.6) #star formation rate density

"SFR scaled by exponential cut off"
FRB_rate_density = psi*np.exp(-z/0.5)

V_sphere = 4/3*np.pi*Mpc**3
V_shell = [y - x for x,y in zip(V_sphere,V_sphere[1:])]
FRB_distribution_in_z = FRB_rate_density[1:]*V_shell

"Normalize PDF"
FRB_distribution_SFR = FRB_distribution_in_z/z_stepsize