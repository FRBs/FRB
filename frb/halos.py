""" Module for DM Halo calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

import warnings

from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from astropy import units
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy import constants

def build_grid(z_FRB=1., ntrial=10, seed=12345, Mlow=1e10, r_max=2., outfile=None, dz_box = 0.1):

    # Hidden here to avoid it becoming a dependency
    import aemHMF
    from pyigm.cgm.models import ModifiedNFW

    # #############3
    # Parameters
    #
    Mhigh = 1e16  # Msun
    Mhigh = 1e16  # Msun
    # mNFW
    y0 = 2.
    alpha = 2.
    fb = 0.75  # constant for all halos (for now)

    warnings.warn("Need to do concentration properly!")
    cgm = ModifiedNFW(alpha=alpha, y0=y0, f_hot=fb)

    # Random numbers
    rstate = np.random.RandomState(seed)

    # Setup HMF
    # https://github.com/astropy/astropy/blob/master/astropy/cosmology/parameters.py
    sigma8 = 0.8159
    ns = 0.9667
    Neff = 3.046
    cosmo_dict = {"om":cosmo.Om0,"ob":cosmo.Ob0,"ol":1.-cosmo.Om0,"ok":0.0,
              "h":cosmo.h,"s8":sigma8,"ns":ns,"w0":-1.0,"Neff":Neff} # "wa":0.0 is assumed internally
    hmf = aemHMF.Aemulus_HMF()
    hmf.set_cosmology(cosmo_dict)

    # Boxes
    nbox = int(z_FRB / dz_box)
    dz_grid = 0.01
    nz = int(z_FRB / dz_grid)
    dX = int(np.sqrt(ntrial))+1
    #
    npad = 6 # Mpc
    base_l = 2*dX + npad
    warnings.warn("Worry about being big enough given cMpc vs pMpc")

    DM_grid = np.zeros((ntrial,nz))

    # Spline distance to z
    D_max = cosmo.comoving_distance(z_FRB)
    D_val = np.linspace(1e-3,D_max.value,200) # IS THIS FINE ENOUGH?
    z_val = np.array([z_at_value(cosmo.comoving_distance, iz) for iz in D_val*units.Mpc])
    D_to_z = IUS(D_val, z_val)

    # Loop me
    prev_zbox = 0.
    for ss in range(nbox):
        zbox = ss*dz_box + dz_box/2.
        print('zbox = {}'.format(zbox))
        a = 1./(1.0 + zbox) # Scale factor
        # Mass function
        M = np.logspace(np.log10(Mlow), np.log10(Mhigh), num=1000)
        lM = np.log(M)
        dndlM = np.array([hmf.dndlM(Mi, a) for Mi in M])
        n_spl = IUS(lM, dndlM)
        cum_n = np.array([n_spl.integral(np.log(Mlow), ilM) for ilM in lM])
        ncum_n = cum_n/cum_n[-1]
        try:
            mhalo_spl = IUS(ncum_n, lM)
        except ValueError:
            # Kludge me
            print("REDUCING Mhigh by 2x")
            Mhigh /= 2.
            M = np.logspace(np.log10(Mlow), np.log10(Mhigh), num=1000)
            lM = np.log(M)
            dndlM = np.array([hmf.dndlM(Mi, a) for Mi in M])
            n_spl = IUS(lM, dndlM)
            cum_n = np.array([n_spl.integral(np.log(Mlow), ilM) for ilM in lM])
            ncum_n = cum_n/cum_n[-1]
            #
            mhalo_spl = IUS(ncum_n, lM)


        # Volume -- Box with base l = 2Mpc
        D_zn = cosmo.comoving_distance(zbox + dz_box/2.) # Full box
        D_zp = cosmo.comoving_distance(ss*dz_box) # Previous
        D_z = D_zn - D_zp
        V = D_z * (base_l*units.Mpc)**2

        # Average N_halo
        avg_n = hmf.n_bin(Mlow, Mhigh, a) * cosmo.h**3 * units.Mpc**-3
        avg_N = (V * avg_n).value

        # Assume Gaussian stats for number of halos
        N_halo = int(np.round(avg_N + np.sqrt(avg_N)*rstate.randn(1)))
        # Random masses
        randM = rstate.random_sample(N_halo)
        rM = np.exp(mhalo_spl(randM))
        # r200
        r200 = (((3*rM*units.M_sun.cgs) / (4*np.pi*200*cosmo.critical_density(zbox)))**(1/3)).to('kpc')
        # Random locations (X,Y,Z)
        X_c = rstate.random_sample(N_halo)*base_l # Mpc
        Y_c = rstate.random_sample(N_halo)*base_l # Mpc
        Z_c = (rstate.random_sample(N_halo)*D_z.to('Mpc') + D_zp).value

        # Redshifts
        z_ran = D_to_z(Z_c)

        # Random loop
        all_DMs = []
        for itrial in range(ntrial):
            # X,Y trial
            X_trial = npad//2 + (2*itrial%dX)  # Step by 2Mpc
            Y_trial = npad//2 + 2*itrial // dX
            # Impact parameters
            try:
                R_com = np.sqrt((X_c-X_trial)**2 + (Y_c-Y_trial)**2)  # Mpc
            except:
                pdb.set_trace()
            R_phys = R_com * 1000. / (1+z_ran) * units.kpc
            # Cut
            intersect = R_phys < r_max*r200
            #print("We hit {} halos".format(np.sum(intersect)))
            if not np.any(intersect):
                continue
            # Loop -- FIND A WAY TO SPEED THIS UP!
            DMs = []
            for iobj in np.where(intersect)[0]:
                cgm.log_Mhalo=np.log10(rM[iobj])
                cgm.M_halo = 10.**cgm.log_Mhalo * constants.M_sun.cgs
                cgm.z = zbox # To be consistent with above;  close enough
                cgm.setup_param(cosmo=cosmo)
                # DM
                DM = cgm.Ne_Rperp(R_phys[iobj], rmax=r_max, add_units=False)
                DMs.append(DM)
            # Save em
            iz = (z_ran[intersect]/dz_grid).astype(int)
            DM_grid[itrial,iz] += DMs
    # Write
    if outfile is not None:
        np.save(outfile, DM_grid, allow_pickle=False)

    return DM_grid

# Command line execution
if __name__ == '__main__':
    build_grid(outfile='z1_mNFW_10000', ntrial=10000)
    #build_grid(outfile='test', ntrial=10)

