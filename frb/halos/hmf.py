import warnings

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from astropy import units
from astropy import constants
from astropy.cosmology import z_at_value
from astropy.table import Table

from frb.halos.models import ModifiedNFW, ICM
from frb.defs import frb_cosmo as cosmo

from IPython import embed

def init_hmf():
    """
    Initialize the Aemulus Halo Mass Function

    WARNING: This uses the original version which codes Tinker+2008
    We may refactor to use the more accurate, new version

    Returns:
        hmfe (hmf_emulator.hmf_emulator): An Aemulus halo mass function emulator.

    """
    # Hidden here to avoid it becoming a dependency
    import hmf_emulator
    # Setup HMF
    # https://github.com/astropy/astropy/blob/master/astropy/cosmology/parameters.py
    #sigma8 = 0.8159
    ns = 0.9667
    Neff = 3.046
    #cosmo_dict = {"om":cosmo.Om0,"ob":cosmo.Ob0,"ol":1.-cosmo.Om0,"ok":0.0,
    #              "h":cosmo.h,"s8":sigma8,"ns":ns,"w0":-1.0,"Neff":Neff} # "wa":0.0 is assumed internally
    cosmo_dict = {"omega_cdm":(cosmo.Om0-cosmo.Ob0)*cosmo.h**2,
                  "omega_b":cosmo.Ob0*cosmo.h**2,"ok":0.0,
                  "ln10As": 3.098, # THIS REPLACES sigma8
                  "H0":cosmo.H0.to('km/(s*Mpc)').value,
                  "n_s":ns,"w0":-1.0,"N_eff":Neff} # "wa":0.0 is assumed internally
    hmfe = hmf_emulator.hmf_emulator()
    hmfe.set_cosmology(cosmo_dict)
    # Return
    return hmfe

# Storing for use
#  Needs to come after def init_hmf()
try:
    import hmf_emulator
except:
    warnings.warn("hmf_emulator not imported.  Hope you are not intending to use the hmf.py module..")
else:
    hmfe = init_hmf()


def frac_in_halos(zvals, Mlow, Mhigh, rmax=1.):
    """
    Calculate the fraction of matter in collapsed halos
     over a mass range and at a given redshift

    Note that the fraction of DM associated with these halos
    will be scaled down by an additional factor of f_diffuse

    Requires Aemulus HMF to be installed

    Args:
        zvals: ndarray
        Mlow: float
          In h^-1 units already so this will be applied for the halo mass function
        Mhigh: float
          In h^-1 units already
        rmax: float
          Extent of the halo in units of rvir
          WARNING: This calculation assumes a single concentration for all halos

    Returns:
        ratios: ndarray
          rho_halo / rho_m
    """

    M = np.logspace(np.log10(Mlow*cosmo.h), np.log10(Mhigh*cosmo.h), num=1000)
    lM = np.log(M)

    ratios = []
    for z in zvals:
        # Setup
        #dndlM = np.array([hmfe.dndlnM(Mi, a)[0] for Mi in M])
        dndlM = M*hmfe.dndM(M, z)
        M_spl = IUS(lM, M * dndlM)

        # Integrate
        rho_tot = M_spl.integral(np.log(Mlow*cosmo.h), np.log(Mhigh*cosmo.h)) * units.M_sun / units.Mpc ** 3
        # Cosmology
        rho_M = cosmo.critical_density(z) * cosmo.Om(z)/(1+z)**3  # Tinker calculations are all mass
        ratio = (rho_tot*cosmo.h**2 / rho_M).decompose()
        #
        ratios.append(ratio)
    ratios = np.array(ratios)
    # Boost halos if extend beyond rvir (homologous in mass, but constant concentration is an approx)
    if rmax != 1.:
        #from pyigm.cgm.models import ModifiedNFW
        c = 7.7
        nfw = ModifiedNFW(c=c)
        M_ratio = nfw.fy_dm(rmax * nfw.c) / nfw.fy_dm(nfw.c)
        ratios *= M_ratio
    # Return
    return np.array(ratios)


def halo_incidence(Mlow, zFRB, radius=None, hmfe=None, 
                   Mhigh=1e16, nsample=20, cumul=False):
    """
    Calculate the (approximate) average number of 
    intersections to halos of a
    given minimum mass to a given zFRB.

    Requires Aemulus HMF to be installed

    Args:
        Mlow: float
          Mass of minimum halo in Solar masses
          The code deals with h^-1 factors so that you do not
          The minimum value is 2e10
        zFRB: float
          Redshift of the FRB
        radius: Quantity, optional
          Physical separation from the sightline for the calculation.
          The calculation will specify this radius as rvir derived from
           Mlow unless this is specified. And this rvir *will* vary with redshift
        hmfe (hmf.hmf_emulator, optional): Halo mass function emulator from Aeumulus
        Mhigh: float, optional
          Mass of maximum halo in Solar masses
        nsammple: int, optional
          Number of samplings in redshift
          20 should be enough
        cumul: bool, optional
          Return the cumulative quantities instead

    Returns:
        If cumul is False
        Navg: float
          Number of average intersections
        elif cumul is True
        zeval: ndarray
        Ncumul: ndarray
    """
    # Mlow limit
    if Mlow < 2e10:
        raise IOError("Calculations are limited to Mlow > 2e10")
    # HMF
    if hmfe is None:
        hmfe = init_hmf()
    #
    zs = np.linspace(0., zFRB, nsample)
    # Mean density
    ns = []
    for iz in zs:
        ins = hmfe.n_in_bins((Mlow * cosmo.h, Mhigh * cosmo.h), iz) 
        ns.append(ins[0]*cosmo.h**3)  # * units.Mpc**-3
    # Interpolate
    ns = units.Quantity(ns*units.Mpc**-3)
    # Radii
    if radius is None:
        rhoc = cosmo.critical_density(zs)
        #https://arxiv.org/pdf/1312.4629.pdf eq5
        q = cosmo.Ode0/(cosmo.Ode0+cosmo.Om0*(1+zs)**3)
        rhovir = (18*np.pi**2-82*q-39*q**2)*rhoc
        r200 = (((3*Mlow*constants.M_sun.cgs) / (4*np.pi*rhovir))**(1/3)).to('kpc')
    else:
        r200 = np.ones_like(zs) * radius
    # Ap
    Ap = np.pi * r200**2

    # l(X)
    loX = ((constants.c/cosmo.H0) * ns * Ap).decompose().value

    # dX
    X = cosmo.absorption_distance(zs)
    dX = X - np.roll(X,1)
    dX[0] = 0.

    # Finish
    if cumul:
        Navg = np.cumsum(loX * dX)
        return zs, Navg
    else:
        Navg = np.sum(loX * dX)
        return Navg


def build_grid(z_FRB=1., ntrial=10, seed=12345, Mlow=1e10, r_max=2., outfile=None, dz_box=0.1,
    dz_grid=0.01, f_hot=0.75, verbose=True):
    """
    Generate a universe of dark matter halos with DM measurements
    Mainly an internal function for generating useful output grids.

    Requires the Aemulus Halo Mass function

    Args:
        z_FRB: float, optional
        ntrial: int, optional
        seed: int, optional
        Mlow: float, optional
          h^-1 mass
        r_max: float, optional
          Extent of the halo in units of rvir
        outfile: str, optional
          Write
        dz_box: float, optional
          Size of the slice of the universe for each sub-calculation
        dz_grid: float, optional
          redshift spacing in the DM grid
        f_hot: float
          Fraction of the cosmic fraction of matter in diffuse gas (for DM)

    Returns:
        DM_grid: ndarray (ntrial, nz)
        halo_tbl: Table
          Table of all the halos intersected

    """
    Mhigh = 1e16  # Msun
    # mNFW
    y0 = 2.
    alpha = 2.

    warnings.warn("Ought to do concentration properly someday!")
    cgm = ModifiedNFW(alpha=alpha, y0=y0, f_hot=f_hot)
    icm = ICM()

    # Random numbers
    rstate = np.random.RandomState(seed)

    # Init HMF
    hmfe = init_hmf()

    # Boxes
    nbox = int(z_FRB / dz_box)
    nz = int(z_FRB / dz_grid)
    dX = int(np.sqrt(ntrial))+1
    #
    npad = 6 # Mpc
    base_l = 2*dX + npad
    print('L_base = {} cMpc'.format(base_l))
    warnings.warn("Worry about being big enough given cMpc vs pMpc")

    DM_grid = np.zeros((ntrial,nz))

    # Spline distance to z
    D_max = cosmo.comoving_distance(z_FRB)
    D_val = np.linspace(1e-3,D_max.value,200) # IS THIS FINE ENOUGH?
    z_val = np.array([z_at_value(cosmo.comoving_distance, iz) for iz in D_val*units.Mpc])
    D_to_z = IUS(D_val, z_val)

    # Save halo info
    #halos = [[] for i in range(ntrial)]
    halo_i, M_i, R_i, DM_i, z_i = [], [], [], [], []

    # Loop me
    prev_zbox = 0.
    #for ss in range(nbox):
    #for ss in [0]:
    for ss in [5]:
        zbox = ss*dz_box + dz_box/2.
        print('zbox = {}'.format(zbox))
        a = 1./(1.0 + zbox) # Scale factor
        # Mass function
        M = np.logspace(np.log10(Mlow*cosmo.h), np.log10(Mhigh*cosmo.h), num=1000)
        lM = np.log(M)
        dndlM = np.array([hmf.dndlM(Mi, a) for Mi in M])
        n_spl = IUS(lM, dndlM)
        cum_n = np.array([n_spl.integral(np.log(Mlow*cosmo.h), ilM) for ilM in lM])
        ncum_n = cum_n/cum_n[-1]
        # As z increases, we have numerical issues at the high mass end (they are too rare)
        try:
            mhalo_spl = IUS(ncum_n, lM)
        except ValueError:
            # Kludge me
            print("REDUCING Mhigh by 2x")
            Mhigh /= 2.
            M = np.logspace(np.log10(Mlow*cosmo.h), np.log10(Mhigh*cosmo.h), num=1000)
            lM = np.log(M)
            dndlM = np.array([hmf.dndlM(Mi, a) for Mi in M])
            n_spl = IUS(lM, dndlM)
            cum_n = np.array([n_spl.integral(np.log(Mlow*cosmo.h), ilM) for ilM in lM])
            ncum_n = cum_n/cum_n[-1]
            #
            mhalo_spl = IUS(ncum_n, lM)

        # Volume -- Box with base l = 2Mpc
        D_zn = cosmo.comoving_distance(zbox + dz_box/2.) # Full box
        D_zp = cosmo.comoving_distance(ss*dz_box) # Previous
        D_z = D_zn - D_zp
        V = D_z * (base_l*units.Mpc)**2

        # Average N_halo
        avg_n = hmf.n_bin(Mlow*cosmo.h, Mhigh*cosmo.h, a) * cosmo.h**3 * units.Mpc**-3
        avg_N = (V * avg_n).value

        # Assume Gaussian stats for number of halos
        N_halo = int(np.round(avg_N + np.sqrt(avg_N)*rstate.randn(1)))

        # Random masses
        randM = rstate.random_sample(N_halo)
        rM = np.exp(mhalo_spl(randM)) / cosmo.h

        # r200
        r200 = (((3*rM*units.M_sun.cgs) / (4*np.pi*200*cosmo.critical_density(zbox)))**(1/3)).to('kpc')

        # Random locations (X,Y,Z)
        X_c = rstate.random_sample(N_halo)*base_l # Mpc
        Y_c = rstate.random_sample(N_halo)*base_l # Mpc
        Z_c = (rstate.random_sample(N_halo)*D_z.to('Mpc') + D_zp).value

        # Check mass fraction
        if verbose:
            Mtot = np.log10(np.sum(rM))
            M_m = (cosmo.critical_density(zbox)*cosmo.Om(zbox) * V/(1+zbox)**3).to('M_sun')
            #print("N_halo: {}  avg_N: {}".format(N_halo, avg_N))
            print("z: {}  Mhalo/M_m = {}".format(zbox, 10**Mtot/M_m.value))
            print(frac_in_halos([zbox], Mlow, Mhigh))

        # Redshifts
        z_ran = D_to_z(Z_c)

        # Loop on trials
        all_DMs = []
        all_nhalo = []
        all_r200 = []
        for itrial in range(ntrial):
            # X,Y trial
            X_trial = npad//2 + (2*itrial%dX)  # Step by 2Mpc
            Y_trial = npad//2 + 2*itrial // dX
            # Impact parameters
            try:
                R_com = np.sqrt((X_c-X_trial)**2 + (Y_c-Y_trial)**2)  # Mpc
            except:
                embed()
            R_phys = R_com * 1000. / (1+z_ran) * units.kpc
            # Cut
            intersect = R_phys < r_max*r200
            print("We hit {} halos".format(np.sum(intersect)))
            all_nhalo.append(np.sum(intersect))
            if not np.any(intersect):
                all_DMs.append(0.)
                continue
            # Loop -- FIND A WAY TO SPEED THIS UP!
            DMs = []
            for iobj in np.where(intersect)[0]:
                # Init
                if rM[iobj] > 1e14: # Use ICM model
                    model = icm
                else:
                    model = cgm
                model.log_Mhalo=np.log10(rM[iobj])
                model.M_halo = 10.**model.log_Mhalo * constants.M_sun.cgs
                model.z = zbox # To be consistent with above;  should be close enough
                model.setup_param(cosmo=cosmo)
                # DM
                DM = model.Ne_Rperp(R_phys[iobj], rmax=r_max, add_units=False)/(1+model.z)
                DMs.append(DM)
                # Save halo info
                halo_i.append(itrial)
                M_i.append(model.M_halo.value)
                R_i.append(R_phys[iobj].value)
                DM_i.append(DM)
                z_i.append(z_ran[iobj])
                all_r200.append(cgm.r200.value)
            # Save em
            iz = (z_ran[intersect]/dz_grid).astype(int)
            DM_grid[itrial,iz] += DMs
            all_DMs.append(np.sum(DMs))
            #print(DMs, np.log10(rM[intersect]), R_phys[intersect])
            if (itrial % 100) == 0:
                embed()

    # Table the halos
    halo_tbl = Table()
    halo_tbl['trial'] = halo_i
    halo_tbl['M'] = M_i
    halo_tbl['R'] = R_i
    halo_tbl['DM'] = DM_i
    halo_tbl['z'] = z_i

    # Write
    if outfile is not None:
        print("Writing to {}".format(outfile))
        np.save(outfile, DM_grid, allow_pickle=False)
        halo_tbl.write(outfile+'.fits', overwrite=True)

    return DM_grid, halo_tbl

