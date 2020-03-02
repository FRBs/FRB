""" Module for basic plots related to FRB host and foreground galaxies"""
import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from IPython import embed

from pkg_resources import resource_filename

from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy import units

from frb.halos import ModifiedNFW
from frb import halos as frb_halos
from frb import igm as frb_igm

from ne2001 import density


def sub_cartoon(ax1, ax2, coord, zFRB, halos=False, host_DM=50., ymax=None,
                IGM_only=True, smin=0.1,
                M31=None, fg_halos=None, dsmx=0.05, FRB_DM=None, yscl = 0.97):
    """
    Cartoon of DM cumulative

    Plot of increasing DM from Earth to the FRB

    Args:
        ax1 (matplotlib.Axis):
            First axis.  Used for Milky Way and local group
        ax2 (matplotlib.Axis):
            Second axis.  Used for Cosmic and Host
        coord (astropy.coord.SkyCoord):
            Coordinate of the FRB used with ne2001
        zFRB (float):
            Redshift of the FRB
        halos (?, optional):
            Not implemented!
        host_DM (float):
            DM to use for the Host
        ymax (tuple or list):
            ymin, ymax values for the y-axis
        IGM_only (bool, optional):
            Use only the IGM for DM_Cosmic, i.e. ignore the presumed average halo contribution
        M31 (bool, optional):
            Include M31 in the calculation?
            NOT IMPLEMENTED RIGHT NOW
        fg_halos (dict or Table):
            Used to add to DM_IGM
            Keys must include 'z' 'DM' 'lbl'
        smin (float, optional):
            Minimum value in axis2 (Gpc)
        dsmx (float, optional): Padding on the x-axis;  Gpc
        FRB_DM (float): Observed value;  sets ymax = FRB_DM+50

    """
    if halos:
        embed()

    gcoord = coord.transform_to('galactic')
    l, b = gcoord.l.value, gcoord.b.value

    ds = []  # kpc
    DM_cumul = []

    # ISM
    ne = density.ElectronDensity()  # **PARAMS)
    for ss in np.linspace(2., 4, 5):  # log pc
        idd = 10. ** ss / 1e3  # kpc
        iDM = ne.DM(l, b, idd)
        # Append
        ds.append(idd)  # kpc
        DM_cumul.append(iDM.value)
        # print(idd, iDM)
    max_ISM = DM_cumul[-1]

    # MW
    Mhalo = np.log10(1.5e12)  # Boylan-Kolchin et al. 2013
    f_hot = 0.75  # Allows for disk + ISM
    c = 7.7
    mnfw_2 = ModifiedNFW(log_Mhalo=Mhalo, f_hot=f_hot, y0=2, alpha=2, c=c)
    # Zero out inner 10kpc
    mnfw_2.zero_inner_ne = 10.  # kpc
    params = dict(F=1., e_density=1.)
    model_ne = density.NEobject(mnfw_2.ne, **params)
    for ss in np.linspace(1., np.log10(mnfw_2.r200.value), 5):  # log kpc
        idd = 10. ** ss  # kpc
        iDM = model_ne.DM(l, b, idd).value
        # Add it in
        if idd == ds[-1]:
            DM_cumul[-1] = DM_cumul[-1] + iDM
        else:
            ds.append(idd)
            DM_cumul.append(max_ISM + iDM)

    DM_ISM_Halo = DM_cumul[-1]

    if M31 is not None:
        raise NotImplemented
        # M31
        m31 = M31()
        a,c =1,0
        x0, y0 = m31.distance.to('kpc').value, 0. # kpc (Riess, A.G., Fliri, J., & Valls - Gabaud, D. 2012, ApJ, 745, 156)
        sep = m31.coord.separation(coord)
        atan = np.arctan(sep.radian)
        b = -1 * a / atan
        M31_Rperp = np.abs(a * x0 + b * y0 + c) / np.sqrt(a ** 2 + b ** 2)  # kpc
        zval, M31_DM_cumul = m31.Ne_Rperp(M31_Rperp * u.kpc, rmax=1., cumul=True)
        # Add em in
        ds += (zval+x0).tolist()
        DM_cumul += (M31_DM_cumul+DM_ISM_Halo).tolist()

    #DM_LG = 0.
    DM_LG = DM_cumul[-1]

    # IGM
    z0 = z_at_value(cosmo.comoving_distance, 1 * units.Mpc)
    zvals = np.linspace(z0, 0.5, 50)
    dz_vals = zvals[1] - zvals[0]
    #
    DM_cosmic_cumul, zeval = frb_igm.average_DM(zvals[-1], cumul=True)
    dzeval = zeval[1] - zeval[0]
    dDM_cosmic = DM_cosmic_cumul - np.roll(DM_cosmic_cumul, 1)
    dDM_cosmic[0] = dDM_cosmic[1]
    #
    DM_interp = IUS(zeval, dDM_cosmic)

    dDM_cosm = DM_interp(zvals) * dz_vals / dzeval
    sub_DM_cosm = np.cumsum(dDM_cosm)
    f_cosm = IUS(zvals, sub_DM_cosm)
    zvals2 = np.linspace(z0, zFRB, 1000)
    DM_cosmic = f_cosm(zvals2)

    # Ignore halos?
    if IGM_only:
        #
        fhalos = frb_halos.frac_in_halos(zvals, 3e10, 1e16, rmax=1.)
        fIGM = 1. - fhalos
        #
        dDM_IGM = DM_interp(zvals) * fIGM * dz_vals / dzeval
        sub_DM_IGM = np.cumsum(dDM_IGM)
        f_IGM = IUS(zvals, sub_DM_IGM)
        DM_IGM = f_IGM(zvals2)
        #
        DM_cosmic = DM_IGM.copy()

    # Halos at last
    if fg_halos is not None:
        for z, halo_DM, lbl in zip(fg_halos['z'], fg_halos['DM'], fg_halos['lbl']):
            iz = np.argmin(np.abs(zvals2 - z))
            DM_cosmic[iz:] += halo_DM
            # Label
            d = cosmo.comoving_distance(z)
            ax1.text(d.to('Gpc').value, DM_cosmic[iz], lbl, color='black',
                     fontsize=13, ha='left', va='top')

    Dc = cosmo.comoving_distance(zvals2).to('kpc')
    ds += Dc.value.tolist()
    DM_cumul += (DM_cosmic + DM_LG).tolist()

    # Host
    if host_DM > 0.:
        ds.append(ds[-1])
        DM_cumul.append(DM_cumul[-1] + host_DM)

    # Plot the DM curve
    ax1.plot(ds, DM_cumul, 'k')

    # max_y = np.max(DM_cumul)
    if FRB_DM is not None:
        ymax = FRB_DM + 50.
    if ymax is not None:
        max_y = ymax
    # Shade me
    lsz = 14.
    ax1.fill_between((0.1, 10.), 0, max_y, color='green', alpha=0.4)  # ISM
    ax1.text(0.15, max_y * yscl, r'\textbf{Galactic}'+'\n'+r'\textbf{ISM}', color='black', fontsize=lsz, ha='left', va='top')
    ax1.fill_between((10., mnfw_2.r200.value), 0, max_y, color='blue', alpha=0.4)  # Galactic Halo
    ax1.text(12., max_y * yscl, r'\textbf{Galactic}'+'\n'+r'\textbf{Halo}', color='black', fontsize=lsz, ha='left', va='top')
    if M31:
        ax1.fill_between((mnfw_2.r200.value, 2e3), 0, max_y, color='red', alpha=0.4)  # Galactic Halo
        ax1.text(300., max_y * yscl, r'\texgbf{M31}', color='black', fontsize=lsz, ha='left', va='top')

    ax1.set_xscale("log", nonposx='clip')
    # ax.set_yscale("log", nonposy='clip')
    if M31:
        ax1.set_xlim(0.1, 2e3)  # kpc
    else:
        ax1.set_xlim(0.1, mnfw_2.r200.value)  # kpc
    ax1.set_ylim(0., max_y)
    ax1.spines['right'].set_visible(False)
    ax1.set_xlabel(r'\textbf{Distance (kpc)}')
    ax1.set_ylabel(r'\textbf{Cumulative DM (pc cm$^{-3}$)}')

    # IGM
    Gds = np.array(ds) / 1e6
    ax2.plot(Gds, DM_cumul, 'k')
    ax2.spines['left'].set_visible(False)
    ax2.yaxis.tick_right()
    ax2.tick_params(labelright='off')
    ax2.minorticks_on()
    ax2.set_xlabel(r'\textbf{Distance (Gpc)}')

    smax = cosmo.comoving_distance(zFRB).to('Gpc').value
    ax2.fill_between((smin, smax-dsmx), 0, max_y, color='gray', alpha=0.4)  # Cosmic
    ilbl = r'\textbf{Cosmic}'
    ax2.text(0.2, max_y * yscl, ilbl, color='black', fontsize=lsz, ha='left', va='top')

    # Host
    if host_DM > 0.:
        ax2.fill_between((smax-dsmx, smax+dsmx), 0, max_y, color='red', alpha=0.4)  # Host
        ax2.set_xlim(smin, smax+dsmx)  # Gpc
    else:
        ax2.set_xlim(smin, smax)  # Gpc

    if FRB_DM is not None:
        ax1.axhline(y=FRB_DM, ls='--', color='black')
        ax2.axhline(y=FRB_DM, ls='--', color='black')

