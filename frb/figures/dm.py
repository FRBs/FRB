""" Module for basic plots related to FRB host and foreground galaxies"""
import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS


from pkg_resources import resource_filename
from matplotlib import pyplot as plt

from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import z_at_value
from astropy import units

from frb.halos.models import ModifiedNFW, M31
from frb.halos import utils as frb_halos
from frb.dm import igm as frb_igm
from frb.dm import cosmic
from frb.figures import utils as ff_utils
from frb import defs

from ne2001 import density

from IPython import embed

def sub_cartoon(ax1, ax2, coord, zFRB, halos=False, host_DM=50., ymax=None,
                IGM_only=True, smin=0.1, cosmo=defs.frb_cosmo,
                show_M31=None, fg_halos=None, dsmx=0.05, FRB_DM=None, yscl = 0.97):
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
        show_M31 (bool, optional):
            Include M31 in the calculation?
            NOT IMPLEMENTED RIGHT NOW
        fg_halos (dict or Table):
            Used to add to DM_IGM
            Keys must include 'z' 'DM' 'lbl'
        smin (float, optional):
            Minimum value in axis2 (Gpc)
        dsmx (float, optional): Padding on the x-axis;  Gpc
            Allows for host.  Set to 0 to ignore host
        FRB_DM (float): Observed value;  sets ymax = FRB_DM+50
        yscl (float, optional):
            Controls placement of labels
        cosmo (astropy.cosmology, optional):
            Defaults to Repo Cosmology

    """
    if halos:
        # NOT READY FOR THIS
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

    if show_M31:
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
    if show_M31:
        ax1.fill_between((mnfw_2.r200.value, 2e3), 0, max_y, color='red', alpha=0.4)  # Galactic Halo
        ax1.text(300., max_y * yscl, r'\texgbf{M31}', color='black', fontsize=lsz, ha='left', va='top')

    ax1.set_xscale("log", nonposx='clip')
    # ax.set_yscale("log", nonposy='clip')
    if show_M31:
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
    #ax2.fill_between((0.1, smax-dsmx), 0, max_y, color='gray', alpha=0.4)  # Galactic Halo
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
        ax1.axhline(y=FRB_DM, ls='--', color='black', lw=3)
        ax2.axhline(y=FRB_DM, ls='--', color='black', lw=3)


def fig_cosmic(frbs, clrs=None, outfile=None, multi_model=False, no_curves=False, 
               F = 0.2,  cosmo=defs.frb_cosmo,
               widen=False, show_nuisance=False, ax=None, zmax=0.75,
               show_sigmaDM=False, cl=(16,84), beta=3., gold_only=True, gold_frbs=None):
    """

    Args:
        frbs (list):
            list of FRB objects
        clrs (list, optional):
        outfile (str, optional):
        multi_model (deprecated):
        no_curves (bool, optional):
            If True, just show the data
        F (float, optional):
            Feedback parameter
        widen (bool, optional):
            If True, make the plot wide
        show_nuisance (bool, optional):
            if True, add a label giving the Nuiscance value
        show_sigmaDM (bool, optional):
            If True, show a model estimate of the scatter in the DM relation
        cl (tuple, optional):
            Confidence limits for the scatter
        beta (float, optional):
            Parameter to the DM scatter estimation
        gold_only (bool, optional):
            If True, limit to the gold standard sample
        gold_frbs (list, optional):
            List of gold standard FRBs
        ax (matplotlib.Axis, optional):
            Use this axis instead of creating one
        zmax (float, optional):
            Max redshift for the MR line
        cosmo (astropy.cosmology, optional):
            Defaults to Repo Cosmology

    Returns:
        dict:
        ax (optional): if outfile is None

    """
    # Init
    if gold_frbs is None:
        gold_frbs = cosmic.gold_frbs

    # Plotting
    ff_utils.set_mplrc()

    bias_clr = 'darkgray'

    return_stuff = {}

    # Start the plot
    if ax is None:
        if widen:
            fig = plt.figure(figsize=(12, 8))
        else:
            fig = plt.figure(figsize=(8, 8))
        plt.clf()
        ax = plt.gca()

    # DM_cosmic from cosmology
    DM_cosmic, zeval = frb_igm.average_DM(zmax, cumul=True)
    return_stuff['zeval'] = zeval
    return_stuff['DM_cosmic'] = DM_cosmic
    DMc_spl = IUS(zeval, DM_cosmic)
    if not no_curves:
        ax.plot(zeval, DM_cosmic, 'k--', label=cosmo.name, lw=2)

    if multi_model:
        # Change Omega_b
        cosmo_highOb = FlatLambdaCDM(Ob0=cosmo.Ob0*1.2, Om0=cosmo.Om0, H0=cosmo.H0)
        DM_cosmic_high, zeval_high = frb_igm.average_DM(zmax, cumul=True, cosmo=cosmo_highOb)
        ax.plot(zeval_high, DM_cosmic_high, '--', color='gray', label=r'DM$_{\rm cosmic} (z) \;\; [1.2 \times \Omega_b]$')
        # Change H0
        cosmo_lowH0 = FlatLambdaCDM(Ob0=cosmo.Ob0, Om0=cosmo.Om0, H0=cosmo.H0/1.2)
        DM_cosmic_lowH0, zeval_lowH0 = frb_igm.average_DM(zmax, cumul=True, cosmo=cosmo_lowH0)
        ax.plot(zeval_lowH0, DM_cosmic_lowH0, ':', color='gray', label=r'DM$_{\rm cosmic} (z) \;\; [H_0/1.2]$')

    if show_sigmaDM:
        #f_C0 = frb_cosmology.build_C0_spline()
        f_C0_3 = cosmic.grab_C0_spline(beta=beta)
        # Updated
        nstep=50
        sigma_DM = F * zeval**(-0.5) #* DM_cosmic.value
        sub_sigma_DM = sigma_DM[::nstep]
        sub_z = zeval[::nstep]
        sub_DM = DM_cosmic.value[::nstep]
        # Loop me
        sigmas, C0s, sigma_lo, sigma_hi = [], [], [], []
        for kk, isigma in enumerate(sub_sigma_DM):
            #res = frb_cosmology.minimize_scalar(frb_cosmology.deviate2, args=(f_C0, isigma))
            #sigmas.append(res.x)
            sigmas.append(isigma)
            C0s.append(float(f_C0_3(isigma)))
            # PDF
            PDF = cosmic.DMcosmic_PDF(cosmic.Delta_values, C0s[-1], sigma=sigmas[-1], 
                                      beta=beta)
            cumsum = np.cumsum(PDF) / np.sum(PDF)
            #if sub_DM[kk] > 200.:
            #    embed(header='131')
            # DO it
            DM = cosmic.Delta_values * sub_DM[kk]
            sigma_lo.append(DM[np.argmin(np.abs(cumsum-cl[0]/100))])
            sigma_hi.append(DM[np.argmin(np.abs(cumsum-cl[1]/100))])
        # Plot
        ax.fill_between(sub_z, sigma_lo, sigma_hi, # DM_cosmic.value-sigma_DM, DM_cosmic.value+sigma_DM,
                        color='gray', alpha=0.3)

    # Do each FRB
    DM_subs = []
    for ifrb in frbs:
        DM_sub = ifrb.DM - ifrb.DMISM
        DM_subs.append(DM_sub.value)
    DM_subs = np.array(DM_subs)

    # chi2
    DMs_MW_host = np.linspace(30., 100., 100)
    zs = np.array([ifrb.z for ifrb in frbs])
    DM_theory = DMc_spl(zs)

    chi2 = np.zeros_like(DMs_MW_host)
    for kk,DM_MW_host in enumerate(DMs_MW_host):
        chi2[kk] = np.sum(((DM_subs-DM_MW_host)-DM_theory)**2)

    imin = np.argmin(chi2)
    DM_MW_host_chisq = DMs_MW_host[imin]
    print("DM_nuisance = {}".format(DM_MW_host))

    # MW + Host term
    def DM_MW_host(z, min_chisq=False):
        if min_chisq:
            return DM_MW_host_chisq
        else:
            return 50. + 50./(1+z)

    # Gold FRBs
    for kk,ifrb in enumerate(frbs):
        if ifrb.frb_name not in gold_frbs:
            continue
        if clrs is not None:
            clr = clrs[kk]
        else:
            clr = None
        ax.scatter([ifrb.z], [DM_subs[kk]-DM_MW_host(ifrb.z)],
                        label=ifrb.frb_name, marker='s', s=90, color=clr)

    # ################################
    # Other FRBs
    s_other = 90

    if not gold_only:
        labeled = False
        for kk, ifrb in enumerate(frbs):
            if ifrb.frb_name in gold_frbs:
                continue
            if not labeled:
                lbl = "Others"
                labeled = True
            else:
                lbl = None
            ax.scatter([ifrb.z], [ifrb.DM.value -
                                      ifrb.DMISM.value - DM_MW_host(ifrb.z)],
                   label=lbl, marker='o', s=s_other, color=bias_clr)


    legend = ax.legend(loc='upper left', scatterpoints=1, borderpad=0.2,
                        handletextpad=0.3, fontsize=19)
    ax.set_xlim(0, 0.7)
    ax.set_ylim(0, 1000.)
    #ax.set_xlabel(r'$z_{\rm FRB}$', fontname='DejaVu Sans')
    ax.set_xlabel(r'$z_{\rm FRB}$', fontname='DejaVu Sans')
    ax.set_ylabel(r'$\rm DM_{cosmic} \; (pc \, cm^{-3})$', fontname='DejaVu Sans')

    #
    if show_nuisance:
        ax.text(0.05, 0.60, r'$\rm DM_{MW,halo} + DM_{host} = $'+' {:02d} pc '.format(int(DM_MW_host))+r'cm$^{-3}$',
            transform=ax.transAxes, fontsize=23, ha='left', color='black')

    ff_utils.set_fontsize(ax, 23.)

    # Layout and save
    if outfile is not None:
        plt.tight_layout(pad=0.2,h_pad=0.1,w_pad=0.1)
        plt.savefig(outfile, dpi=400)
        print('Wrote {:s}'.format(outfile))
        plt.close()
        return return_stuff
    else:
        return ax, return_stuff

