"""
generalized_halo_plotter.py

Helper module to sweep over k (closure radius in units of r_vir) for the 
generalized g × M_DM halo model and generate diagnostic plots.

The g × M_DM Framework
----------------------
This module visualizes the baryon profile framework where:

    M_B(<r) = g(y) × (Ω_B / (Ω_M - Ω_B)) × M_DM(<r)

where y = c × r / r_vir is the dimensionless radius, and g(y) → 1 as y → ∞,
ensuring the enclosed baryon fraction approaches the cosmic mean at large radii.

The parameter k controls the closure radius r_close = k × r_vir, which is 
the radius at which the enclosed baryon fraction equals the cosmic mean.
Internally, this sets y_out = k × c in the model.

Usage
-----
from generalized_halo_plotter import make_plots
from generalized_halo_model import TanhHaloProfile

kvals = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5]
make_plots(
    k=kvals,
    model_class=TanhHaloProfile,
    model_kwargs=dict(log_Mhalo=12.18, c=7.67),
    save_prefix="g_times_Mdm_profiles",
    show=True
)

References
----------
- Ayromlou et al. 2023 (closure radius concept)
- Bryan & Norman 1998 (virial overdensity)
- Notebook: generalized_halo_profile.ipynb (derivation)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Union, Sequence, Dict, Any

import numpy as np
import matplotlib.pyplot as plt

from chimefrbm31.data.NFW_model.mmNFW.sunil_models import generalised_halo_profile


KCloseLike = Union[float, int, Sequence[Union[float, int]], np.ndarray]


def save_fig(fig, save_file=None, format='pdf', dpi=400, make_transparent=False):
    """
    Save the figure to a file.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Matplotlib figure object
    save_file : str, optional
        Path to save the figure (without extension). If None, will not save.
    format : str
        Format to save the figure, default is 'pdf'
    dpi : int
        Dots per inch for the saved figure
    make_transparent : bool
        If True, also save a transparent PNG version
    """
    if save_file is not None:
        save_file_full = save_file + f".{format}"
        fig.savefig(save_file_full, format=format, dpi=dpi)
        print(f"Figure saved to {save_file_full}")

        if make_transparent:
            save_file_png = save_file + '_transparent.png'
            fig.savefig(save_file_png, format='png', dpi=dpi, transparent=True)
            print(f"Figure saved with transparency to {save_file_png}")
    else:
        print("No file specified, figure not saved.")


@dataclass
class PlotOutputs:
    """
    Container for plot outputs.
    
    Attributes
    ----------
    k_values : np.ndarray
        Array of k values (closure radius in units of r_vir) used
    figures : Dict[str, plt.Figure]
        Dictionary mapping plot names to figure objects
    """
    k_values: np.ndarray
    figures: Dict[str, plt.Figure]


def _as_k_array(k: KCloseLike) -> np.ndarray:
    """
    Normalize k input to a 1D numpy array of floats.
    
    Parameters
    ----------
    k : float, int, or array-like
        Closure radius parameter(s) in units of r_vir
    
    Returns
    -------
    np.ndarray
        1D array of k values
    """
    if isinstance(k, (float, int, np.floating, np.integer)):
        return np.array([float(k)], dtype=float)

    arr = np.array(list(k), dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError("k must be a float/int or a non-empty 1D list/array.")
    return arr


def _set_fontsize(ax, base: float = 14.0) -> None:
    """
    Set consistent font sizes for axis labels and ticks.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object to modify
    base : float
        Base font size for tick labels; titles and axis labels are base + 2
    """
    ax.title.set_fontsize(base + 2)
    ax.xaxis.label.set_fontsize(base + 2)
    ax.yaxis.label.set_fontsize(base + 2)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontsize(base)


def make_plots(
    k: KCloseLike,
    model_class=None,
    add_mNFW_model_class: bool = False,
    model_kwargs: Optional[Dict[str, Any]] = None,
    r_over_rvir: Optional[np.ndarray] = None,
    save_prefix: Optional[str] = None,
    show: bool = True,
) -> PlotOutputs:
    """
    Build one or more generalized halo models and generate diagnostic plots.
    
    This function creates models using the g × M_DM framework for different
    closure radii and produces plots comparing:
    1. Density profiles (ρ_B and ρ_DM vs radius)
    2. Enclosed mass ratios (M_B/M_DM normalized to cosmic expectation)

    Parameters
    ----------
    k : float or array-like
        Closure radius parameter(s) in units of r_vir. This is the radius
        at which the enclosed baryon fraction equals the cosmic mean.
        Internally converted to y_out = k × c for the model.
    model_class : class, optional
        Model class to use. Must accept k=... in __init__ and expose:
        - Attributes: k, c, rvir, cosmo, cosmo_frac
        - Methods: rho_b(xyz), rho_dm(xyz), M_B(y), M_DM(y)
        Defaults to TanhHaloProfile.
    add_mNFW_model_class : bool
        If True, add the standard mNFW profile (from frb.halos.models) 
        for comparison on the density plot.
    model_kwargs : dict, optional
        Extra keyword arguments passed to model_class(...).
        Common options: log_Mhalo, c, z
    r_over_rvir : np.ndarray, optional
        Radius grid in units of r_vir. If None, uses logspace from 
        1e-3 to 20 with 1000 points.
    save_prefix : str, optional
        If provided, saves figures to "<save_prefix>__<name>.pdf".
    show : bool
        If True, call plt.show() at the end.

    Returns
    -------
    PlotOutputs
        Dataclass containing:
        - k_values: array of k values used
        - figures: dict mapping plot names to Figure objects
        
    Examples
    --------
    >>> from generalized_halo_plotter import make_plots
    >>> 
    >>> # Single closure radius
    >>> make_plots(k=1.5, model_kwargs=dict(log_Mhalo=12.18))
    >>> 
    >>> # Multiple closure radii for comparison
    >>> make_plots(
    ...     k=[0.5, 1.0, 1.5, 2.0],
    ...     model_kwargs=dict(log_Mhalo=12.18, c=7.67),
    ...     save_prefix="m31_halo_models",
    ...     show=True
    ... )
    
    Notes
    -----
    The closure radius r_close = k × r_vir is defined as the radius where:
    
        M_B(<r_close) / M_total(<r_close) = Ω_B / Ω_M
    
    For k < 1, baryons are more concentrated than dark matter.
    For k > 1, baryons extend beyond the virial radius before reaching
    the cosmic fraction.
    """
    if model_class is None:
        print("No model_class provided, defaulting to TanhHaloProfile.")
        model_class = generalised_halo_profile.TanhHaloProfile

    model_kwargs = {} if model_kwargs is None else dict(model_kwargs)
    k_values = _as_k_array(k)

    # Build models for each k value
    models = [model_class(k=float(kval), **model_kwargs) for kval in k_values]

    # Optionally load mNFW for comparison
    mnfw = None
    if add_mNFW_model_class:
        try:
            from frb.halos.models import ModifiedNFW
            mnfw = ModifiedNFW(alpha=2, y0=2)
        except ImportError:
            print("Warning: Could not import ModifiedNFW from frb.halos.models")
            add_mNFW_model_class = False

    # Shared radius grid (dimensionless: r / r_vir)
    if r_over_rvir is None:
        r_over_rvir = np.logspace(-3, np.log10(20.0), 1000)
    r_over_rvir = np.asarray(r_over_rvir, dtype=float)

    figures: Dict[str, plt.Figure] = {}

    # Set up colors for different k values
    n_colors = len(k_values)
    if n_colors <= 10:
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i) for i in range(n_colors)]
    else:
        cmap = plt.get_cmap("cividis")
        colors = cmap(np.linspace(0, 1, n_colors))

    # =========================================================================
    # Plot 1: Density Profiles
    # =========================================================================
    # Shows ρ_DM (NFW), ρ_DM × Ω_b/Ω_dm (cosmic scaling), and ρ_B for each k
    
    fig1 = plt.figure(figsize=(7, 7))
    ax = plt.gca()
    
    for i, m in enumerate(models):
        # Convert r/r_vir to physical radius in kpc
        rval_kpc = r_over_rvir * m.rvir.to("kpc").value
        xyz = np.zeros((3, rval_kpc.size))
        xyz[2, :] = rval_kpc
        
        # Get density profiles
        rho_b = m.rho_b(xyz)
        rho_dm = m.rho_dm(xyz)

        # Plot reference profiles only once
        if i == 0:
            ax.plot(r_over_rvir, rho_dm, ls="--", lw=2, 
                    label=r"$\rho_{\rm DM}$ (NFW)", color="black")
            ax.plot(r_over_rvir, rho_dm * (m.cosmo.Ob0 / m.cosmo.Odm0), 
                    ls=":", lw=2, 
                    label=r"$\rho_{\rm DM} \times \Omega_b/\Omega_{\rm dm}$", 
                    color="black")
            
            # Add mNFW comparison if available
            if add_mNFW_model_class and mnfw is not None:
                rho_b_mnfw = mnfw.rho_b(xyz)
                ax.plot(r_over_rvir, rho_b_mnfw, ls="-.", lw=2, 
                        label=r"$\rho_b$ (mNFW)", color="green")
        
        # Plot baryon density for this closure radius
        ax.plot(r_over_rvir, rho_b, lw=2, 
                label=rf"$\rho_b,\ r_{{\rm close}}={m.k:.2f}\,r_{{\rm vir}}$", 
                color=colors[i])

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$r / r_{\rm vir}$")
    ax.set_ylabel(r"Density (g cm$^{-3}$)")
    ax.legend(fontsize=10, loc="best")
    _set_fontsize(ax, 14)
    plt.tight_layout()
    figures["densities"] = fig1

    # =========================================================================
    # Plot 2: Enclosed Mass Ratios
    # =========================================================================
    # Shows (M_B/M_DM) / (Ω_B/(Ω_M - Ω_B)) vs radius
    # This ratio should approach 1 at r = r_close = k × r_vir
    
    fig2 = plt.figure(figsize=(7, 7))
    ax = plt.gca()
    
    # Reference line at cosmic ratio
    ax.axhline(1.0, color="k", ls="--", 
               label=r"Cosmic: $M_b/M_{\rm DM} = \Omega_b/(\Omega_m - \Omega_b)$")

    for i, m in enumerate(models):
        # Convert r/r_vir to dimensionless y = c × r/r_vir
        y_vals = r_over_rvir * m.c
        
        # Get enclosed masses
        M_dm_r = m.M_DM(y_vals)
        M_b_r = m.M_B(y_vals)
        
        # Ratio normalized to cosmic expectation
        # Should equal 1 when M_B/M_DM = Ω_B/(Ω_M - Ω_B)
        ratio = (M_b_r / M_dm_r) / m.cosmo_frac

        ax.plot(r_over_rvir, ratio, lw=2,
                label=rf"$r_{{\rm close}}={m.k:.2f}\,r_{{\rm vir}}$",
                color=colors[i])

        # Vertical line marking closure radius
        ax.axvline(m.k, color=colors[i], ls=":", lw=1, alpha=0.7)

    ax.set_xscale("log")
    ax.set_xlabel(r"$r / r_{\rm vir}$")
    ax.set_ylabel(r"$(M_b/M_{\rm DM}) / (\Omega_b/(\Omega_m - \Omega_b))$")
    ax.set_ylim(0, 1.2)
    ax.legend(fontsize=10, loc="best")
    _set_fontsize(ax, 14)
    plt.tight_layout()
    figures["enclosed_mass_ratios"] = fig2

    # =========================================================================
    # Save figures if requested
    # =========================================================================
    if save_prefix:
        for name, fig in figures.items():
            save_fig(fig=fig, save_file=f"{save_prefix}_{name}")

    if show:
        plt.show()

    return PlotOutputs(
        k_values=k_values,
        figures=figures,
    )


if __name__ == "__main__":
    # Example: sweep over closure radii from 0.5 to 2.5 r_vir
    kvals = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5]
    
    make_plots(
        k=kvals,
        add_mNFW_model_class=True,
        model_kwargs=dict(log_Mhalo=12.18, c=7.67),
        save_prefix="generalised_halo_profiles",
        show=False
    )