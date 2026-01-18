"""
mmnfw_plotter.py

Helper module to sweep over k for the m^2NFW "closure" model
and generate the standard diagnostic plots.

Usage
-----
from mmnfw_plotter import make_plots
from old_m31_models import mClosureModifiedNFW   # or your module path

kvals = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5]
make_plots(
    k=kvals,
    model_class=mClosureModifiedNFW,
    model_kwargs=dict(y0=2.0,log_Mhalo=12.2, c=7.67),
    save_prefix="m2NFW_y0_2p0",
    show=True
)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Union, Sequence, Dict, Any

import numpy as np
import matplotlib.pyplot as plt

from chimefrbm31.data.NFW_model.mmNFW.sunil_models import generalised_halo_profile


KCloseLike = Union[float, int, Sequence[Union[float, int]], np.ndarray]



def save_fig(fig,save_file=None,format='pdf',dpi=400,make_transparent=False):
    """
    Save the figure to a file.
    :param fig: Matplotlib figure object
    :param save_file: Path to save the figure, if None, will not save
    :param format: Format to save the figure, default is 'pdf'
    :param dpi: Dots per inch for the saved figure
    """
    if save_file is not None:
        save_file = save_file + f".{format}"
        fig.savefig(save_file, format=format, dpi=dpi)
        print(f"Figure saved to {save_file}")

        if make_transparent:
            save_file = save_file.replace('.' + format, '_transparent.' + 'png')
            fig.savefig(save_file, format='png', dpi=dpi, transparent=True)
            print(f"Figure saved with transparency to {save_file}")
    else:
        print("No file specified, figure not saved.")


@dataclass
class PlotOutputs:
    """Returned handles and computed arrays."""
    k_values: np.ndarray
    figures: Dict[str, plt.Figure]


def _as_k_array(k: KCloseLike) -> np.ndarray:
    """Normalize k input to a 1D numpy array of floats."""
    if isinstance(k, (float, int, np.floating, np.integer)):
        return np.array([float(k)], dtype=float)

    arr = np.array(list(k), dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError("k must be a float/int or a non-empty 1D list/array.")
    return arr


def _set_fontsize(ax, base: float = 14.0) -> None:
    """Simple fontsize helper (replacement for fig_utils.set_fontsize)."""
    ax.title.set_fontsize(base + 2)
    ax.xaxis.label.set_fontsize(base + 2)
    ax.yaxis.label.set_fontsize(base + 2)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontsize(base)


def make_plots(
    k: KCloseLike,
    model_class = None,
    add_mNFW_model_class: bool = False,
    model_kwargs: Optional[Dict[str, Any]] = None,
    r_over_rvir: Optional[np.ndarray] = None,
    save_prefix: Optional[str] = None,
    show: bool = True,
) -> PlotOutputs:
    """
    Build one or more m^2NFW closure models and generate diagnostic plots.

    Parameters
    ----------
    k:
        float or list/array of floats representing the k parameter.
    model_class:
        Your model class (e.g., mClosureModifiedNFW). Must accept k=... in __init__.
        The object must expose: k, y0, c, rvir, rho0, rho0_b, rho_b_A,
        plus methods: rho_b(xyz), rho_dm(xyz), fy_b(y), fy_dm(y).
    model_kwargs:
        Extra keyword arguments passed into model_class(...).
    r_over_rvir:
        Optional radius grid in units of r_vir. If None, uses
        logspace from 1e-3 to 20 with 1000 points.
    save_prefix:
        If provided, saves figures to "<save_prefix>__<name>.png" at 300 dpi.
    show:
        If True, plt.show() figures at the end.

    Returns
    -------
    PlotOutputs with computed arrays and figure handles.
    """
    if model_class is None:
        print("No model_class provided, defaulting to mClosureModifiedNFW.")
        model_class=generalised_halo_profile.TanhHaloProfile

    model_kwargs = {} if model_kwargs is None else dict(model_kwargs)
    k_values = _as_k_array(k)

    # Build models
    models = [model_class(k=float(k), **model_kwargs) for k in k_values]

    if add_mNFW_model_class: 
        from frb.halos.models import ModifiedNFW
        mnfw = ModifiedNFW(alpha=2, y0=2, )

    # Extract sweep outputs

    # Shared radius grid (dimensionless)
    if r_over_rvir is None:
        r_over_rvir = np.logspace(-3, np.log10(20.0), 1000)
    r_over_rvir = np.asarray(r_over_rvir, dtype=float)

    figures: Dict[str, plt.Figure] = {}

    # Let's set up colors
    n_colors = len(k_values)
    if n_colors <= 10:
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i) for i in range(n_colors)]
    else:
        cmap = plt.get_cmap("cividis")
        colors = cmap(np.linspace(0, 1, n_colors))


    # ---------- Plot 1: densities ----------
    fig1 = plt.figure(figsize=(7, 7))
    ax = plt.gca()
    for i, m in enumerate(models):
        rval_kpc = r_over_rvir * m.rvir.to("kpc").value
        xyz = np.zeros((3, rval_kpc.size))
        xyz[2, :] = rval_kpc
        rho_b = m.rho_b(xyz)
        rho_dm = m.rho_dm(xyz)

        if i == 0:
            ax.plot(r_over_rvir, rho_dm, ls="--", lw=2, label=r"$\rho_{dm}$ (NFW)", color="black")
            ax.plot(r_over_rvir, rho_dm*(m.cosmo.Ob0/m.cosmo.Odm0), ls=":", lw=2, label=r"$\rho_{dm}\,\Omega_b/\Omega_{dm}$", color="black")
            if add_mNFW_model_class:
                rho_b_mnfw = mnfw.rho_b(xyz)
                ax.plot(r_over_rvir, rho_b_mnfw, ls="-.", lw=2, label=r"$\rho_b$ (mNFW)")
        ax.plot(r_over_rvir, rho_b, lw=2, label=rf"$\rho_b,\ r_{{close}}={m.k:.2f} r_{{vir}}$", color=colors[i])

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$r/r_{\rm vir}$")
    ax.set_ylabel(r"Density (g cm$^{-3}$)")
    ax.legend(fontsize=10, loc="best")
    _set_fontsize(ax, 14)
    plt.tight_layout()
    figures["densities"] = fig1

    # ---------- Plot 2: enclosed mass ratios ----------
    omega_ref = None
    try:
        cosmo = models[0].cosmo
        omega_ref = (cosmo.Ob0 / cosmo.Odm0).value if hasattr(cosmo.Ob0, "value") else (cosmo.Ob0 / cosmo.Odm0)
    except Exception:
        omega_ref = None

    fig2 = plt.figure(figsize=(7, 7))
    ax = plt.gca()
    ax.axhline(1.0, color="k", ls="--", label=r"cosmic: $(M_b/M_{dm})=\Omega_B/(\Omega_{M} - \Omega_{B})$")

    for i, m in enumerate(models):
        M_dm_r = m.M_DM(r_over_rvir * m.c)
        M_b_r  = m.M_B(r_over_rvir * m.c)
        ratio = (M_b_r / M_dm_r)/m.cosmo_frac

        ax.plot(r_over_rvir, ratio, lw=2,
                label=rf"$r_{{close}}={m.k:.2f} r_{{vir}}$",
                color=colors[i]
        )

        # Adding thin vertical line at k
        ax.axvline(m.k, color=colors[i], ls=":", lw=1, alpha=0.7)

    ax.set_xscale("log")
    ax.set_xlabel(r"$r / r_{\rm vir}$", fontsize=14)
    ax.set_ylabel(r"$(M_b/M_{dm})/(\Omega_b/\Omega_{dm})$" if omega_ref is not None else r"$M_b/M_{dm}$")
    ax.legend(fontsize=10, loc="best")
    _set_fontsize(ax, 14)
    plt.tight_layout()
    figures["enclosed_mass_ratios"] = fig2

    # Save if requested
    if save_prefix:
        for name, fig in figures.items():
            save_fig(fig=fig,
                    save_file=f"{save_prefix}__{name}",
                                )

    if show:
        plt.show()

    return PlotOutputs(
        k_values=k_values,
        figures=figures,
    )


if __name__ == "__main__":
    kvals = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5]
    make_plots(
    k=kvals,
    add_mNFW_model_class=True,
    model_kwargs=dict(log_Mhalo=12.18, c=7.67),
    save_prefix="generalised_halo_profiles",
    show=False
)