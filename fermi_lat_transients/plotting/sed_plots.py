"""
SED (Spectral Energy Distribution) visualization.

Publication-quality plots of individual and stacked SEDs with
upper limits shown as downward arrows/caps.
"""

import numpy as np
import matplotlib.pyplot as plt
from ..constants import TS_DETECTION


# ─── Publication style ──────────────────────────────────────────────────────
PLOT_STYLE = {
    'font.family': 'serif',
    'font.size': 12,
    'axes.linewidth': 1.2,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
}


def plot_sed(sed_points, title=None, ax=None, save_path=None,
             color='royalblue', label=None, ts_threshold=None,
             show_detections=True, show_upper_limits=True):
    """
    Plot an SED (E^2 dN/dE vs E) with detections and/or upper limits.

    Points with TS > ts_threshold are shown as detections (circles),
    others as upper limits (downward arrows).

    Parameters
    ----------
    sed_points : dict
        Output of likelihood.sed.compute_sed_points() or
        stacking.sed_stack.stacked_sed_points().
        Must contain: 'eref', 'e2dnde_mle', 'e2dnde_ul95', 'ts_max'.
    title : str, optional
    ax : matplotlib Axes, optional
    save_path : str, optional
    color : str
        Color for the data points.
    label : str, optional
        Legend label.
    ts_threshold : float, optional
        TS threshold for detection vs UL (default: TS_DETECTION=25).
    show_detections : bool
        Show detection points.
    show_upper_limits : bool
        Show upper limit arrows.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if ts_threshold is None:
        ts_threshold = TS_DETECTION

    with plt.rc_context(PLOT_STYLE):
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))
        else:
            fig = ax.figure

        eref = np.asarray(sed_points['eref'])
        e2_mle = np.asarray(sed_points['e2dnde_mle'])
        e2_ul = np.asarray(sed_points['e2dnde_ul95'])
        ts = np.asarray(sed_points['ts_max'])
        emin = np.asarray(sed_points.get('emin', eref))
        emax = np.asarray(sed_points.get('emax', eref))

        detected = ts >= ts_threshold
        upper_lim = ~detected

        # Energy error bars (bin width)
        xerr_lo = eref - emin
        xerr_hi = emax - eref

        # Detections
        if show_detections and detected.any():
            ax.errorbar(
                eref[detected], e2_mle[detected],
                xerr=[xerr_lo[detected], xerr_hi[detected]],
                fmt='o', color=color, ms=6, lw=1.5, capsize=3,
                label=label or 'Detection',
            )

        # Upper limits
        if show_upper_limits and upper_lim.any():
            ax.errorbar(
                eref[upper_lim], e2_ul[upper_lim],
                xerr=[xerr_lo[upper_lim], xerr_hi[upper_lim]],
                fmt='v', color=color, ms=7, lw=1.5, capsize=4,
                uplims=True,
                label=(label + ' (UL)') if label else '95% CL UL',
            )

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Energy [MeV]')
        ax.set_ylabel('$E^2\\,dN/dE$ [MeV cm$^{-2}$ s$^{-1}$]')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, which='both')

        if title:
            ax.set_title(title)

        if save_path:
            fig.savefig(save_path)

    return fig


def plot_stacked_sed(individual_seds, stacked_sed, title=None,
                     ax=None, save_path=None, source_names=None):
    """
    Plot stacked SED overlaid with individual source SEDs.

    Parameters
    ----------
    individual_seds : list of dict
        List of individual SED point dicts.
    stacked_sed : dict
        Stacked SED points.
    title : str, optional
    ax : matplotlib Axes, optional
    save_path : str, optional
    source_names : list of str, optional
        Names for individual sources.

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(PLOT_STYLE):
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(9, 6))
        else:
            fig = ax.figure

        # Individual sources (thin gray)
        for i, sed in enumerate(individual_seds):
            name = source_names[i] if source_names and i < len(source_names) else None
            eref = np.asarray(sed['eref'])
            e2_ul = np.asarray(sed['e2dnde_ul95'])
            ax.plot(eref, e2_ul, 'v', color='gray', ms=4, alpha=0.4,
                    label=name if i < 3 else None)  # Only label first few

        # Stacked (bold)
        eref = np.asarray(stacked_sed['eref'])
        e2_ul = np.asarray(stacked_sed['e2dnde_ul95'])
        emin = np.asarray(stacked_sed.get('emin', eref))
        emax = np.asarray(stacked_sed.get('emax', eref))
        xerr_lo = eref - emin
        xerr_hi = emax - eref

        ax.errorbar(
            eref, e2_ul,
            xerr=[xerr_lo, xerr_hi],
            fmt='v', color='crimson', ms=9, lw=2, capsize=5,
            uplims=True,
            label='Stacked 95% CL UL',
        )

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Energy [MeV]')
        ax.set_ylabel('$E^2\\,dN/dE$ [MeV cm$^{-2}$ s$^{-1}$]')
        ax.legend(fontsize=9, loc='upper right')
        ax.grid(True, alpha=0.3, which='both')

        if title:
            ax.set_title(title)

        if save_path:
            fig.savefig(save_path)

    return fig
