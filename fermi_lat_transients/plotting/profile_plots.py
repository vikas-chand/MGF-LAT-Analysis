"""
Profile likelihood and upper limit visualization.
"""

import numpy as np
import matplotlib.pyplot as plt


# ─── Publication style ──────────────────────────────────────────────────────
PLOT_STYLE = {
    'font.family': 'serif',
    'font.size': 12,
    'axes.linewidth': 1.2,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
}


def plot_profile(flux_grid, test_statistic, ul_result=None, title=None,
                 ax=None, save_path=None, **kwargs):
    """
    Plot a profile likelihood (TS vs flux) curve with optional UL marker.

    Parameters
    ----------
    flux_grid : array
        Flux values [ph/cm2/s].
    test_statistic : array
        TS values.
    ul_result : dict, optional
        Output from upper_limit_from_profile() with 'ul95', 'ts_max', etc.
    title : str, optional
        Plot title.
    ax : matplotlib Axes, optional
        Axes to plot on. If None, creates new figure.
    save_path : str, optional
        Path to save the figure.

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(PLOT_STYLE):
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(7, 4.5))
        else:
            fig = ax.figure

        # Main TS curve
        ax.plot(flux_grid, test_statistic, 'b-', lw=2, label='TS profile')

        # Mark UL
        if ul_result is not None:
            ul = ul_result['ul95']
            ts_max = ul_result['ts_max']
            threshold = ul_result.get('threshold', ts_max - 2.71)

            # Horizontal threshold line
            ax.axhline(threshold, color='gray', ls='--', lw=1,
                       label=f'$\\Delta$TS = 2.71')

            # Vertical UL line
            ax.axvline(ul, color='r', ls='--', lw=1.5,
                       label=f'UL$_{{95\\%}}$ = {ul:.2e}')

            # Best fit marker
            ax.axvline(ul_result['best_fit'], color='green', ls=':',
                       lw=1, alpha=0.7, label=f'Best fit = {ul_result["best_fit"]:.2e}')

        ax.set_xlabel('Flux [ph cm$^{-2}$ s$^{-1}$]')
        ax.set_ylabel('TS')
        ax.set_xscale('log')

        if title:
            ax.set_title(title)
        ax.legend(fontsize=9, loc='upper right')
        ax.grid(True, alpha=0.3)

        if save_path:
            fig.savefig(save_path)

    return fig


def plot_stacked_profile(stacked_result, ul_result=None,
                         show_individual=True, title=None,
                         ax=None, save_path=None):
    """
    Plot stacked profile likelihood curve with optional individual curves.

    Parameters
    ----------
    stacked_result : dict
        Output of stacking.stack.stack_profiles().
    ul_result : dict, optional
        Output of stacking.stack.stacked_upper_limit().
    show_individual : bool
        If True, show individual TS curves in gray.
    title : str, optional
        Plot title.
    ax : matplotlib Axes, optional
    save_path : str, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    with plt.rc_context(PLOT_STYLE):
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(8, 5))
        else:
            fig = ax.figure

        flux_grid = stacked_result['flux_grid']

        # Individual curves
        if show_individual and 'ts_individual' in stacked_result:
            for i, ts_i in enumerate(stacked_result['ts_individual']):
                ax.plot(flux_grid, ts_i, color='gray', lw=0.8, alpha=0.5,
                        label='Individual' if i == 0 else None)

        # Stacked curve
        ax.plot(flux_grid, stacked_result['ts_stacked'], 'b-', lw=2.5,
                label=f'Stacked (N={stacked_result["n_sources"]})')

        # UL marker
        if ul_result is not None:
            ax.axvline(ul_result['ul95'], color='r', ls='--', lw=1.5,
                       label=f'Stacked UL$_{{95\\%}}$ = {ul_result["ul95"]:.2e}')
            ax.axhline(ul_result['threshold'], color='gray', ls='--',
                       lw=1, alpha=0.6)

        ax.set_xlabel('Flux [ph cm$^{-2}$ s$^{-1}$]')
        ax.set_ylabel('TS')
        ax.set_xscale('log')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        if title:
            ax.set_title(title)

        if save_path:
            fig.savefig(save_path)

    return fig
