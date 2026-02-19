"""
Profile likelihood stacking.

Combines profile likelihood curves from multiple sources by summing
their TS (test statistic) profiles on a common flux grid, producing
a combined constraint that is more sensitive than any individual source.
"""

import numpy as np
from ..likelihood.profile import upper_limit_from_profile


def stack_profiles(profile_list, flux_grid=None):
    """
    Stack profile likelihood curves from multiple sources.

    The stacked TS is the sum of individual TS profiles evaluated
    on a common flux grid. Each profile is interpolated onto the
    common grid if needed.

    Parameters
    ----------
    profile_list : list of dict
        Each dict must contain 'flux_grid' and 'ts' arrays
        (output of likelihood.profile.profile_scan).
    flux_grid : array-like, optional
        Common flux grid for stacking. If None, uses the grid from
        the first profile.

    Returns
    -------
    dict
        'flux_grid': common flux array,
        'ts_stacked': summed TS array,
        'ts_individual': list of individual TS arrays (on common grid),
        'n_sources': number of sources stacked.
    """
    if len(profile_list) == 0:
        raise ValueError("profile_list is empty")

    if flux_grid is None:
        flux_grid = np.asarray(profile_list[0]['flux_grid'])
    else:
        flux_grid = np.asarray(flux_grid)

    ts_individual = []
    for prof in profile_list:
        fg = np.asarray(prof['flux_grid'])
        ts = np.asarray(prof['ts'])

        if np.array_equal(fg, flux_grid):
            ts_interp = ts
        else:
            ts_interp = np.interp(flux_grid, fg, ts, left=0.0, right=0.0)

        ts_individual.append(ts_interp)

    ts_stacked = np.sum(ts_individual, axis=0)

    return {
        'flux_grid': flux_grid,
        'ts_stacked': ts_stacked,
        'ts_individual': ts_individual,
        'n_sources': len(profile_list),
    }


def stacked_upper_limit(stacked_result, delta_ts=2.71):
    """
    Extract the upper limit from a stacked profile.

    Parameters
    ----------
    stacked_result : dict
        Output of stack_profiles().
    delta_ts : float
        Delta-TS threshold (2.71 for 95% CL).

    Returns
    -------
    dict
        'ul95': stacked upper limit [ph/cm2/s],
        'best_fit': best-fit flux,
        'ts_max': maximum stacked TS,
        'n_sources': number of sources stacked.
    """
    ul = upper_limit_from_profile(
        stacked_result['flux_grid'],
        stacked_result['ts_stacked'],
        delta_ts=delta_ts,
    )
    ul['n_sources'] = stacked_result['n_sources']
    return ul


def evolution_plot_data(profile_list, flux_grid=None):
    """
    Compute stacking evolution: how the stacked TS grows as sources
    are added one by one.

    Parameters
    ----------
    profile_list : list of dict
        Individual profile scan results.
    flux_grid : array-like, optional
        Common flux grid.

    Returns
    -------
    dict
        'n_sources': array [1, 2, ..., N],
        'ts_max_cumulative': max stacked TS after adding each source,
        'ul95_cumulative': stacked UL after adding each source.
    """
    if flux_grid is None:
        flux_grid = np.asarray(profile_list[0]['flux_grid'])

    n_arr = []
    ts_max_arr = []
    ul_arr = []

    for n in range(1, len(profile_list) + 1):
        result = stack_profiles(profile_list[:n], flux_grid=flux_grid)
        ul = stacked_upper_limit(result)
        n_arr.append(n)
        ts_max_arr.append(ul['ts_max'])
        ul_arr.append(ul['ul95'])

    return {
        'n_sources': np.array(n_arr),
        'ts_max_cumulative': np.array(ts_max_arr),
        'ul95_cumulative': np.array(ul_arr),
    }
