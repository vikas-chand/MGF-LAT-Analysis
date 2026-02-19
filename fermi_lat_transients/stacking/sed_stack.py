"""
SED stacking: combine per-energy-bin profile likelihoods across sources.

For each energy bin, the TS profiles from individual sources are summed,
producing a stacked SED with combined constraints per band.
"""

import numpy as np
from ..likelihood.profile import upper_limit_from_profile
from ..physics.conversions import e2dnde_from_band_flux


def stack_sed_profiles(sed_results_list, energy_bins=None):
    """
    Stack SED profile likelihoods from multiple sources.

    For each energy bin, sums the TS profiles across all sources.

    Parameters
    ----------
    sed_results_list : list of dict
        Each element is the output of likelihood.sed.sed_profile_scan().
    energy_bins : list, optional
        Energy bin edges. If None, taken from the first result.

    Returns
    -------
    dict
        'energy_bins': bin edges,
        'bins': list of dicts per bin with stacked profiles and fluxes.
    """
    if len(sed_results_list) == 0:
        raise ValueError("sed_results_list is empty")

    if energy_bins is None:
        energy_bins = sed_results_list[0]['energy_bins']

    nbins = len(energy_bins) - 1
    stacked_bins = []

    for i in range(nbins):
        emin = energy_bins[i]
        emax = energy_bins[i + 1]
        eref = np.sqrt(emin * emax)

        # Collect TS profiles for this bin from all sources
        flux_grid = None
        ts_list = []

        for sed_result in sed_results_list:
            if i < len(sed_result['bins']):
                b = sed_result['bins'][i]
                if flux_grid is None:
                    flux_grid = np.asarray(b['flux_grid'])
                ts_i = np.asarray(b['ts'])

                if len(ts_i) == len(flux_grid):
                    ts_list.append(ts_i)
                else:
                    ts_interp = np.interp(
                        flux_grid, b['flux_grid'], ts_i,
                        left=0.0, right=0.0
                    )
                    ts_list.append(ts_interp)

        if flux_grid is None or len(ts_list) == 0:
            stacked_bins.append({
                'emin': emin, 'emax': emax, 'eref': eref,
                'ts_stacked': np.array([]),
                'flux_grid': np.array([]),
                'flux_mle': 0.0, 'flux_ul95': 0.0,
                'ts_max': 0.0,
                'e2dnde_mle': 0.0, 'e2dnde_ul95': 0.0,
            })
            continue

        ts_stacked = np.sum(ts_list, axis=0)

        # Best fit and UL from stacked profile
        i_max = np.argmax(ts_stacked)
        ts_max = ts_stacked[i_max]
        flux_mle = flux_grid[i_max]

        ul_result = upper_limit_from_profile(flux_grid, ts_stacked)
        flux_ul95 = ul_result['ul95']

        # Convert to E^2 dN/dE
        _, e2dnde_mle = e2dnde_from_band_flux(flux_mle, emin, emax, -2.0)
        _, e2dnde_ul95 = e2dnde_from_band_flux(flux_ul95, emin, emax, -2.0)

        stacked_bins.append({
            'emin': emin,
            'emax': emax,
            'eref': eref,
            'flux_grid': flux_grid,
            'ts_stacked': ts_stacked,
            'flux_mle': flux_mle,
            'flux_ul95': flux_ul95,
            'ts_max': ts_max,
            'e2dnde_mle': e2dnde_mle,
            'e2dnde_ul95': e2dnde_ul95,
            'n_sources': len(ts_list),
        })

    return {
        'energy_bins': energy_bins,
        'bins': stacked_bins,
    }


def stacked_sed_points(stacked_sed):
    """
    Extract a summary table from stacked SED results.

    Parameters
    ----------
    stacked_sed : dict
        Output of stack_sed_profiles().

    Returns
    -------
    dict of arrays
        'emin', 'emax', 'eref', 'flux_mle', 'flux_ul95',
        'e2dnde_mle', 'e2dnde_ul95', 'ts_max'.
    """
    bins = stacked_sed['bins']
    return {
        'emin': np.array([b['emin'] for b in bins]),
        'emax': np.array([b['emax'] for b in bins]),
        'eref': np.array([b['eref'] for b in bins]),
        'flux_mle': np.array([b['flux_mle'] for b in bins]),
        'flux_ul95': np.array([b['flux_ul95'] for b in bins]),
        'e2dnde_mle': np.array([b['e2dnde_mle'] for b in bins]),
        'e2dnde_ul95': np.array([b['e2dnde_ul95'] for b in bins]),
        'ts_max': np.array([b['ts_max'] for b in bins]),
    }
