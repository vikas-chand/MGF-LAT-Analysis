"""
SED (Spectral Energy Distribution) construction via per-bin profile scans.

For each energy bin, a profile likelihood scan is performed to determine
the flux (or upper limit) in that band, building a broadband SED.
"""

import numpy as np
from .unbinned import freeze_all_params, free_diffuse_norms, fit
from .profile import upper_limit_from_profile
from ..physics.conversions import e2dnde_from_band_flux


def _set_pl2_band(like, src_name, index, emin_bin, emax_bin, flux):
    """
    Configure a PowerLaw2 source for a specific energy band.

    Parameters
    ----------
    like : UnbinnedAnalysis
        Likelihood object.
    src_name : str
        Source name.
    index : float
        Fixed photon index.
    emin_bin, emax_bin : float
        Band energy bounds [MeV].
    flux : float
        Integral flux to set [ph/cm2/s].
    """
    src = like.model[src_name]
    spec = src.funcs['Spectrum']

    # Set energy bounds
    for pname, val in [('LowerLimit', emin_bin), ('UpperLimit', emax_bin)]:
        try:
            par = spec.getParam(pname)
            par.setBounds(max(val * 0.01, 1.0), val * 100.0)
            par.setValue(val)
            idx = like.par_index(src_name, pname)
            like.freeze(idx)
        except Exception:
            pass

    # Set and freeze index
    idx_par = like.par_index(src_name, 'Index')
    par = spec.getParam('Index')
    lo, hi = par.getBounds()
    if index < lo:
        par.setBounds(index * 1.5, hi)
    if index > hi:
        par.setBounds(lo, index * 0.5)
    like[idx_par] = index
    like.freeze(idx_par)

    # Set flux
    flux_par = like.par_index(src_name, 'Integral')
    par = spec.getParam('Integral')
    lo, hi = par.getBounds()
    safe_flux = max(flux, 1e-30)
    if safe_flux < lo:
        par.setBounds(safe_flux * 0.01, hi)
    if safe_flux > hi:
        par.setBounds(lo, safe_flux * 100.0)
    like[flux_par] = safe_flux
    like.freeze(flux_par)


def sed_profile_scan(like, src_name, energy_bins, flux_grid=None,
                     index=-2.0, free_diffuse=True, delta_ts=2.71,
                     optimizer='NEWMINUIT', verbosity=0):
    """
    Perform profile likelihood scans in multiple energy bins to build an SED.

    Parameters
    ----------
    like : UnbinnedAnalysis
        Pre-fit likelihood object with transient source.
    src_name : str
        Transient source name.
    energy_bins : list of float
        Energy bin edges [MeV], e.g., [100, 316, 1000, 3162, 10000].
    flux_grid : array-like, optional
        Flux values to scan per bin. Default: logspace(-12, -4, 100).
    index : float
        Fixed photon index for each bin.
    free_diffuse : bool
        Free diffuse norms during scan.
    delta_ts : float
        Delta-TS for upper limit (2.71 for 95% CL).
    optimizer : str
        Fit optimizer.
    verbosity : int
        Verbosity level.

    Returns
    -------
    dict
        'energy_bins': bin edges,
        'bins': list of dicts per bin, each with:
            'emin', 'emax', 'eref': energy info [MeV],
            'flux_grid', 'logL', 'ts': scan results,
            'flux_mle': best-fit flux [ph/cm2/s],
            'flux_ul95': 95% CL upper limit [ph/cm2/s],
            'ts_max': maximum TS,
            'e2dnde_mle', 'e2dnde_ul95': SED values [MeV/cm2/s].
    """
    if flux_grid is None:
        flux_grid = np.concatenate(([0.0], np.logspace(-12, -4, 100)))
    flux_grid = np.asarray(flux_grid)

    nbins = len(energy_bins) - 1
    results = {
        'energy_bins': energy_bins,
        'bins': [],
    }

    for i in range(nbins):
        emin_bin = energy_bins[i]
        emax_bin = energy_bins[i + 1]
        eref = np.sqrt(emin_bin * emax_bin)

        logL_arr = np.full_like(flux_grid, np.nan, dtype=float)

        for j, f in enumerate(flux_grid):
            try:
                freeze_all_params(like)
                _set_pl2_band(like, src_name, index, emin_bin, emax_bin, max(f, 1e-30))

                if free_diffuse:
                    free_diffuse_norms(like)

                logL_val = fit(like, optimizer=optimizer, verbosity=verbosity)
                logL_arr[j] = -logL_val
            except Exception:
                logL_arr[j] = np.nan

        # Compute TS
        logL_null = logL_arr[0] if flux_grid[0] == 0.0 else np.nanmin(logL_arr)
        ts = 2.0 * (logL_arr - logL_null)
        ts = np.where(np.isnan(ts), 0.0, ts)

        # Best fit and UL
        i_max = np.nanargmax(logL_arr) if np.any(~np.isnan(logL_arr)) else 0
        ts_max = ts[i_max]
        flux_mle = flux_grid[i_max]

        ul_result = upper_limit_from_profile(flux_grid, ts, delta_ts=delta_ts)
        flux_ul95 = ul_result['ul95']

        # Convert to E^2 dN/dE
        _, e2dnde_mle = e2dnde_from_band_flux(flux_mle, emin_bin, emax_bin, index)
        _, e2dnde_ul95 = e2dnde_from_band_flux(flux_ul95, emin_bin, emax_bin, index)

        results['bins'].append({
            'emin': emin_bin,
            'emax': emax_bin,
            'eref': eref,
            'flux_grid': flux_grid,
            'logL': logL_arr,
            'ts': ts,
            'flux_mle': flux_mle,
            'flux_ul95': flux_ul95,
            'ts_max': ts_max,
            'e2dnde_mle': e2dnde_mle,
            'e2dnde_ul95': e2dnde_ul95,
        })

    return results


def compute_sed_points(sed_results):
    """
    Extract a summary SED table from sed_profile_scan results.

    Parameters
    ----------
    sed_results : dict
        Output from sed_profile_scan().

    Returns
    -------
    dict of arrays
        'emin', 'emax', 'eref': energy arrays [MeV],
        'flux_mle', 'flux_ul95': flux arrays [ph/cm2/s],
        'e2dnde_mle', 'e2dnde_ul95': SED arrays [MeV/cm2/s],
        'ts_max': TS per bin.
    """
    bins = sed_results['bins']
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
