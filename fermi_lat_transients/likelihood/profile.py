"""
Profile likelihood scanning and upper limit extraction.

Provides functions to scan the likelihood as a function of source flux
and extract confidence intervals and upper limits.
"""

import numpy as np
from .unbinned import freeze_all_params, free_diffuse_norms, set_flux, fit


def profile_scan(like, src_name, flux_grid=None, free_diffuse=True,
                 optimizer='NEWMINUIT', verbosity=0):
    """
    Perform a profile likelihood scan over source flux.

    For each flux value in the grid:
    1. Set source flux to the grid value and freeze it
    2. Optionally free diffuse normalizations
    3. Re-fit the model
    4. Record the log-likelihood

    Parameters
    ----------
    like : UnbinnedAnalysis
        Configured likelihood object (should be pre-fit).
    src_name : str
        Name of the transient source.
    flux_grid : array-like, optional
        Flux values to scan [ph/cm2/s].
        Default: [0] + logspace(-11, -3, 150).
    free_diffuse : bool
        If True, free diffuse norms at each scan point.
    optimizer : str
        Optimizer for re-fitting.
    verbosity : int
        Verbosity level.

    Returns
    -------
    dict
        'flux_grid': array of flux values,
        'logL': array of log-likelihood at each flux,
        'logL_null': log-likelihood at flux=0 (null hypothesis),
        'ts': array of TS = 2*(logL - logL_null),
        'ts_max': maximum TS value,
        'best_fit_flux': flux at maximum TS.
    """
    if flux_grid is None:
        flux_grid = np.concatenate(([0.0], np.logspace(-11, -3, 150)))
    flux_grid = np.asarray(flux_grid)

    logL = np.full_like(flux_grid, np.nan, dtype=float)

    for i, f in enumerate(flux_grid):
        try:
            # Freeze everything
            freeze_all_params(like)

            # Set source flux
            set_flux(like, src_name, max(f, 1e-30))

            # Free diffuse norms
            if free_diffuse:
                free_diffuse_norms(like)

            # Fit
            logL_val = fit(like, optimizer=optimizer, verbosity=verbosity)
            logL[i] = -logL_val  # fit() returns -logL, we want logL
        except Exception as e:
            logL[i] = np.nan

    # Null hypothesis (flux = 0 or first point)
    logL_null = logL[0] if flux_grid[0] == 0.0 else np.nanmin(logL)

    # Test statistic
    ts = 2.0 * (logL - logL_null)
    ts = np.where(np.isnan(ts), 0.0, ts)

    # Best fit
    valid = ~np.isnan(logL)
    if valid.any():
        i_max = np.nanargmax(logL)
        ts_max = ts[i_max]
        best_fit_flux = flux_grid[i_max]
    else:
        ts_max = 0.0
        best_fit_flux = 0.0

    return {
        'flux_grid': flux_grid,
        'logL': logL,
        'logL_null': logL_null,
        'ts': ts,
        'ts_max': ts_max,
        'best_fit_flux': best_fit_flux,
    }


def upper_limit_from_profile(flux_grid, test_statistic, delta_ts=2.71,
                             method='right_side'):
    """
    Extract upper limit from a profile likelihood (TS) curve.

    The 95% CL upper limit is the flux value where the TS drops by
    delta_ts from its maximum (one-sided).

    Parameters
    ----------
    flux_grid : array
        Flux values [ph/cm2/s].
    test_statistic : array
        TS values at each flux point.
    delta_ts : float
        Delta-TS threshold (2.71 for 95% CL, 6.63 for 99% CL).
    method : str
        'right_side' - find where TS drops below (ts_max - delta_ts) on
        the right side of the peak.

    Returns
    -------
    dict
        'ul95': upper limit flux [ph/cm2/s],
        'best_fit': best-fit flux at TS peak,
        'ts_max': maximum TS value,
        'threshold': TS threshold used.
    """
    flux_grid = np.asarray(flux_grid)
    ts = np.asarray(test_statistic)

    # Find peak
    i_max = np.argmax(ts)
    ts_max = ts[i_max]
    best_fit = flux_grid[i_max]

    # Threshold
    ts_threshold = ts_max - delta_ts

    if method == 'right_side':
        # Search right of peak for crossing
        right_flux = flux_grid[i_max:]
        right_ts = ts[i_max:]

        # Find first crossing below threshold
        below = right_ts < ts_threshold
        if below.any():
            i_cross = np.argmax(below)
            # Linear interpolation for precision
            if i_cross > 0:
                f1, f2 = right_flux[i_cross - 1], right_flux[i_cross]
                t1, t2 = right_ts[i_cross - 1], right_ts[i_cross]
                if abs(t2 - t1) > 1e-30:
                    ul = f1 + (f2 - f1) * (ts_threshold - t1) / (t2 - t1)
                else:
                    ul = f2
            else:
                ul = right_flux[0]
        else:
            # TS never drops below threshold - UL is at edge of grid
            ul = flux_grid[-1]
    else:
        raise ValueError(f"Unknown method: {method}")

    return {
        'ul95': ul,
        'best_fit': best_fit,
        'ts_max': ts_max,
        'threshold': ts_threshold,
    }
