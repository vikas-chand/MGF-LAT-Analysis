"""
File I/O utilities for Fermi-LAT transient analysis.

Functions for discovering data files, validating GTIs,
and saving/loading profile likelihood results.
"""

import numpy as np
from pathlib import Path
from astropy.io import fits


def find_ft1_ft2(data_dir, ft1_pattern="*_EV*.fits", ft2_pattern="*_SC*.fits"):
    """
    Discover FT1 (event) and FT2 (spacecraft) FITS files in a directory.

    Parameters
    ----------
    data_dir : str or Path
        Directory containing Fermi data files.
    ft1_pattern : str
        Glob pattern for event files.
    ft2_pattern : str
        Glob pattern for spacecraft files.

    Returns
    -------
    dict
        {'ft1': Path or None, 'ft2': Path or None}
    """
    data_dir = Path(data_dir)
    ft1_files = sorted(data_dir.glob(ft1_pattern))
    ft2_files = sorted(data_dir.glob(ft2_pattern))

    return {
        'ft1': ft1_files[0] if ft1_files else None,
        'ft2': ft2_files[0] if ft2_files else None,
    }


def ensure_nonempty_gti(gti_file):
    """
    Check that a GTI FITS file contains at least one good time interval.

    Parameters
    ----------
    gti_file : str or Path
        Path to the GTI FITS file.

    Returns
    -------
    bool
        True if GTI has at least one interval.

    Raises
    ------
    ValueError
        If the GTI file is empty.
    """
    with fits.open(gti_file) as hdul:
        gti_ext = None
        for ext in hdul:
            if ext.name == 'GTI':
                gti_ext = ext
                break
        if gti_ext is None:
            raise ValueError(f"No GTI extension found in {gti_file}")
        n_gti = len(gti_ext.data)
        if n_gti == 0:
            raise ValueError(f"GTI file has zero intervals: {gti_file}")
    return True


def save_profile(filepath, flux_grid, logL_with, logL_null=None,
                 test_statistic=None, **metadata):
    """
    Save profile likelihood scan results to an NPZ file.

    Parameters
    ----------
    filepath : str or Path
        Output file path (.npz).
    flux_grid : array
        Flux values scanned [ph/cm2/s].
    logL_with : array
        Log-likelihood values with source at each flux.
    logL_null : float, optional
        Log-likelihood of null (no source) model.
    test_statistic : array, optional
        TS = 2 * (logL_with - logL_null) at each flux.
    **metadata
        Additional metadata to save (e.g., src_name, emin, emax).
    """
    save_dict = {
        'flux_grid': np.asarray(flux_grid),
        'logL_with': np.asarray(logL_with),
    }
    if logL_null is not None:
        save_dict['logL_null'] = np.asarray(logL_null)
    if test_statistic is not None:
        save_dict['test_statistic'] = np.asarray(test_statistic)

    for key, val in metadata.items():
        save_dict[key] = np.asarray(val)

    np.savez(filepath, **save_dict)


def load_profile(filepath):
    """
    Load profile likelihood scan results from an NPZ file.

    Parameters
    ----------
    filepath : str or Path
        Path to the .npz file.

    Returns
    -------
    dict
        Dictionary with keys: flux_grid, logL_with, logL_null,
        test_statistic, and any additional metadata.
    """
    data = np.load(filepath, allow_pickle=True)
    return dict(data)


def write_single_interval_gti(outpath, tstart, tstop, mjdref=51910.0):
    """
    Write a FITS GTI file with a single time interval.

    Parameters
    ----------
    outpath : str or Path
        Output file path.
    tstart, tstop : float
        Start and stop times [MET seconds].
    mjdref : float
        MJD reference for the Fermi mission (default 51910.0 = 2001-01-01).
    """
    col_start = fits.Column(name='START', format='D', array=[tstart])
    col_stop = fits.Column(name='STOP', format='D', array=[tstop])
    gti_hdu = fits.BinTableHDU.from_columns([col_start, col_stop])
    gti_hdu.name = 'GTI'
    gti_hdu.header['MJDREFI'] = int(mjdref)
    gti_hdu.header['MJDREFF'] = mjdref - int(mjdref)
    gti_hdu.header['TIMESYS'] = 'TT'
    gti_hdu.header['TIMEUNIT'] = 's'

    hdul = fits.HDUList([fits.PrimaryHDU(), gti_hdu])
    hdul.writeto(str(outpath), overwrite=True)
