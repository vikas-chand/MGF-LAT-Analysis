"""
Target/catalog handling for transient source analysis.

Provides functions to load source catalogs (tab, CSV, FITS),
build target dictionaries, and convert between time systems.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from astropy.time import Time

from ..constants import FERMI_EPOCH


def utc_to_met(t0, time_kind="auto"):
    """
    Convert a time value to Fermi Mission Elapsed Time (MET).

    Parameters
    ----------
    t0 : str or float
        Time value. If string, interpreted as UTC (ISO format).
        If float, interpreted as MET directly (unless time_kind='utc').
    time_kind : str
        'auto' (default), 'utc', or 'met'.
        - 'auto': string -> UTC, number -> MET
        - 'utc': force UTC interpretation
        - 'met': force MET interpretation

    Returns
    -------
    float
        Fermi MET [seconds].
    """
    if time_kind == "met":
        return float(t0)

    if time_kind == "auto":
        if isinstance(t0, str):
            time_kind = "utc"
        else:
            return float(t0)

    if time_kind == "utc":
        t = Time(t0, scale='utc')
        met = (t - FERMI_EPOCH).sec
        return met

    raise ValueError(f"Unknown time_kind: {time_kind!r}")


def build_target(name, ra_deg, dec_deg, t0, dt_pre=0.0, dt_post=500.0,
                 time_kind="auto", **extra):
    """
    Build a target dictionary for analysis.

    Parameters
    ----------
    name : str
        Source name (e.g., 'GRB_200415A').
    ra_deg, dec_deg : float
        Source coordinates [degrees].
    t0 : str or float
        Trigger time (UTC string or MET float).
    dt_pre : float
        Seconds before trigger to include (default 0).
    dt_post : float
        Seconds after trigger to include (default 500).
    time_kind : str
        Time interpretation ('auto', 'utc', 'met').
    **extra
        Additional metadata (e.g., host_galaxy, distance_Mpc).

    Returns
    -------
    dict
        Target dictionary with keys: name, ra, dec, t0_met, tmin, tmax,
        plus any extra metadata.
    """
    met = utc_to_met(t0, time_kind)
    target = {
        'name': name,
        'ra': float(ra_deg),
        'dec': float(dec_deg),
        't0_met': met,
        'tmin': met - dt_pre,
        'tmax': met + dt_post,
    }
    target.update(extra)
    return target


def load_catalog(filepath, fmt="auto", column_name="GRB_name",
                 column_ra="gal_ra_deg", column_dec="gal_dec_deg",
                 column_trigger="MET_trig_time",
                 dt_pre=0.0, dt_post=500.0, **extra_columns):
    """
    Load a source catalog and return a list of target dictionaries.

    Parameters
    ----------
    filepath : str or Path
        Path to catalog file.
    fmt : str
        File format: 'auto' (detect from extension), 'tab', 'csv', 'fits'.
    column_name : str
        Column name for source name.
    column_ra, column_dec : str
        Column names for RA, Dec [degrees].
    column_trigger : str
        Column name for trigger time.
    dt_pre, dt_post : float
        Time window before/after trigger [seconds].
    **extra_columns
        Mapping of target_key -> column_name for additional metadata.

    Returns
    -------
    list of dict
        List of target dictionaries.
    """
    filepath = Path(filepath)

    if fmt == "auto":
        ext = filepath.suffix.lower()
        if ext in ('.csv',):
            fmt = "csv"
        elif ext in ('.fits', '.fit'):
            fmt = "fits"
        else:
            fmt = "tab"

    if fmt == "csv":
        df = pd.read_csv(filepath)
    elif fmt == "fits":
        from astropy.table import Table
        t = Table.read(filepath)
        df = t.to_pandas()
    else:  # tab-separated
        df = pd.read_csv(filepath, sep=r'\s+')

    targets = []
    for _, row in df.iterrows():
        extra = {}
        for target_key, col_name in extra_columns.items():
            if col_name in row:
                extra[target_key] = row[col_name]

        target = build_target(
            name=str(row[column_name]),
            ra_deg=float(row[column_ra]),
            dec_deg=float(row[column_dec]),
            t0=row[column_trigger],
            dt_pre=dt_pre,
            dt_post=dt_post,
            time_kind="auto",
            **extra,
        )
        targets.append(target)

    return targets
