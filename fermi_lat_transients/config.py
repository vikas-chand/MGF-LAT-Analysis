"""
Configuration management for Fermi-LAT transient analysis.

Loads analysis parameters from a YAML file and provides defaults.
"""

import os
import yaml
import copy
from pathlib import Path

# ─── Default configuration ──────────────────────────────────────────────────

DEFAULTS = {
    # Source catalog
    "catalog_file": None,
    "catalog_format": "tab",          # tab, csv, fits
    "column_name": "GRB_name",
    "column_ra": "gal_ra_deg",
    "column_dec": "gal_dec_deg",
    "column_trigger_time": "MET_trig_time",

    # Fermi data paths
    "ft1": None,                      # Event file or filelist
    "ft2": None,                      # Spacecraft file
    "galdiff": None,                  # Galactic diffuse model FITS
    "isodiff": None,                  # Isotropic diffuse template
    "ltcube": None,                   # Pre-computed livetime cube (optional)
    "catalog_fits": None,             # 4FGL catalog FITS file

    # LAT analysis parameters
    "emin": 100.0,                    # MeV
    "emax": 100000.0,                 # MeV
    "roi": 12.0,                      # degrees
    "zmax": 100.0,                    # degrees
    "evclass": 8,                     # 8=Transient, 128=Source
    "evtype": 3,                      # FRONT+BACK
    "irfs": "P8R3_TRANSIENT020E_V3",

    # Time window relative to trigger
    "dt_pre": 0.0,                    # seconds before trigger
    "dt_post": 500.0,                 # seconds after trigger

    # Profile likelihood scanning
    "spectral_index": -2.0,
    "flux_npts": 150,
    "flux_min": 1.0e-11,             # ph/cm2/s
    "flux_max": 1.0e-3,              # ph/cm2/s
    "delta_ts_95": 2.71,

    # SED energy bins (MeV, bin edges)
    "sed_bins": [100, 316, 1000, 3162, 10000],

    # Transient source model
    "transient_src_name": "Transient",
    "transient_spectrum": "PowerLaw2",

    # Output
    "output_dir": "output",

    # Optimizer
    "optimizer": "NEWMINUIT",
    "fit_tolerance": 1e-5,
}


class Config(dict):
    """
    Configuration container that behaves like a dict but also supports
    attribute access (cfg.emin instead of cfg['emin']).

    Parameters
    ----------
    yaml_path : str or Path, optional
        Path to a YAML configuration file. Values override defaults.
    **overrides : dict
        Additional keyword overrides applied after YAML loading.
    """

    def __init__(self, yaml_path=None, **overrides):
        super().__init__(copy.deepcopy(DEFAULTS))

        if yaml_path is not None:
            yaml_path = Path(yaml_path)
            if not yaml_path.exists():
                raise FileNotFoundError(f"Config file not found: {yaml_path}")
            with open(yaml_path, 'r') as f:
                user_cfg = yaml.safe_load(f) or {}
            self.update(user_cfg)

        # Apply any programmatic overrides
        self.update(overrides)

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"Config has no parameter '{key}'")

    def __setattr__(self, key, value):
        self[key] = value

    def save(self, path):
        """Save current configuration to a YAML file."""
        path = Path(path)
        with open(path, 'w') as f:
            yaml.dump(dict(self), f, default_flow_style=False, sort_keys=False)

    def validate(self):
        """Check that required paths are set."""
        issues = []
        if self.get("ft2") is None:
            issues.append("ft2 (spacecraft file) must be specified")
        if self.get("galdiff") is None:
            issues.append("galdiff (Galactic diffuse model) must be specified")
        if self.get("isodiff") is None:
            issues.append("isodiff (isotropic diffuse template) must be specified")
        if issues:
            raise ValueError("Configuration issues:\n  " + "\n  ".join(issues))
        return True


def load_config(yaml_path=None, **overrides):
    """
    Convenience function to load a Config object.

    Parameters
    ----------
    yaml_path : str or Path, optional
        Path to YAML config file.
    **overrides
        Additional parameter overrides.

    Returns
    -------
    Config
    """
    return Config(yaml_path=yaml_path, **overrides)
