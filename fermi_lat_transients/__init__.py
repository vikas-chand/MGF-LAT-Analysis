"""
fermi_lat_transients
====================

A Python package for Fermi-LAT analysis of transient sources
(Magnetar Giant Flares, Fast X-ray Transients, GRBs, etc.).

Provides a modular pipeline for:
- Event selection and data reduction (pipeline)
- Unbinned likelihood analysis (likelihood)
- Profile likelihood scanning and upper limits (likelihood.profile)
- SED construction (likelihood.sed)
- Stacking analysis (stacking)
- Physics calculations and flux conversions (physics)
- Publication-quality visualization (plotting)

Quick Start
-----------
>>> from fermi_lat_transients.config import load_config
>>> from fermi_lat_transients.data import load_catalog
>>> cfg = load_config('inputs.yaml')
>>> targets = load_catalog(cfg.catalog_file)
"""

__version__ = "0.1.0"
__author__ = "Vikas Chand"

from .config import load_config, Config
