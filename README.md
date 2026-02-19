# fermi_lat_transients

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

**Fermi-LAT unbinned likelihood analysis and stacking pipeline for transient sources.**

A reusable Python package for GeV gamma-ray analysis of any transient source population (Magnetar Giant Flares, Fast X-ray Transients, GRBs, etc.) using Fermi-LAT data.

## Overview

This package provides a complete pipeline for:
- **Event selection & data reduction** via Fermitools wrappers
- **Unbinned likelihood analysis** with 4FGL catalog models
- **Profile likelihood scanning** and 95% CL upper limit extraction
- **SED construction** via per-energy-bin profile scans
- **Stacking analysis** combining profile likelihoods across multiple sources
- **Physics calculations**: flux conversions, fireball parameters, cosmological rates
- **Publication-quality visualization**: TS maps, profile plots, SED plots

Currently applied to 7 extragalactic MGF candidates detected by Fermi-GBM.

## Installation

```bash
# From source (recommended for development)
git clone https://github.com/vikas-chand/MGF-LAT-Analysis.git
cd MGF-LAT-Analysis
pip install -e .

# Or install dependencies only
pip install -r requirements.txt
```

**Note:** [Fermitools](https://github.com/fermi-lat/Fermitools-conda/) must be installed separately (conda-based).

## Quick Start

### Using the configuration-driven pipeline

```bash
# 1. Edit inputs.yaml with your data paths and parameters
# 2. Run the full pipeline
lat-transients-run inputs.yaml
```

### Using the package in Python

```python
from fermi_lat_transients.config import load_config
from fermi_lat_transients.data import load_catalog, build_target
from fermi_lat_transients.pipeline import run_pipeline
from fermi_lat_transients.likelihood import profile_scan, upper_limit_from_profile
from fermi_lat_transients.stacking import stack_profiles, stacked_upper_limit
from fermi_lat_transients.physics import isotropic_energy, K_from_photon_flux

# Load config and catalog
cfg = load_config('inputs.yaml')
targets = load_catalog(cfg.catalog_file)

# Or build a single target manually
target = build_target('MySource', ra=180.0, dec=45.0, t0=600000000, dt_post=1000)
```

## Configuration

All analysis parameters are controlled via `inputs.yaml`:

```yaml
# Source catalog
catalog_file: "data/catalogs/my_sources.csv"
catalog_format: "csv"

# Fermi data paths
ft2: "/path/to/spacecraft.fits"
galdiff: "/path/to/gll_iem_v07.fits"
isodiff: "/path/to/iso_P8R3_SOURCE_V3_v1.txt"
catalog_fits: "/path/to/gll_psc_v35.fit"

# Analysis parameters
emin: 100.0         # MeV
emax: 100000.0      # MeV
roi: 12.0           # degrees
irfs: "P8R3_TRANSIENT020E_V3"
dt_post: 500.0      # seconds after trigger
spectral_index: -2.0
```

See [`fermi_lat_transients/inputs.yaml`](fermi_lat_transients/inputs.yaml) for the full template.

## Repository Structure

```
MGF-LAT-Analysis/
├── fermi_lat_transients/       # Python package
│   ├── pipeline/               #   Fermitools wrappers (gtselect, gtmktime, ...)
│   ├── likelihood/             #   Unbinned likelihood, profile scans, SED
│   ├── stacking/               #   Profile & SED stacking
│   ├── physics/                #   Flux conversions, fireball params, rates
│   ├── data/                   #   Catalog loading, file I/O
│   ├── plotting/               #   TS maps, profiles, SEDs
│   ├── client.py               #   Main pipeline orchestration
│   ├── config.py               #   YAML configuration system
│   ├── constants.py            #   Physical constants
│   └── inputs.yaml             #   Default config template
│
├── notebooks/                  # Analysis notebooks (examples)
│   └── mgf/                    #   MGF-specific analysis
│
├── data/                       # Source catalogs & results
│   ├── catalogs/               #   Source lists
│   ├── results/                #   Per-source analysis outputs
│   ├── stacking/               #   Stacking results
│   └── stacking_offset/        #   Control region results
│
├── figures/                    # Publication figures
├── setup.py                    # Package installation
├── requirements.txt            # Dependencies
└── README.md
```

## Package API

| Module | Description |
|--------|-------------|
| `fermi_lat_transients.pipeline` | Fermitools wrappers: `gtselect`, `gtmktime`, `gtltcube`, `gtexpmap`, `gtdiffrsp` |
| `fermi_lat_transients.likelihood` | `setup_unbinned`, `fit`, `profile_scan`, `upper_limit_from_profile`, `sed_profile_scan` |
| `fermi_lat_transients.stacking` | `stack_profiles`, `stacked_upper_limit`, `stack_sed_profiles` |
| `fermi_lat_transients.physics` | `K_from_photon_flux`, `energy_flux_from_K`, `isotropic_energy`, `L0_from_Lgamma`, fireball params |
| `fermi_lat_transients.data` | `load_catalog`, `build_target`, `utc_to_met`, `find_ft1_ft2` |
| `fermi_lat_transients.plotting` | `plot_profile`, `plot_stacked_profile`, `plot_sed`, `plot_stacked_sed` |
| `fermi_lat_transients.config` | `load_config`, `Config` class |

## MGF Candidates

| GRB | Host Galaxy | Trigger MET |
|-----|------------|-------------|
| GRB 081213A | NGC 0253 | 250834182 |
| GRB 120616A | IC 0342 | 361552012 |
| GRB 180128A | NGC 0253 | 538809001 |
| GRB 200415A | NGC 0253 | 608633290 |
| GRB 200423A | NGC 6946 | 609342856 |
| GRB 231024A | NGC 0253 | 719846439 |
| GRB 231115A | M82 | 721755386 |

## Analysis Pipeline

1. **Data Selection** - Fermi-LAT Pass 8 events, 100 MeV - 100 GeV
2. **Likelihood Analysis** - Unbinned likelihood with 4FGL catalog sources
3. **TS Maps** - Test Statistic maps at each source position
4. **Profile Likelihood** - Flux scanning and 95% CL upper limits
5. **SED Construction** - Multi-band spectral energy distributions
6. **Stacking** - Combined profile likelihood across all sources
7. **Physics** - Fireball parameters, isotropic energies, rate constraints

## Dependencies

- **Fermitools** (gtselect, gtmktime, gtltcube, gtexpmap, gtlike) - install via conda
- Python 3.8+, numpy, pandas, astropy, matplotlib, scipy, pyyaml

## Data

Large FITS data files are not tracked by git. Download Fermi-LAT data from the [FSSC Data Server](https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi).

## License

This project is for scientific research purposes.
