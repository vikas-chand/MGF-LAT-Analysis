# MGF LAT Analysis

Fermi-LAT GeV analysis of extragalactic Magnetar Giant Flare (MGF) candidates.

## Overview

This repository contains the Fermi-LAT analysis pipeline and results for 7 extragalactic MGF candidates detected by Fermi-GBM. The analysis includes individual source likelihood analysis, TS map generation, spectral energy distribution (SED) construction, and a stacking analysis combining all candidates to derive the most constraining upper limits on delayed/extended GeV emission from MGFs.

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

## Repository Structure

```
MGF_LAT_Analysis/
├── trigger_analysis/           # Individual GRB LAT likelihood analysis
│   ├── MGF_analysis_LAT.ipynb  # Master analysis notebook
│   ├── MGF_GRB_200415A_LAT_v1.ipynb  # Detailed GRB 200415A analysis
│   ├── MGF_Discussion.ipynb    # Fireball parameter calculations
│   ├── GBM_eMGF_candidates.txt # Source catalog
│   ├── make4FGLxml.py          # 4FGL XML model builder
│   ├── GRB_*/                  # Per-source analysis & results
│   │   ├── MGF_GRB_*_LAT_v1.ipynb
│   │   ├── *_Analysis/         # XML models, TS maps, ROI regions
│   │   ├── source_TS.txt       # TS values
│   │   └── trigcat.txt         # Trigger catalog info
│   ├── gtburst_results/        # GTBURST upper limits & flux plots
│   └── exposure_check/         # LAT exposure verification
│
├── stacking_analysis/          # Combined stacking of all GRBs
│   ├── MGF_Stacking_Analysis.ipynb
│   ├── MGF_Stacking_Analysis_with_IndividualULs.ipynb
│   ├── MGF_SED_Stacking_Analysis.ipynb
│   ├── MGF_Analysis_SED.ipynb
│   ├── MGF_Discussion_Stacking.ipynb
│   ├── Calculate_TstartGTI.ipynb
│   ├── upper_limits_summary.csv
│   ├── GRB_*/                  # Per-source stacking data
│   └── *.csv / *.pdf           # SED points & plots
│
├── stacking_analysis_offset/   # Control regions (offset stacking)
│
├── rates_analysis/             # MGF rate calculations
│   └── GRB_Rates_Analytical.ipynb
│
├── utils/                      # Shared utility scripts
│   ├── ts_maps/                # TS map generation & publication plots
│   │   └── plot_tsmap_publication.py
│   └── sed_scripts/            # SED construction (likeSED v13.1)
│
└── figures/                    # Publication-ready figures
```

## Analysis Pipeline

1. **Data Selection** - Fermi-LAT Pass 8 SOURCE class events, 100 MeV - 300 GeV
2. **Likelihood Analysis** - Unbinned likelihood with 4FGL catalog sources
3. **TS Maps** - Test Statistic maps at each GRB position
4. **Upper Limits** - 95% CL flux upper limits per source
5. **SED Construction** - Multi-band spectral energy distributions
6. **Stacking** - Joint likelihood combining all candidates
7. **Fireball Physics** - Derived parameters (L0, T0, eta) from scaling relations

## Dependencies

- Fermitools (gtselect, gtmktime, gtltcube, gtexpmap, gtlike, gttsmap, gtsrcprob)
- Python 3, astropy, numpy, scipy, pandas, matplotlib

## Data

Large FITS data files are not tracked by git. To reproduce the analysis, download Fermi-LAT data from the [FSSC Data Server](https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi) using the coordinates and time windows specified in each notebook.
