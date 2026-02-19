# Analysis Notebooks

Jupyter notebooks demonstrating the `fermi_lat_transients` package applied to
extragalactic Magnetar Giant Flare (MGF) candidates.

## MGF Analysis Notebooks

| Notebook | Description |
|----------|-------------|
| `MGF_analysis_LAT.ipynb` | Master single-source LAT likelihood analysis template |
| `MGF_GRB_200415A_LAT_v1.ipynb` | Detailed analysis of GRB 200415A (detected source) with TS maps |
| `MGF_Discussion.ipynb` | Fireball scaling relations (L0, T0, eta) for all candidates |
| `MGF_Stacking_Analysis.ipynb` | Profile likelihood stacking of all 7 MGF candidates |
| `MGF_Stacking_Analysis_with_IndividualULs.ipynb` | Individual + stacked 95% CL upper limits |
| `MGF_SED_Stacking_Analysis.ipynb` | Per-energy-bin SED stacking analysis |
| `MGF_Analysis_SED.ipynb` | Single-event SED construction (GRB 200415A) |
| `MGF_Discussion_Stacking.ipynb` | Physics discussion: flux conversions, isotropic energies |
| `MGF_UpperLimits_Plots.ipynb` | Publication upper limit comparison plots |
| `Calculate_TstartGTI.ipynb` | GTI extraction and exposure verification |
| `GRB_Rates_Analytical.ipynb` | Cosmological volume and MGF rate calculations |

## Usage

These notebooks serve as examples for the `fermi_lat_transients` package.
To adapt for a different source population (e.g., FXTs), create a new
subdirectory and modify the catalog and configuration.
