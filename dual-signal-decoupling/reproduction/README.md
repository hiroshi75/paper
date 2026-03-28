# Reproduction Guide

## Data
All pollen data from Neotoma Paleoecology Database (public). Cached JSON files used for reproducibility.

## Requirements
R >= 4.3 with neotoma2, vegan, mgcv, segmented

## Scripts
- `01_lag_estimation.R` — Bootstrap lag estimation (1000 iterations)
- `02_site_replication.R` — Site-level classification (indicator-first vs AP-first)
- `03_threshold_gam.R` — GAM and segmented regression threshold analysis
- `run_all_cached.R` — Complete pipeline using cached data
