# Reproduction Guide

## Overview
Reproduces analyses in "Detecting Transformation, Not Use: Pollen Exceedance Distinguishes European Pastoral Agriculture from Indigenous Agroforestry Across 442 Sites."

## Requirements
- R >= 4.3
- R packages: neotoma2, vegan, parallel

## Data
All pollen data from Neotoma Paleoecology Database (public). Pre-computed results in `data/`.

## Scripts (code/)
- `check_eastern_na_pollen.R` — Data availability assessment (209 sites, Zea mays detection)
- `pilot_ena_exceedance.R` — Pilot analysis (96 sites)
- `full_ena_exceedance.R` — Full analysis (111 sites, all results)
- `ena_regional_stats.R` — Regional stratification + statistical tests
- `ena_rpp_uncertainty.R` — RPP Monte Carlo (10,000 iterations)
- `ena_charcoal_analysis.R` — Charcoal corroboration (22 sites)

## Execution
```bash
Rscript code/full_ena_exceedance.R    # Main analysis
Rscript code/ena_regional_stats.R      # Regional + stats
Rscript code/ena_rpp_uncertainty.R     # RPP Monte Carlo
```
