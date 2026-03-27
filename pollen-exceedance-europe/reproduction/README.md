# Reproduction Guide

## Overview

This package reproduces the analyses in "Pastoral Indicator Exceedance, Not First Appearance, Distinguishes Anthropogenic from Natural Pollen Signals: Evidence from 331 Sites in Northwestern and Central Europe."

## Requirements

- R >= 4.3
- R packages: neotoma2, vegan, spdep, lme4, lmerTest, metafor, pwr, sf

## Data

All pollen data are publicly available from the Neotoma Paleoecology Database (https://www.neotomadb.org/). The analysis scripts use the `neotoma2` R package to query data directly from the API.

Pre-computed intermediate data files (signal phase detection, decomposition, classification results) are available in `data/`.

## Analysis Scripts

Scripts in `code/` are organized by analysis phase:

### Phase 1: Core analyses
- `phase1_circularity_removal.R` — Time-window-free reclassification + independent Neolithic concordance test
- `phase1_trigger_falsification.R` — H1 vs H3 pastoral BC% comparison
- `phase1_spatial_autocorrelation.R` — Moran's I tests at 50km and 100km
- `phase1_amplification_investigation.R` — Component trade-off analysis (replacement pattern)

### Phase 2: Robustness and sensitivity
- `phase2_power_analysis.R` — Power analysis for magnitude comparisons (Cohen's d scenarios)
- `phase2_independent_classification.R` — Time-based and geographic classification tests
- `phase2_poaceae_sensitivity.R` — Poaceae exclusion sensitivity
- `phase2_mixed_effects.R` — LMM and meta-analysis with region random effect
- `phase2_lag_investigation.R` — Lag distribution and decomposition
- `phase2_alt_classifications.R` — Alternative classification approaches

## Execution

```bash
# Install R dependencies
Rscript -e 'install.packages(c("neotoma2","vegan","spdep","lme4","lmerTest","metafor","pwr","sf"), repos="https://cloud.r-project.org")'

# Run individual analyses (each script is self-contained)
Rscript code/phase1_circularity_removal.R
Rscript code/phase2_power_analysis.R
# etc.
```

## Expected Output

Each script produces CSV results and/or console output matching the statistics reported in the paper (Tables 2-9, Section 3 results).
