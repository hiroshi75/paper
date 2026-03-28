# Reproduction Guide

## Shared Holocene Trends Inflate Correlations Between Archaeological Population Proxies and Paleoecological Variables

### Data Sources

| Database | Version | Access | Usage |
|----------|---------|--------|-------|
| Neotoma Paleoecology DB | current | https://www.neotomadb.org / R neotoma2 | Pollen and charcoal records |
| p3k14c | v2024.1 | https://github.com/people3k/p3k14c | Radiocarbon dates for SPD |

### Regional Definitions

#### Pollen regions (Neotoma)

| Region | N sites | N samples | Criteria |
|--------|---------|-----------|----------|
| British Isles | 184 | 8,161 | UK + Ireland bounding box |
| Scandinavia | 290 | 16,369 | Norway, Sweden, Denmark, Finland |
| Central Europe | 57 | 4,074 | Germany, Austria, Switzerland, Czech Republic |

#### Charcoal regions (Neotoma)

| Region | N sites | Criteria |
|--------|---------|----------|
| Europe | 86 | European continent |
| North America | 55 | Continental US + southern Canada |

#### Radiocarbon regions (p3k14c)

Matching bounding boxes for each pollen/charcoal region above.

### Software Requirements

- R >= 4.0
- rcarbon >= 1.5 (Crema and Bevan, 2021)
- neotoma2 (R package for Neotoma API)
- mgcv (Wood, 2017) for GAM analyses
- IntCal20 calibration curve (Reimer et al., 2020)

### Analysis Pipeline

1. **Pollen data acquisition**: Query Neotoma via neotoma2 R package for each region. Filter for terrestrial pollen samples with count >= 150 grains.

2. **Rarefaction**: Compute rarefied taxonomic richness at n = 150 grains per sample using rarefaction (Birks and Line, 1992).

3. **Charcoal data acquisition**: Query Neotoma for charcoal accumulation records. Z-score normalize each site record independently.

4. **Radiocarbon SPD**: Extract p3k14c dates for matching regions. Calibrate with IntCal20. Apply site binning (h = 200 years) via `binPrep()`. Compute SPD via `spd()` in rcarbon.

5. **Binning**: Aggregate all variables into 500-year bins spanning 0-10,000 cal BP (20 bins).

6. **Raw Spearman correlation**: Compute Spearman's rho and p-value for each region-proxy pair.

7. **First-difference correlation**: First-difference both series (delta_t = x_{t+1} - x_t), yielding 19 pairs. Compute Spearman's rho on differenced series.

8. **GAM-residual correlation**: Fit y ~ s(time) and x ~ s(time) using mgcv. Extract residuals. Compute Spearman's rho on residual pairs.

9. **Additive GAM**: Fit y ~ s(time) + s(SPD) and compare to y ~ s(time) by AIC and deviance explained.

10. **Sensitivity test**: Repeat British Isles analysis with 200-year bins (50 bins, 49 first-difference pairs).

### Key Parameters

- Temporal window: 0-10,000 cal BP
- Bin width: 500 years (primary), 200 years (sensitivity)
- Rarefaction threshold: 150 grains
- Site binning: h = 200 years
- Calibration curve: IntCal20

### Expected Key Results

| Region | Proxy | Raw rho | Raw p | FD rho | FD p |
|--------|-------|---------|-------|--------|------|
| British Isles | Pollen richness | 0.735 | 0.0003 | 0.086 | 0.726 |
| Scandinavia | Pollen richness | 0.850 | <0.0001 | 0.454 | 0.052 |
| Central Europe | Pollen richness | -0.261 | 0.253 | 0.062 | 0.797 |
| Europe | Charcoal | -0.515 | 0.031 | 0.050 | 0.857 |
| North America | Charcoal | -0.743 | 0.0002 | -0.042 | 0.864 |
| British Isles (200yr) | Pollen richness | 0.691 | <0.001 | -0.003 | 0.984 |

### Supplementary Results

- British Isles GAM-residual rho: 0.005
- Scandinavia GAM-residual rho: 0.423
- British Isles additive GAM: SPD adds 0% deviance beyond s(time)
