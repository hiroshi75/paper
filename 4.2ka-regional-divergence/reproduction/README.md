# Reproduction Guide

## Regional Divergence in Population Responses to the 4.2 ka Climate Event

### Data Sources

| Database | Version | Access | N dates |
|----------|---------|--------|---------|
| p3k14c | v2024.1 | https://github.com/people3k/p3k14c | 180,070 |
| JOAD | v1.1.0 | https://doi.org/10.5334/joad.115 | 44,425 |
| NERD | — | Palmisano et al. (2021) supplement | — |

### Paleoclimate Proxy Records

| Record | Source |
|--------|--------|
| GISP2 ice core | NOAA Paleoclimatology (Alley, 2000) |
| Dongge Cave | Wang et al. (2005) |
| Gol-e-Zard Cave | Carolin et al. (2019) |

### Software Requirements

- R >= 4.0
- rcarbon >= 1.5 (Crema and Bevan, 2021)
- IntCal20 calibration curve (Reimer et al., 2020)

### Analysis Pipeline

1. **Data extraction**: Filter p3k14c by regional bounding boxes (Table 1 in paper). Filter JOAD by Japanese regional classification (Table 2).
2. **Quality control**: Exclude dates with error > 500 years, age <= 0, marine samples.
3. **Calibration**: All dates calibrated with IntCal20 via `calibrate()` in rcarbon.
4. **Site binning**: 200-year bins per site using `binPrep()` in rcarbon.
5. **SPD computation**: `spd()` with site-binned dates, normalized to unity.
6. **Population change**: Compare mean normalized SPD in pre-4.2 ka (5000-4500 cal BP) vs at-4.2 ka (4400-4000 cal BP) windows.
7. **Sensitivity analysis**: Repeat with centered, narrow, wide, and extended baseline windows.
8. **Statistical testing**:
   - `modelTest()`: Single-region test against exponential null (200 simulations).
   - `permTest()`: Mark permutation test between region pairs (500 permutations).
9. **Seed stability**: All permutation tests repeated with seeds 42, 123, 456, 789, 2024.

### Regional Bounding Boxes (p3k14c)

| Region | Lat (N) | Lon (E) |
|--------|---------|---------|
| Near East | 25-40 | 30-60 |
| Western Europe | 42-55 | -10-15 |
| Northern Europe | 55-70 | -10-30 |
| British Isles | 50-60 | -10-2 |
| East Asia (China) | 20-45 | 100-125 |
| Yellow River sub-region | >34 | 100-125 |
| Yangtze sub-region | 26-34 | 100-125 |

### Key Parameters

- Calibration curve: IntCal20
- Bin width: 200 years
- permTest permutations: 500
- modelTest simulations: 200
- Significance envelope: 95%

### Expected Key Results

| Comparison | p-value |
|------------|---------|
| Near East vs W. Europe (p3k14c) | 0.003 |
| Chubu vs Kyushu (JOAD) | 0.002 |
| Near East vs W. Europe (NERD) | 0.003 |
| Yellow River vs Yangtze | 0.71 (NS) |
