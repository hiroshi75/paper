# Reproduction Guide: Thermal Niche Asymmetry and Climate Refugia

This guide walks you through reproducing the complete climate refugia analysis from scratch. Each step explains what the code does, why it matters, and what you should see when it finishes.

---

## Prerequisites

### Software

| Requirement | Version | Check command |
|-------------|---------|---------------|
| Python | 3.10 or later | `python3 --version` |
| pip | any recent | `pip --version` |
| git | any recent | `git --version` |

### Install dependencies

Using conda:
```bash
conda env create -f environment.yml
conda activate refugia-asymmetry
```

Or using pip:
```bash
pip install pandas numpy scipy scikit-learn matplotlib rasterio elapid requests
```

### Estimated time and disk space

| Step | Time | Disk |
|------|------|------|
| GBIF data download | ~10 min | ~5 MB |
| WorldClim download | ~5 min | ~200 MB |
| SDM analysis (all 5 pairs) | ~20-30 min | ~2 MB |
| Figure generation | ~2 min | ~5 MB |
| **Total** | **~40-50 min** | **~300 MB** |

---

## Quick Start

```bash
# 1. Download data
cd reproduction/data/
python download_data.py --output-dir .

# 2. Run the main analysis scripts
cd ../code/

# Build ensemble SDMs for European pairs
python refugia_maxent_ensemble.py

# Run SSP scenario projections
python refugia_ssp_ensemble.py

# Run global validation pairs
python refugia_global_expansion.py
```

---

## Step-by-Step Walkthrough

### Step 1: Data Acquisition

#### 1.1 GBIF Species Occurrence Data

Downloads observation records for 9 European tree species from GBIF (https://www.gbif.org/). Each record contains species name, latitude, longitude, and year of observation.

```bash
cd reproduction/data/
python download_data.py
```

**Expected output:** 9 CSV files (one per species), each with 2,800-4,000+ records.

#### 1.2 WorldClim Climate Data

Downloads WorldClim 2.1 bioclimatic variables at 10 arc-minute resolution. Four variables are used:
- BIO1: Mean Annual Temperature (C)
- BIO5: Max Temperature of Warmest Month (C)
- BIO6: Min Temperature of Coldest Month (C)
- BIO12: Annual Precipitation (mm)

### Step 2: Species Distribution Modelling

The analysis uses a 10-algorithm ensemble SDM approach:
- GLM, GAM, GBM, CTA, ANN, SRE, FDA, MARS, RF, MaxEnt
- AUC-weighted consensus projections
- Spatial block cross-validation (5 longitudinal blocks)
- 1:1 presence:background ratio

```bash
python refugia_maxent_ensemble.py
```

**Expected model performance (Ensemble AUC):**

| Species | AUC |
|---------|-----|
| *Quercus robur* | ~0.815 |
| *Fagus sylvatica* | ~0.829 |
| *Picea abies* | ~0.882 |
| *Prunus serotina* | ~0.910 |
| *Quercus ilex* | ~0.849 |

### Step 3: Climate Projections

CMIP6 SSP2-4.5 and SSP5-8.5 scenarios (2061-2080) with three GCMs:
- ACCESS-CM2
- IPSL-CM6A-LR
- MRI-ESM2-0

```bash
python refugia_ssp_ensemble.py
```

### Step 4: Refugia Quantification

Refugia are defined as grid cells where:
- Native species suitability >= 0.5
- Invasive species suitability < 0.5

**Expected key results under SSP5-8.5:**

| Pair | Refugia loss |
|------|-------------|
| *F. sylvatica* vs *A. altissima* | -89.9% |
| *Q. robur* vs *R. pseudoacacia* | -83.9% |
| *P. abies* vs *P. serotina* | -65.9% |
| *Q. ilex* vs *P. americana* | -41.2% |
| *Q. ilex* vs *I. glandulifera* | **+45.6%** |

### Step 5: Global Validation

Five additional pairs from four continents:

```bash
python refugia_global_expansion.py
```

---

## Key Scripts

| Script | Purpose |
|--------|---------|
| `refugia_sdm_v2.py` | Core SDM methodology (RF + MaxEnt) |
| `refugia_multi_species.py` | Multi-species pair analysis |
| `refugia_quadrant_completion.py` | Complete all 4 quadrants |
| `refugia_maxent_ensemble.py` | 10-algorithm ensemble SDM |
| `refugia_ssp_ensemble.py` | CMIP6 SSP scenario projections |
| `sdm_ensemble.py` | Ensemble model utilities |
| `refugia_global_expansion.py` | Global validation (5 continents) |

---

## Troubleshooting

### `ModuleNotFoundError: No module named 'elapid'`
```bash
pip install elapid
```
If this fails, the analysis will still run with Random Forest only.

### `ModuleNotFoundError: No module named 'rasterio'`
Requires GDAL. On Ubuntu/Debian:
```bash
sudo apt-get install gdal-bin libgdal-dev
pip install rasterio
```

### GBIF download timeout
If GBIF is temporarily down, wait and retry. You can also manually download from https://www.gbif.org/.

### Small numerical differences
Differences of +/-5% on loss percentages are normal due to random seed implementation differences across platforms.

---

## Data Sources

- **GBIF**: https://www.gbif.org (species occurrence records)
- **WorldClim 2.1**: https://www.worldclim.org (current climate)
- **CMIP6**: WorldClim 2.1 CMIP6 archive (future projections)
