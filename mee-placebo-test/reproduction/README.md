# Reproduction Guide: Placebo Test Framework for Citizen Science

This guide explains how to reproduce all six case studies from the placebo test framework paper.

---

## Prerequisites

### Software

| Requirement | Version | Check command |
|-------------|---------|---------------|
| Python | 3.10 or later | `python3 --version` |
| pip | any recent | `pip --version` |

### Install dependencies

Using conda:
```bash
conda env create -f environment.yml
conda activate mee-placebo-test
```

Or using pip:
```bash
pip install pandas numpy scipy scikit-learn statsmodels matplotlib seaborn requests rasterio elapid
```

---

## Quick Start

```bash
# 1. Download GBIF data for all case studies
cd reproduction/data/
python download_data.py --output-dir .

# 2. Run individual case studies
cd ../code/

# Case 1: Pollinator phenology
python pollinator_phenology_v2.py

# Case 2: Dawn chorus and ALAN
python soundscape_alan_analysis.py
python soundscape_bias_tests.py

# Case 3: Plant traits and drought
python plant_phenology_mismatch.py

# Case 4: Nightlight and insect communities
python nightlight_insects_analysis.py

# Case 5: Bird richness and ALAN
python ebird_alan_placebo.py

# Case 6: Forest fragmentation
python fragmentation_v2.py
```

---

## Case Studies

### Case 1: Pollinator Phenological Shift

**Script:** `pollinator_phenology_v2.py`

**Data:** 49,137 Apidae occurrence records from GBIF (Europe, 2000--2024)

**Analysis:**
- Calculate annual first-appearance day-of-year per latitude band
- Test latitudinal gradient in phenological shift
- Apply three effort correction methods

**Expected results:**
- Uncorrected: significant latitudinal gradient (R^2 = 0.79, p = 0.007)
- Corrected: gradient disappears (p = 0.60)
- Overall shift survives: -3.4 days/year (p < 10^-6)

### Case 2: Dawn Chorus and Artificial Light

**Scripts:** `soundscape_alan_analysis.py`, `soundscape_bias_tests.py`

**Data:** 8,073 dawn recordings from Xeno-canto + VIIRS nightlight radiance

**Analysis:**
- Correlate recording time vs. ALAN for dawn (04:00--08:00) and midday (10:00--14:00) windows
- Test observer fixed effects

**Expected results:**
- Dawn: r = -0.142 (p < 10^-37)
- Midday placebo: r = -0.208 (p < 10^-56) -- **stronger than dawn**
- With observer fixed effects: p = 0.90

### Case 3: Plant Traits and Drought

**Script:** `plant_phenology_mismatch.py`

**Data:** GBIF occurrence + WOODIV trait data for European trees

**Analysis:**
- Pilot (12 sites): CWM traits vs. drought response
- Expansion (46 sites): same model

**Expected results:**
- Pilot: R^2 = 0.895 (wood density positive)
- Expanded: R^2 = 0.483 (wood density sign reversal)

### Case 4: Nightlight and Insect Communities

**Script:** `nightlight_insects_analysis.py`

**Data:** GBIF Lepidoptera + Odonata (Europe) + VIIRS nightlight

**Analysis (placebo-first design):**
1. First run placebo: Odonata richness vs. ALAN
2. Then run target: Lepidoptera richness vs. ALAN

**Expected results:**
- Placebo (Odonata): beta = 3.21 (p < 0.0001)
- Target (Lepidoptera): beta = 1.32 (positive, opposite to prediction)

### Case 5: Bird Richness and ALAN

**Script:** `ebird_alan_placebo.py`

**Data:** GBIF/eBird Passeriformes + Strigiformes (Europe) + VIIRS

**Expected results:**
- eBird diurnal placebo: beta = 2.61 (p = 0.0005)
- eBird nocturnal target: beta = 0.66 (p < 0.00001)
- Placebo/main ratio: ~4x

### Case 6: Forest Fragmentation Thresholds

**Script:** `fragmentation_v2.py`

**Data:** GBIF birds + Hansen Global Forest Change + VIIRS

**Expected results:**
- Uncorrected: both target and placebo significant
- Effort-corrected: placebo non-significant; target retains tropical breakpoint at 13.6%

---

## Data Sources

| Source | URL | Used in |
|--------|-----|---------|
| GBIF | https://www.gbif.org | All cases |
| Xeno-canto | https://www.xeno-canto.org | Case 2 |
| eBird | https://ebird.org | Case 5 |
| VIIRS nightlight | https://eogdata.mines.edu/products/vnl/ | Cases 2, 4, 5, 6 |
| WorldClim 2.1 | https://www.worldclim.org | Cases 3, 4 |
| MODIS NDVI | https://lpdaac.usgs.gov | Case 3 |
| Hansen Forest Change | https://glad.earthengine.app | Case 6 |
| WOODIV traits | via WOODIV database | Case 3 |

---

## Key Script Descriptions

| Script | Case | Purpose |
|--------|------|---------|
| `pollinator_phenology_v2.py` | 1 | Pollinator phenological shift analysis with effort correction |
| `soundscape_alan_analysis.py` | 2 | Dawn chorus timing vs. ALAN analysis |
| `soundscape_bias_tests.py` | 2 | Placebo tests (midday window) and observer fixed effects |
| `plant_phenology_mismatch.py` | 3 | Plant trait--drought analysis with sample expansion |
| `nightlight_insects_analysis.py` | 4 | ALAN--insect richness with Odonata placebo |
| `ebird_alan_placebo.py` | 5 | Cross-platform bird richness--ALAN placebo |
| `fragmentation_v2.py` | 6 | Forest fragmentation thresholds with effort correction |

---

## Troubleshooting

### GBIF rate limits
The GBIF API has rate limits. If downloads fail, wait a few minutes and retry. Use `--case N` to download data for a single case study.

### Xeno-canto data
Audio recordings and metadata must be downloaded separately from https://www.xeno-canto.org using their API. The download script handles GBIF data only.

### VIIRS nightlight data
Annual composites must be downloaded manually from https://eogdata.mines.edu/products/vnl/. Select the VNL V2 annual composite for the study years.

### Small numerical differences
Results may vary slightly (+/-5%) across platforms due to floating-point arithmetic and random seed implementation differences.
