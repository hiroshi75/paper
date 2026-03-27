#!/usr/bin/env python3
"""
Fat-Tailed Biodiversity Trends — Full-Scale Analysis

Extends pilot (n=20) to n=100 European grid cells with:
1. Effort correction (rarefied richness change)
2. Subbotin distribution fitting (Normal and Laplace are special cases)
3. Spatial predictors of extreme changes (quantile analysis)
4. Multi-panel publication figure

Tests McGill & Magurran (2026, Ecology Letters) finding that biodiversity
trends are leptokurtic, and extends by asking WHAT PREDICTS extreme changes.
"""

import json
import os
import time
import urllib.request
import urllib.error
import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

GBIF_BASE = "https://api.gbif.org/v1/occurrence/search"
TAXON_KEY = 212  # Aves
LIMIT = 300

PERIOD_EARLY = ("2000-01-01", "2010-12-31")
PERIOD_LATE = ("2014-01-01", "2024-12-31")

# 100 grid cells across Europe (lat 38-62N, lon -10 to 30E, every 2.5 degrees)
GRID_CELLS = []
for lat in np.arange(38, 62, 2.5):
    for lon in np.arange(-10, 30, 4):
        GRID_CELLS.append((float(lat), float(lon)))

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
CACHE_PATH = os.path.join(OUTPUT_DIR, "fat_tailed_full_cache.csv")
FIG_PATH = os.path.join(OUTPUT_DIR, "fig_fat_tailed_full.png")
RESULTS_PATH = os.path.join(OUTPUT_DIR, "fat_tailed_full_results.md")
DATA_PATH = os.path.join(OUTPUT_DIR, "fat_tailed_full_data.csv")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def query_gbif(lat_min, lat_max, lon_min, lon_max, event_start, event_end):
    """Query GBIF occurrence API. Returns list of records."""
    params = (
        f"taxonKey={TAXON_KEY}"
        f"&hasCoordinate=true"
        f"&decimalLatitude={lat_min},{lat_max}"
        f"&decimalLongitude={lon_min},{lon_max}"
        f"&eventDate={event_start},{event_end}"
        f"&limit={LIMIT}"
    )
    url = f"{GBIF_BASE}?{params}"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode())
            return data.get("results", []), data.get("count", 0)
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
            print(f"    [warn] attempt {attempt+1} failed: {exc}")
            time.sleep(2 * (attempt + 1))
    return [], 0


def unique_species(records):
    """Return set of unique speciesKey values."""
    return {r["speciesKey"] for r in records if r.get("speciesKey") is not None}


def rarefied_richness(records, n_subsample=100, n_iter=20):
    """Estimate rarefied richness by subsampling n records, repeated n_iter times."""
    species_list = [r.get("speciesKey") for r in records if r.get("speciesKey") is not None]
    if len(species_list) <= n_subsample:
        return len(set(species_list))
    richnesses = []
    for _ in range(n_iter):
        sub = np.random.choice(species_list, size=n_subsample, replace=False)
        richnesses.append(len(set(sub)))
    return np.mean(richnesses)


def subbotin_logpdf(x, mu, sigma, beta):
    """Log-PDF of Subbotin (generalized Normal) distribution.

    beta=2 -> Normal, beta=1 -> Laplace, beta<1 -> heavier tails.
    """
    from scipy.special import gamma as gamma_fn
    z = np.abs(x - mu) / sigma
    log_norm = np.log(beta) - np.log(2 * sigma) - np.log(gamma_fn(1 + 1/beta))
    return log_norm - z**beta


def fit_subbotin(data):
    """Fit Subbotin distribution via MLE. Returns (mu, sigma, beta, loglik)."""
    mu_init = np.median(data)
    sigma_init = np.std(data)

    best_ll = -np.inf
    best_params = (mu_init, sigma_init, 2.0)

    # Grid search over beta, optimize mu and sigma analytically
    for beta in np.arange(0.3, 4.1, 0.1):
        # For given beta, MLE of mu is the value minimizing sum |x-mu|^beta
        def neg_ll_mu(mu):
            return np.sum(np.abs(data - mu)**beta)

        res = minimize_scalar(neg_ll_mu, bounds=(np.min(data), np.max(data)), method='bounded')
        mu = res.x
        sigma = (np.mean(np.abs(data - mu)**beta))**(1/beta)

        if sigma <= 0:
            continue

        ll = np.sum(subbotin_logpdf(data, mu, sigma, beta))
        if ll > best_ll:
            best_ll = ll
            best_params = (mu, sigma, beta)

    return (*best_params, best_ll)


def aic(k, ll):
    return 2 * k - 2 * ll


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    np.random.seed(42)

    print(f"Fat-Tailed Biodiversity Trends — Full Analysis")
    print(f"Grid cells: {len(GRID_CELLS)}")
    print(f"Periods: {PERIOD_EARLY} vs {PERIOD_LATE}")
    print()

    # Check for cache
    results = []
    if os.path.exists(CACHE_PATH):
        print(f"Loading cache from {CACHE_PATH}")
        with open(CACHE_PATH) as f:
            reader = csv.DictReader(f)
            for row in reader:
                results.append({k: float(v) if k != 'cell_id' else int(v)
                               for k, v in row.items()})
        print(f"Loaded {len(results)} cached cells")

    cached_coords = {(r['lat'], r['lon']) for r in results}
    remaining = [(lat, lon) for lat, lon in GRID_CELLS if (lat, lon) not in cached_coords]

    if remaining:
        print(f"Querying {len(remaining)} remaining cells...")

    for idx, (lat, lon) in enumerate(remaining):
        cell_id = len(results) + 1
        label = f"Cell {cell_id}/{len(GRID_CELLS)} (lat={lat:.1f}, lon={lon:.1f})"
        print(f"[{label}] querying...")

        early_recs, early_count = query_gbif(lat, lat + 1, lon, lon + 1, *PERIOD_EARLY)
        time.sleep(0.3)
        late_recs, late_count = query_gbif(lat, lat + 1, lon, lon + 1, *PERIOD_LATE)
        time.sleep(0.3)

        sp_early = unique_species(early_recs)
        sp_late = unique_species(late_recs)
        s1 = len(sp_early)
        s2 = len(sp_late)

        # Effort-corrected (rarefied) richness
        s1_rare = rarefied_richness(early_recs, n_subsample=100) if len(early_recs) >= 10 else s1
        s2_rare = rarefied_richness(late_recs, n_subsample=100) if len(late_recs) >= 10 else s2

        # Raw and rarefied changes
        raw_change = (s2 - s1) / s1 if s1 > 0 else np.nan
        rare_change = (s2_rare - s1_rare) / s1_rare if s1_rare > 0 else np.nan

        results.append({
            'cell_id': cell_id, 'lat': lat, 'lon': lon,
            's1': s1, 's2': s2, 'n_early': early_count, 'n_late': late_count,
            's1_rare': s1_rare, 's2_rare': s2_rare,
            'raw_change': raw_change, 'rare_change': rare_change
        })

        if s1 > 0:
            print(f"  S1={s1}, S2={s2}, raw={raw_change:+.3f}, rarefied={rare_change:+.3f}")
        else:
            print(f"  S1=0 (empty cell)")

        # Save cache periodically
        if cell_id % 10 == 0:
            _save_cache(results)

    _save_cache(results)

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------
    raw_changes = np.array([r['raw_change'] for r in results if not np.isnan(r['raw_change'])])
    rare_changes = np.array([r['rare_change'] for r in results if not np.isnan(r['rare_change'])])
    lats = np.array([r['lat'] for r in results if not np.isnan(r['raw_change'])])
    lons = np.array([r['lon'] for r in results if not np.isnan(r['raw_change'])])
    efforts_early = np.array([r['n_early'] for r in results if not np.isnan(r['raw_change'])])
    efforts_late = np.array([r['n_late'] for r in results if not np.isnan(r['raw_change'])])

    n = len(raw_changes)
    print(f"\n=== ANALYSIS (n={n} valid cells) ===\n")

    if n < 10:
        print("Too few valid cells. STOP.")
        return

    # --- Summary stats ---
    for label, data in [("Raw", raw_changes), ("Rarefied", rare_changes)]:
        print(f"\n{label} changes:")
        print(f"  Mean={np.mean(data):.4f}, SD={np.std(data, ddof=1):.4f}")
        print(f"  Skewness={stats.skew(data):.3f}")
        print(f"  Excess kurtosis={stats.kurtosis(data, fisher=True):.3f}")
        sw, sw_p = stats.shapiro(data) if len(data) >= 3 else (np.nan, np.nan)
        print(f"  Shapiro-Wilk: W={sw:.4f}, p={sw_p:.6f}")

    # --- Distribution fitting (on rarefied changes) ---
    data = rare_changes

    # Normal
    norm_p = stats.norm.fit(data)
    ll_norm = np.sum(stats.norm.logpdf(data, *norm_p))

    # Laplace
    lap_p = stats.laplace.fit(data)
    ll_lap = np.sum(stats.laplace.logpdf(data, *lap_p))

    # Subbotin
    sub_mu, sub_sigma, sub_beta, ll_sub = fit_subbotin(data)

    aic_norm = aic(2, ll_norm)
    aic_lap = aic(2, ll_lap)
    aic_sub = aic(3, ll_sub)

    print(f"\n--- Distribution Fits (rarefied, n={n}) ---")
    print(f"Normal:   mu={norm_p[0]:.4f}, sigma={norm_p[1]:.4f}, AIC={aic_norm:.2f}")
    print(f"Laplace:  loc={lap_p[0]:.4f}, scale={lap_p[1]:.4f}, AIC={aic_lap:.2f}")
    print(f"Subbotin: mu={sub_mu:.4f}, sigma={sub_sigma:.4f}, beta={sub_beta:.2f}, AIC={aic_sub:.2f}")
    print(f"  (beta=2 is Normal, beta=1 is Laplace, beta<1 is super-heavy tails)")

    best_model = min([("Normal", aic_norm), ("Laplace", aic_lap), ("Subbotin", aic_sub)], key=lambda x: x[1])
    print(f"Best model: {best_model[0]} (AIC={best_model[1]:.2f})")

    # --- Spatial predictors of extreme changes ---
    print(f"\n--- Spatial Predictors of Extreme Change ---")

    # Define tails: |change| > 90th percentile
    threshold = np.percentile(np.abs(data), 90)
    is_extreme = np.abs(data) > threshold

    print(f"Extreme threshold (90th pctile of |change|): {threshold:.3f}")
    print(f"Extreme cells: {np.sum(is_extreme)}/{n}")

    # Latitude correlation with |change|
    r_lat, p_lat = stats.pearsonr(lats, np.abs(data))
    print(f"Latitude vs |change|: r={r_lat:.3f}, p={p_lat:.4f}")

    # Longitude correlation
    r_lon, p_lon = stats.pearsonr(lons, np.abs(data))
    print(f"Longitude vs |change|: r={r_lon:.3f}, p={p_lon:.4f}")

    # Effort ratio as predictor
    effort_ratio = np.log10(efforts_late + 1) - np.log10(efforts_early + 1)
    r_eff, p_eff = stats.pearsonr(effort_ratio, data)
    print(f"log(effort ratio) vs change: r={r_eff:.3f}, p={p_eff:.4f}")

    # --- Figure ---
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))

    # (a) Histogram + fits
    ax = axes[0, 0]
    ax.hist(data, bins='auto', density=True, alpha=0.5, color='steelblue', edgecolor='white')
    x = np.linspace(data.min() - 0.5, data.max() + 0.5, 300)
    ax.plot(x, stats.norm.pdf(x, *norm_p), 'r-', lw=2, label='Normal')
    ax.plot(x, stats.laplace.pdf(x, *lap_p), 'g--', lw=2, label='Laplace')
    # Subbotin PDF
    sub_pdf = np.exp([subbotin_logpdf(xi, sub_mu, sub_sigma, sub_beta) for xi in x])
    ax.plot(x, sub_pdf, 'b:', lw=2, label=f'Subbotin (β={sub_beta:.1f})')
    ax.set_xlabel('Rarefied richness change')
    ax.set_ylabel('Density')
    ax.set_title(f'(a) Distribution of changes (n={n})')
    ax.legend(fontsize=7)
    annot = f"Kurtosis={stats.kurtosis(data, fisher=True):.2f}\nBest: {best_model[0]}"
    ax.annotate(annot, xy=(0.97, 0.97), xycoords='axes fraction', ha='right', va='top',
                fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

    # (b) Q-Q plot
    ax = axes[0, 1]
    stats.probplot(data, dist='norm', plot=ax)
    ax.set_title('(b) Q-Q plot vs Normal')
    ax.get_lines()[0].set(marker='o', markersize=4, color='steelblue')
    ax.get_lines()[1].set(color='red', linewidth=1.5)

    # (c) Raw vs Rarefied comparison
    ax = axes[0, 2]
    ax.scatter(raw_changes, rare_changes, alpha=0.6, s=20, c='steelblue')
    lims = [min(raw_changes.min(), rare_changes.min()) - 0.1,
            max(raw_changes.max(), rare_changes.max()) + 0.1]
    ax.plot(lims, lims, 'k--', alpha=0.3)
    r_rr, p_rr = stats.pearsonr(raw_changes, rare_changes)
    ax.set_xlabel('Raw change')
    ax.set_ylabel('Rarefied change')
    ax.set_title(f'(c) Raw vs Rarefied (r={r_rr:.2f})')

    # (d) Map of changes
    ax = axes[1, 0]
    scatter = ax.scatter(lons, lats, c=data, cmap='RdYlBu_r', s=40,
                         edgecolors='k', linewidths=0.5, vmin=-1, vmax=1)
    plt.colorbar(scatter, ax=ax, label='Rarefied change')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('(d) Spatial pattern of changes')

    # (e) Latitude vs |change|
    ax = axes[1, 1]
    ax.scatter(lats, np.abs(data), alpha=0.6, s=20, c='steelblue')
    ax.set_xlabel('Latitude')
    ax.set_ylabel('|Richness change|')
    ax.set_title(f'(e) Latitude vs |change| (r={r_lat:.2f}, p={p_lat:.3f})')
    # Add trend line
    z = np.polyfit(lats, np.abs(data), 1)
    ax.plot(sorted(lats), np.polyval(z, sorted(lats)), 'r-', alpha=0.5)

    # (f) Effort ratio vs change
    ax = axes[1, 2]
    ax.scatter(effort_ratio, data, alpha=0.6, s=20, c='steelblue')
    ax.set_xlabel('log10(effort late/early)')
    ax.set_ylabel('Rarefied change')
    ax.set_title(f'(f) Effort ratio vs change (r={r_eff:.2f}, p={p_eff:.3f})')
    z2 = np.polyfit(effort_ratio, data, 1)
    ax.plot(sorted(effort_ratio), np.polyval(z2, sorted(effort_ratio)), 'r-', alpha=0.5)

    plt.suptitle('Fat-Tailed Biodiversity Trends: European Birds (GBIF)', fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(FIG_PATH, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nFigure saved: {FIG_PATH}")

    # --- Save data CSV ---
    with open(DATA_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['cell_id','lat','lon','s1','s2',
                                                'n_early','n_late','s1_rare','s2_rare',
                                                'raw_change','rare_change'])
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    # --- Save results markdown ---
    kurt_raw = stats.kurtosis(raw_changes, fisher=True)
    kurt_rare = stats.kurtosis(rare_changes, fisher=True)
    sw_raw, swp_raw = stats.shapiro(raw_changes)
    sw_rare, swp_rare = stats.shapiro(rare_changes)

    md = f"""# Fat-Tailed Biodiversity Trends — Full Analysis Results

## Overview
- **Grid cells**: {len(GRID_CELLS)} queried, {n} valid
- **Taxon**: Aves (European birds)
- **Periods**: 2000-2010 vs 2014-2024
- **Method**: GBIF occurrence → species richness → proportional change
- **Effort correction**: Rarefaction to 100 records (20 iterations)

## Summary Statistics

| Metric | Raw | Rarefied |
|--------|-----|----------|
| Mean change | {np.mean(raw_changes):.4f} | {np.mean(rare_changes):.4f} |
| SD | {np.std(raw_changes, ddof=1):.4f} | {np.std(rare_changes, ddof=1):.4f} |
| Skewness | {stats.skew(raw_changes):.3f} | {stats.skew(rare_changes):.3f} |
| **Excess kurtosis** | **{kurt_raw:.3f}** | **{kurt_rare:.3f}** |
| Shapiro-Wilk W | {sw_raw:.4f} | {sw_rare:.4f} |
| Shapiro-Wilk p | {swp_raw:.6f} | {swp_rare:.6f} |

## Distribution Fits (rarefied changes)

| Distribution | Parameters | Log-lik | AIC | ΔAIC |
|-------------|-----------|---------|-----|------|
| Normal | μ={norm_p[0]:.4f}, σ={norm_p[1]:.4f} | {ll_norm:.2f} | {aic_norm:.2f} | {aic_norm - best_model[1]:.2f} |
| Laplace | loc={lap_p[0]:.4f}, scale={lap_p[1]:.4f} | {ll_lap:.2f} | {aic_lap:.2f} | {aic_lap - best_model[1]:.2f} |
| Subbotin | μ={sub_mu:.4f}, σ={sub_sigma:.4f}, β={sub_beta:.2f} | {ll_sub:.2f} | {aic_sub:.2f} | {aic_sub - best_model[1]:.2f} |

**Best model**: {best_model[0]} (AIC={best_model[1]:.2f})
- β={sub_beta:.2f} (2.0=Normal, 1.0=Laplace, <1.0=super-heavy tails)

## Spatial Predictors of Extreme Change

Extreme change threshold (90th percentile of |change|): {threshold:.3f}

| Predictor | r | p-value | Interpretation |
|-----------|---|---------|---------------|
| Latitude vs |change| | {r_lat:.3f} | {p_lat:.4f} | {'Higher latitudes show more extreme changes' if r_lat > 0 else 'Lower latitudes show more extreme changes'} |
| Longitude vs |change| | {r_lon:.3f} | {p_lon:.4f} | {'Eastern cells more extreme' if r_lon > 0 else 'Western cells more extreme'} |
| log(effort ratio) vs change | {r_eff:.3f} | {p_eff:.4f} | {'Effort increase → apparent richness increase' if r_eff > 0 else 'Effort-independent'} |
| Raw vs rarefied correlation | {r_rr:.3f} | {p_rr:.6f} | Rarefaction {'preserves' if r_rr > 0.8 else 'modifies'} the signal |

## Key Findings

1. **Kurtosis**: Excess kurtosis = {kurt_rare:.2f} ({'leptokurtic/fat-tailed' if kurt_rare > 0 else 'platykurtic/thin-tailed'})
2. **Best distribution**: {best_model[0]} (Subbotin β={sub_beta:.2f})
3. **Effort confound**: r={r_eff:.3f}, p={p_eff:.4f} — {'⚠️ effort drives signal' if p_eff < 0.05 else '✅ effort-independent'}
4. **Spatial pattern**: {'Latitude predicts extremes' if p_lat < 0.05 else 'No clear spatial pattern in extremes'}

## Output Files
- Script: `shared/scripts/fat_tailed_full.py`
- Data: `shared/data/fat_tailed_full_data.csv`
- Figure: `shared/data/fig_fat_tailed_full.png`
- Results: `shared/data/fat_tailed_full_results.md`

## Reference
McGill & Magurran (2026) Ecology Letters — biodiversity trends are leptokurtic
"""

    with open(RESULTS_PATH, 'w') as f:
        f.write(md)
    print(f"Results saved: {RESULTS_PATH}")


def _save_cache(results):
    """Save intermediate results to CSV cache."""
    with open(CACHE_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['cell_id','lat','lon','s1','s2',
                                                'n_early','n_late','s1_rare','s2_rare',
                                                'raw_change','rare_change'])
        writer.writeheader()
        for r in results:
            writer.writerow(r)


if __name__ == "__main__":
    main()
