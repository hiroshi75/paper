#!/usr/bin/env python3
"""
Fat-Tailed Biodiversity Trends — European Plants (Tracheophyta)

Replicates the bird fat-tailed analysis (fat_tailed_full.py) for vascular plants.
Uses the same 100 European grid cells, same two periods, same rarefaction method.

Fits Normal, Laplace, Student-t, and Asymmetric Laplace distributions to
assess whether fat-tailed biodiversity trends are taxon-general.
"""

import json
import os
import time
import urllib.request
import urllib.error
import numpy as np
from scipy import stats
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

GBIF_BASE = "https://api.gbif.org/v1/occurrence/search"
TAXON_KEY = 7707728  # Tracheophyta (vascular plants)
LIMIT = 300

PERIOD_EARLY = ("2000-01-01", "2010-12-31")
PERIOD_LATE = ("2014-01-01", "2024-12-31")

# Same 100 grid cells as bird analysis
GRID_CELLS = []
for lat in np.arange(38, 62, 2.5):
    for lon in np.arange(-10, 30, 4):
        GRID_CELLS.append((float(lat), float(lon)))

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
CACHE_PATH = os.path.join(OUTPUT_DIR, "fat_tailed_plants_cache.csv")
FIG_PATH = os.path.join(OUTPUT_DIR, "fig_fat_tailed_plants.png")
RESULTS_PATH = os.path.join(OUTPUT_DIR, "fat_tailed_plants_results.md")

# Filtering thresholds
MIN_S1 = 30          # minimum early-period species richness
MIN_N_EARLY = 50     # minimum early-period record count


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def query_gbif(lat_min, lat_max, lon_min, lon_max, event_start, event_end):
    """Query GBIF occurrence API. Returns list of records and total count."""
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


# ---------------------------------------------------------------------------
# Distribution fitting
# ---------------------------------------------------------------------------

def fit_normal(data):
    """Fit Normal distribution. Returns (params, loglik, aic)."""
    mu, sigma = stats.norm.fit(data)
    ll = np.sum(stats.norm.logpdf(data, mu, sigma))
    return {'mu': mu, 'sigma': sigma}, ll, 2 * 2 - 2 * ll


def fit_laplace(data):
    """Fit Laplace distribution. Returns (params, loglik, aic)."""
    loc, scale = stats.laplace.fit(data)
    ll = np.sum(stats.laplace.logpdf(data, loc, scale))
    return {'loc': loc, 'scale': scale}, ll, 2 * 2 - 2 * ll


def fit_student_t(data):
    """Fit Student-t distribution. Returns (params, loglik, aic)."""
    df, loc, scale = stats.t.fit(data)
    ll = np.sum(stats.t.logpdf(data, df, loc, scale))
    return {'df': df, 'loc': loc, 'scale': scale}, ll, 2 * 3 - 2 * ll


def asymmetric_laplace_logpdf(x, mu, b1, b2):
    """Log-PDF of Asymmetric Laplace distribution.

    f(x) = 1/(b1+b2) * exp(-|x-mu|/b1) if x < mu
    f(x) = 1/(b1+b2) * exp(-|x-mu|/b2) if x >= mu
    """
    x = np.asarray(x)
    logpdf = np.where(
        x < mu,
        -np.log(b1 + b2) - np.abs(x - mu) / b1,
        -np.log(b1 + b2) - np.abs(x - mu) / b2
    )
    return logpdf


def fit_asymmetric_laplace(data):
    """Fit Asymmetric Laplace via MLE (Nelder-Mead). Returns (params, loglik, aic)."""
    data = np.asarray(data)

    def neg_ll(params):
        mu, log_b1, log_b2 = params
        b1 = np.exp(log_b1)
        b2 = np.exp(log_b2)
        ll = np.sum(asymmetric_laplace_logpdf(data, mu, b1, b2))
        return -ll

    # Initial guess: median for mu, MAD for scale
    mu0 = np.median(data)
    mad = np.median(np.abs(data - mu0))
    if mad < 1e-6:
        mad = np.std(data)

    x0 = [mu0, np.log(mad + 1e-6), np.log(mad + 1e-6)]
    res = minimize(neg_ll, x0, method='Nelder-Mead',
                   options={'maxiter': 5000, 'xatol': 1e-8, 'fatol': 1e-8})

    mu = res.x[0]
    b1 = np.exp(res.x[1])
    b2 = np.exp(res.x[2])
    ll = -res.fun
    k = 3  # mu, b1, b2
    return {'mu': mu, 'b1': b1, 'b2': b2}, ll, 2 * k - 2 * ll


def asymmetric_laplace_pdf(x, mu, b1, b2):
    """PDF of Asymmetric Laplace distribution."""
    return np.exp(asymmetric_laplace_logpdf(x, mu, b1, b2))


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

def _save_cache(results):
    """Save intermediate results to CSV cache."""
    with open(CACHE_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['cell_id','lat','lon','s1','s2',
                                                'n_early','n_late','s1_rare','s2_rare',
                                                'raw_change','rare_change'])
        writer.writeheader()
        for r in results:
            writer.writerow(r)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    np.random.seed(42)

    print(f"Fat-Tailed Biodiversity Trends — European Plants (Tracheophyta)")
    print(f"Grid cells: {len(GRID_CELLS)}")
    print(f"Periods: {PERIOD_EARLY} vs {PERIOD_LATE}")
    print(f"Filtering: S1>={MIN_S1}, n_early>={MIN_N_EARLY}")
    print()

    # Check for cache
    results = []
    if os.path.exists(CACHE_PATH):
        print(f"Loading cache from {CACHE_PATH}")
        with open(CACHE_PATH) as f:
            reader = csv.DictReader(f)
            for row in reader:
                results.append({k: float(v) if k != 'cell_id' else int(float(v))
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
            print(f"  S1={s1}, S2={s2}, n_early={early_count}, raw={raw_change:+.3f}, rarefied={rare_change:+.3f}")
        else:
            print(f"  S1=0 (empty cell)")

        # Save cache periodically
        if cell_id % 10 == 0:
            _save_cache(results)

    _save_cache(results)

    # ------------------------------------------------------------------
    # Filter cells: S1 >= MIN_S1 and n_early >= MIN_N_EARLY
    # ------------------------------------------------------------------
    filtered = [r for r in results
                if not np.isnan(r['raw_change'])
                and r['s1'] >= MIN_S1
                and r['n_early'] >= MIN_N_EARLY]

    raw_changes = np.array([r['raw_change'] for r in filtered])
    rare_changes = np.array([r['rare_change'] for r in filtered])
    lats = np.array([r['lat'] for r in filtered])
    lons = np.array([r['lon'] for r in filtered])
    efforts_early = np.array([r['n_early'] for r in filtered])
    efforts_late = np.array([r['n_late'] for r in filtered])

    n_total = len(results)
    n = len(filtered)
    print(f"\n=== ANALYSIS ===")
    print(f"Total cells queried: {n_total}")
    print(f"Valid cells after filtering (S1>={MIN_S1}, n_early>={MIN_N_EARLY}): {n}")

    if n < 10:
        print("Too few valid cells after filtering. STOP.")
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

    norm_params, ll_norm, aic_norm = fit_normal(data)
    lap_params, ll_lap, aic_lap = fit_laplace(data)
    t_params, ll_t, aic_t = fit_student_t(data)
    al_params, ll_al, aic_al = fit_asymmetric_laplace(data)

    models = [
        ("Normal", norm_params, ll_norm, aic_norm),
        ("Laplace", lap_params, ll_lap, aic_lap),
        ("Student-t", t_params, ll_t, aic_t),
        ("Asym. Laplace", al_params, ll_al, aic_al),
    ]

    best_model_name, _, _, best_aic = min(models, key=lambda x: x[3])

    print(f"\n--- Distribution Fits (rarefied, n={n}) ---")
    print(f"Normal:         mu={norm_params['mu']:.4f}, sigma={norm_params['sigma']:.4f}, "
          f"LL={ll_norm:.2f}, AIC={aic_norm:.2f}")
    print(f"Laplace:        loc={lap_params['loc']:.4f}, scale={lap_params['scale']:.4f}, "
          f"LL={ll_lap:.2f}, AIC={aic_lap:.2f}")
    print(f"Student-t:      df={t_params['df']:.2f}, loc={t_params['loc']:.4f}, "
          f"scale={t_params['scale']:.4f}, LL={ll_t:.2f}, AIC={aic_t:.2f}")
    print(f"Asym. Laplace:  mu={al_params['mu']:.4f}, b1={al_params['b1']:.4f}, "
          f"b2={al_params['b2']:.4f}, LL={ll_al:.2f}, AIC={aic_al:.2f}")
    print(f"\nBest model: {best_model_name} (AIC={best_aic:.2f})")

    # --- Spatial predictors ---
    threshold = np.percentile(np.abs(data), 90)
    is_extreme = np.abs(data) > threshold
    r_lat, p_lat = stats.pearsonr(lats, np.abs(data))
    r_lon, p_lon = stats.pearsonr(lons, np.abs(data))
    effort_ratio = np.log10(efforts_late + 1) - np.log10(efforts_early + 1)
    r_eff, p_eff = stats.pearsonr(effort_ratio, data)
    r_rr, p_rr = stats.pearsonr(raw_changes, rare_changes)

    print(f"\n--- Spatial Predictors ---")
    print(f"Extreme threshold (90th pctile |change|): {threshold:.3f}")
    print(f"Extreme cells: {np.sum(is_extreme)}/{n}")
    print(f"Latitude vs |change|: r={r_lat:.3f}, p={p_lat:.4f}")
    print(f"Longitude vs |change|: r={r_lon:.3f}, p={p_lon:.4f}")
    print(f"log(effort ratio) vs change: r={r_eff:.3f}, p={p_eff:.4f}")

    # --- Figure ---
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))

    # (a) Histogram + fits
    ax = axes[0, 0]
    ax.hist(data, bins='auto', density=True, alpha=0.5, color='forestgreen', edgecolor='white')
    x = np.linspace(data.min() - 0.5, data.max() + 0.5, 300)
    ax.plot(x, stats.norm.pdf(x, norm_params['mu'], norm_params['sigma']),
            'r-', lw=2, label='Normal')
    ax.plot(x, stats.laplace.pdf(x, lap_params['loc'], lap_params['scale']),
            'g--', lw=2, label='Laplace')
    ax.plot(x, stats.t.pdf(x, t_params['df'], t_params['loc'], t_params['scale']),
            'b:', lw=2, label=f"Student-t (df={t_params['df']:.1f})")
    ax.plot(x, asymmetric_laplace_pdf(x, al_params['mu'], al_params['b1'], al_params['b2']),
            'm-.', lw=2, label='Asym. Laplace')
    ax.set_xlabel('Rarefied richness change')
    ax.set_ylabel('Density')
    ax.set_title(f'(a) Distribution of changes (n={n})')
    ax.legend(fontsize=7)
    kurt_val = stats.kurtosis(data, fisher=True)
    skew_val = stats.skew(data)
    annot = f"Kurtosis={kurt_val:.2f}\nSkewness={skew_val:.2f}\nBest: {best_model_name}"
    ax.annotate(annot, xy=(0.97, 0.97), xycoords='axes fraction', ha='right', va='top',
                fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

    # (b) Q-Q plot
    ax = axes[0, 1]
    stats.probplot(data, dist='norm', plot=ax)
    ax.set_title('(b) Q-Q plot vs Normal')
    ax.get_lines()[0].set(marker='o', markersize=4, color='forestgreen')
    ax.get_lines()[1].set(color='red', linewidth=1.5)

    # (c) Raw vs Rarefied comparison
    ax = axes[0, 2]
    ax.scatter(raw_changes, rare_changes, alpha=0.6, s=20, c='forestgreen')
    lims = [min(raw_changes.min(), rare_changes.min()) - 0.1,
            max(raw_changes.max(), rare_changes.max()) + 0.1]
    ax.plot(lims, lims, 'k--', alpha=0.3)
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
    ax.scatter(lats, np.abs(data), alpha=0.6, s=20, c='forestgreen')
    ax.set_xlabel('Latitude')
    ax.set_ylabel('|Richness change|')
    ax.set_title(f'(e) Latitude vs |change| (r={r_lat:.2f}, p={p_lat:.3f})')
    z = np.polyfit(lats, np.abs(data), 1)
    ax.plot(sorted(lats), np.polyval(z, sorted(lats)), 'r-', alpha=0.5)

    # (f) Effort ratio vs change
    ax = axes[1, 2]
    ax.scatter(effort_ratio, data, alpha=0.6, s=20, c='forestgreen')
    ax.set_xlabel('log10(effort late/early)')
    ax.set_ylabel('Rarefied change')
    ax.set_title(f'(f) Effort ratio vs change (r={r_eff:.2f}, p={p_eff:.3f})')
    z2 = np.polyfit(effort_ratio, data, 1)
    ax.plot(sorted(effort_ratio), np.polyval(z2, sorted(effort_ratio)), 'r-', alpha=0.5)

    plt.suptitle('Fat-Tailed Biodiversity Trends: European Plants — Tracheophyta (GBIF)',
                 fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(FIG_PATH, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nFigure saved: {FIG_PATH}")

    # --- Compare with birds ---
    # Load bird results if available
    bird_cache = os.path.join(OUTPUT_DIR, "fat_tailed_full_cache.csv")
    bird_comparison = ""
    if os.path.exists(bird_cache):
        bird_results = []
        with open(bird_cache) as f:
            reader = csv.DictReader(f)
            for row in reader:
                bird_results.append({k: float(v) for k, v in row.items()})
        bird_rare = np.array([r['rare_change'] for r in bird_results
                              if not np.isnan(r['rare_change'])
                              and r['s1'] >= MIN_S1
                              and r['n_early'] >= MIN_N_EARLY])
        if len(bird_rare) > 10:
            bird_kurt = stats.kurtosis(bird_rare, fisher=True)
            bird_skew = stats.skew(bird_rare)
            plant_kurt = stats.kurtosis(data, fisher=True)
            plant_skew = stats.skew(data)

            print(f"\n=== CROSS-TAXON COMPARISON ===")
            print(f"Birds (Aves):        kurtosis={bird_kurt:.3f}, skewness={bird_skew:.3f}, n={len(bird_rare)}")
            print(f"Plants (Tracheophyta): kurtosis={plant_kurt:.3f}, skewness={plant_skew:.3f}, n={n}")

            if bird_kurt > 0 and plant_kurt > 0:
                print("=> Both taxa show leptokurtic (fat-tailed) distributions — pattern is taxon-general")
            elif bird_kurt > 0 and plant_kurt <= 0:
                print("=> Birds are leptokurtic but plants are not — pattern may be taxon-specific")
            elif bird_kurt <= 0 and plant_kurt > 0:
                print("=> Plants are leptokurtic but birds are not — unexpected reversal")
            else:
                print("=> Neither taxon shows clear leptokurtosis")

            bird_comparison = f"""
## Cross-Taxon Comparison (Birds vs Plants)

| Metric | Birds (Aves) | Plants (Tracheophyta) |
|--------|-------------|----------------------|
| n (filtered cells) | {len(bird_rare)} | {n} |
| Excess kurtosis | {bird_kurt:.3f} | {plant_kurt:.3f} |
| Skewness | {bird_skew:.3f} | {plant_skew:.3f} |

**Interpretation**: {"Both taxa show leptokurtic (fat-tailed) distributions, supporting the hypothesis that fat-tailed biodiversity trends are taxon-general." if bird_kurt > 0 and plant_kurt > 0 else "The pattern differs between taxa — fat-tailed trends may be taxon-specific." if bird_kurt > 0 and plant_kurt <= 0 else "Results are mixed — further investigation needed."}
"""
    else:
        print("\n[Note] Bird cache not found — skipping cross-taxon comparison")
        bird_comparison = "\n## Cross-Taxon Comparison\nBird data cache not available for comparison.\n"

    # --- Save results markdown ---
    kurt_raw = stats.kurtosis(raw_changes, fisher=True)
    kurt_rare = stats.kurtosis(rare_changes, fisher=True)
    skew_raw = stats.skew(raw_changes)
    skew_rare = stats.skew(rare_changes)
    sw_raw, swp_raw = stats.shapiro(raw_changes)
    sw_rare, swp_rare = stats.shapiro(rare_changes)

    md = f"""# Fat-Tailed Biodiversity Trends — European Plants (Tracheophyta)

## Overview
- **Grid cells**: {len(GRID_CELLS)} queried, {n} valid (after filtering S1>={MIN_S1}, n_early>={MIN_N_EARLY})
- **Taxon**: Tracheophyta (vascular plants), taxonKey={TAXON_KEY}
- **Periods**: 2000-2010 vs 2014-2024
- **Method**: GBIF occurrence -> species richness -> proportional change
- **Effort correction**: Rarefaction to 100 records (20 iterations)

## Summary Statistics

| Metric | Raw | Rarefied |
|--------|-----|----------|
| Mean change | {np.mean(raw_changes):.4f} | {np.mean(rare_changes):.4f} |
| SD | {np.std(raw_changes, ddof=1):.4f} | {np.std(rare_changes, ddof=1):.4f} |
| Skewness | {skew_raw:.3f} | {skew_rare:.3f} |
| **Excess kurtosis** | **{kurt_raw:.3f}** | **{kurt_rare:.3f}** |
| Shapiro-Wilk W | {sw_raw:.4f} | {sw_rare:.4f} |
| Shapiro-Wilk p | {swp_raw:.6f} | {swp_rare:.6f} |

## Distribution Fits (rarefied changes)

| Distribution | Parameters | Log-lik | AIC | dAIC |
|-------------|-----------|---------|-----|------|
| Normal | mu={norm_params['mu']:.4f}, sigma={norm_params['sigma']:.4f} | {ll_norm:.2f} | {aic_norm:.2f} | {aic_norm - best_aic:.2f} |
| Laplace | loc={lap_params['loc']:.4f}, scale={lap_params['scale']:.4f} | {ll_lap:.2f} | {aic_lap:.2f} | {aic_lap - best_aic:.2f} |
| Student-t | df={t_params['df']:.2f}, loc={t_params['loc']:.4f}, scale={t_params['scale']:.4f} | {ll_t:.2f} | {aic_t:.2f} | {aic_t - best_aic:.2f} |
| Asym. Laplace | mu={al_params['mu']:.4f}, b1={al_params['b1']:.4f}, b2={al_params['b2']:.4f} | {ll_al:.2f} | {aic_al:.2f} | {aic_al - best_aic:.2f} |

**Best model**: {best_model_name} (AIC={best_aic:.2f})

## Spatial Predictors of Extreme Change

Extreme change threshold (90th percentile of |change|): {threshold:.3f}

| Predictor | r | p-value | Interpretation |
|-----------|---|---------|---------------|
| Latitude vs |change| | {r_lat:.3f} | {p_lat:.4f} | {'Higher latitudes show more extreme changes' if r_lat > 0 else 'Lower latitudes show more extreme changes'} |
| Longitude vs |change| | {r_lon:.3f} | {p_lon:.4f} | {'Eastern cells more extreme' if r_lon > 0 else 'Western cells more extreme'} |
| log(effort ratio) vs change | {r_eff:.3f} | {p_eff:.4f} | {'Effort increase -> apparent richness increase' if r_eff > 0 else 'Effort-independent'} |
| Raw vs rarefied correlation | {r_rr:.3f} | {p_rr:.6f} | Rarefaction {'preserves' if r_rr > 0.8 else 'modifies'} the signal |

## Key Findings

1. **Kurtosis**: Excess kurtosis = {kurt_rare:.2f} ({'leptokurtic/fat-tailed' if kurt_rare > 0 else 'platykurtic/thin-tailed'})
2. **Skewness**: {skew_rare:.2f} ({'positive skew — more extreme gains' if skew_rare > 0.3 else 'negative skew — more extreme losses' if skew_rare < -0.3 else 'approximately symmetric'})
3. **Best distribution**: {best_model_name}
4. **Effort confound**: r={r_eff:.3f}, p={p_eff:.4f} — {'effort drives signal' if p_eff < 0.05 else 'effort-independent'}
5. **Spatial pattern**: {'Latitude predicts extremes' if p_lat < 0.05 else 'No clear spatial pattern in extremes'}
{bird_comparison}
## Output Files
- Script: `shared/scripts/fat_tailed_plants.py`
- Cache: `shared/data/fat_tailed_plants_cache.csv`
- Figure: `shared/data/fig_fat_tailed_plants.png`
- Results: `shared/data/fat_tailed_plants_results.md`

## Reference
McGill & Magurran (2026) Ecology Letters — biodiversity trends are leptokurtic
"""

    with open(RESULTS_PATH, 'w') as f:
        f.write(md)
    print(f"Results saved: {RESULTS_PATH}")


if __name__ == "__main__":
    main()
