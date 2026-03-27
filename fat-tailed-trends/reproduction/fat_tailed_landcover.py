#!/usr/bin/env python3
"""
fat_tailed_landcover.py
=======================
Adds MODIS land-cover change proxies as predictors of extreme biodiversity
changes in the 72 reliable European bird grid cells.

Approach
--------
1. Load fat_tailed_full_data.csv, filter to reliable cells (S1>=30, n_early>=50).
2. Compute effort growth rate  = log(n_late / n_early).
3. Assign categorical land-cover type based on lat/lon (CORINE-style proxy).
4. Query GBIF occurrence counts for a set of known urban-adapted bird species
   per cell to compute an urbanisation index.
5. Multiple regression:  |rare_change| ~ lat + lon + landcover + effort_growth + S1
6. ANOVA: |rare_change| ~ landcover category
7. 4-panel figure saved to fig_fat_tailed_landcover.png.
8. Results summary saved to fat_tailed_landcover_results.md.
"""

import sys, os, time, json, math, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests

warnings.filterwarnings("ignore")

DATA_DIR = "/home/ayu/ecolab2/shared/data"
INPUT_CSV = os.path.join(DATA_DIR, "fat_tailed_full_data.csv")
OUT_FIG = os.path.join(DATA_DIR, "fig_fat_tailed_landcover.png")
OUT_MD = os.path.join(DATA_DIR, "fat_tailed_landcover_results.md")

# Urban-adapted bird species (European) with GBIF taxon keys
URBAN_SPECIES = {
    "Passer domesticus": 5231190,
    "Columba livia": 2495414,
    "Sturnus vulgaris": 9515886,
    "Pica pica": 2482520,
    "Corvus corone": 2482082,
    "Turdus merula": 9703694,
    "Streptopelia decaocto": 2495347,
    "Apus apus": 5228513,
    "Delichon urbicum": 9598250,
    "Hirundo rustica": 9598129,
}

SESSION = requests.Session()


def gbif_count(lat, lon, taxon_key=None):
    """Return GBIF occurrence count within +/-1 degree of lat/lon."""
    params = {
        "decimalLatitude": f"{lat-1},{lat+1}",
        "decimalLongitude": f"{lon-1},{lon+1}",
        "hasCoordinate": "true",
        "limit": 0,
    }
    if taxon_key:
        params["taxonKey"] = taxon_key
    url = "https://api.gbif.org/v1/occurrence/search"
    for attempt in range(3):
        try:
            r = SESSION.get(url, params=params, timeout=30)
            r.raise_for_status()
            return r.json().get("count", 0)
        except Exception as e:
            if attempt == 2:
                return np.nan
            time.sleep(1)


def query_cell(lat, lon):
    """Query GBIF for total bird count and urban species counts for one cell."""
    total = gbif_count(lat, lon, taxon_key=212)  # class Aves
    urban = 0
    for sp_name, tk in URBAN_SPECIES.items():
        c = gbif_count(lat, lon, taxon_key=tk)
        if c is not None and not (isinstance(c, float) and np.isnan(c)):
            urban += c
        time.sleep(0.02)
    return lat, lon, total, urban


def assign_landcover(lat, lon):
    """Assign CORINE-style land-cover category based on lat/lon."""
    if lat < 43:
        return "Mediterranean"
    elif lat <= 52 and lon < 5:
        return "W-Europe mixed"
    elif lat <= 52 and lon <= 15:
        return "C-Europe agri-forest"
    else:
        return "N/E-Europe forest"


def main():
    print("=" * 60)
    print("Fat-tailed bird analysis: Land-cover predictor")
    print("=" * 60)
    sys.stdout.flush()

    # 1. Load data
    df = pd.read_csv(INPUT_CSV)
    print(f"Loaded {len(df)} cells from {INPUT_CSV}")

    # Filter to reliable cells
    df = df[(df["s1"] >= 30) & (df["n_early"] >= 50)].copy()
    df = df.reset_index(drop=True)
    print(f"Reliable cells (S1>=30, n_early>=50): {len(df)}")

    # 2. Derived variables
    df["effort_growth"] = np.log(df["n_late"] / df["n_early"])
    df["abs_rare_change"] = df["rare_change"].abs()
    df["landcover"] = df.apply(lambda r: assign_landcover(r["lat"], r["lon"]), axis=1)

    print("\nLand-cover category counts:")
    print(df["landcover"].value_counts().to_string())
    sys.stdout.flush()

    # 3. Query GBIF for urban species index per cell (parallel)
    print(f"\nQuerying GBIF for {len(df)} cells (parallel, 4 workers)...")
    sys.stdout.flush()

    results = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {}
        for i, row in df.iterrows():
            f = pool.submit(query_cell, row["lat"], row["lon"])
            futures[f] = i
        done_count = 0
        for f in as_completed(futures):
            done_count += 1
            lat, lon, total, urban = f.result()
            idx = futures[f]
            results[idx] = (total, urban)
            frac = urban / total if total and total > 0 else 0
            print(f"  [{done_count}/{len(df)}] Cell ({lat:.0f},{lon:.0f}): total={total:,}, urban={urban:,}, frac={frac:.3f}")
            sys.stdout.flush()

    # Assign results back
    df["gbif_total"] = df.index.map(lambda i: results[i][0])
    df["gbif_urban"] = df.index.map(lambda i: results[i][1])
    df["urban_fraction"] = df["gbif_urban"] / df["gbif_total"]
    df["urban_fraction"] = df["urban_fraction"].replace([np.inf, -np.inf], np.nan)

    # Save enriched data
    enriched_path = os.path.join(DATA_DIR, "fat_tailed_landcover_enriched.csv")
    df.to_csv(enriched_path, index=False)
    print(f"\nEnriched data saved to {enriched_path}")
    sys.stdout.flush()

    # ── Statistical analyses ───────────────────────────────────────────
    print("\n" + "=" * 60)
    print("STATISTICAL ANALYSES")
    print("=" * 60)
    sys.stdout.flush()

    # ANOVA: |rare_change| ~ landcover
    print("\n--- ANOVA: |rare_change| ~ landcover ---")
    groups = [g["abs_rare_change"].dropna().values for _, g in df.groupby("landcover")]
    group_labels = [name for name, _ in df.groupby("landcover")]
    if len(groups) >= 2:
        f_stat, p_anova = stats.f_oneway(*groups)
        print(f"F = {f_stat:.3f}, p = {p_anova:.4f}")
        h_stat, p_kw = stats.kruskal(*groups)
        print(f"Kruskal-Wallis H = {h_stat:.3f}, p = {p_kw:.4f}")
    else:
        f_stat, p_anova, h_stat, p_kw = np.nan, np.nan, np.nan, np.nan

    # Per-group summary
    print("\nPer land-cover group summary (|rare_change|):")
    group_summary = df.groupby("landcover")["abs_rare_change"].agg(["count", "mean", "median", "std"])
    print(group_summary.to_string())
    sys.stdout.flush()

    # Multiple regression
    print("\n--- Multiple regression ---")
    df_reg = df.dropna(subset=["abs_rare_change", "effort_growth", "urban_fraction", "lat", "lon"]).copy()
    df_reg = pd.get_dummies(df_reg, columns=["landcover"], drop_first=True, dtype=float)

    lc_cols = [c for c in df_reg.columns if c.startswith("landcover_")]
    pred_cols = ["lat", "lon", "effort_growth", "s1", "urban_fraction"] + lc_cols

    X = df_reg[pred_cols].copy()
    X = sm.add_constant(X)
    y = df_reg["abs_rare_change"]

    if len(df_reg) > len(pred_cols) + 2:
        model = sm.OLS(y, X).fit()
        print(model.summary())
        reg_summary = model.summary().as_text()
        r2 = model.rsquared
        adj_r2 = model.rsquared_adj
        f_pvalue = model.f_pvalue
    else:
        print("Not enough data for full regression")
        reg_summary = "Insufficient data for full model"
        r2 = adj_r2 = f_pvalue = np.nan
        model = None

    # Simpler regression: effort_growth + urban_fraction
    print("\n--- Reduced regression (effort_growth + urban_fraction) ---")
    X_simple = df_reg[["effort_growth", "urban_fraction"]].copy()
    X_simple = sm.add_constant(X_simple)
    y_simple = df_reg["abs_rare_change"]
    if len(df_reg) > 5:
        model_simple = sm.OLS(y_simple, X_simple).fit()
        print(model_simple.summary())
        simple_summary = model_simple.summary().as_text()
    else:
        simple_summary = "Insufficient data"
        model_simple = None

    # Correlation matrix
    print("\n--- Correlation matrix ---")
    corr_cols = ["abs_rare_change", "effort_growth", "urban_fraction", "lat", "lon", "s1"]
    corr_df = df[corr_cols].dropna()
    corr_matrix = corr_df.corr()
    print(corr_matrix.round(3).to_string())
    sys.stdout.flush()

    # Bivariate correlations for reporting
    valid = df.dropna(subset=["effort_growth", "abs_rare_change"])
    r_eff, p_eff = stats.pearsonr(valid["effort_growth"], valid["abs_rare_change"])
    valid_u = df.dropna(subset=["urban_fraction", "abs_rare_change"])
    if len(valid_u) > 2:
        r_urb, p_urb = stats.pearsonr(valid_u["urban_fraction"], valid_u["abs_rare_change"])
    else:
        r_urb, p_urb = np.nan, np.nan

    # ── Figure ─────────────────────────────────────────────────────────
    print("\nGenerating 4-panel figure...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel A: |rare_change| by land-cover category (box plot)
    ax = axes[0, 0]
    lc_order = ["Mediterranean", "W-Europe mixed", "C-Europe agri-forest", "N/E-Europe forest"]
    lc_data = [df[df["landcover"] == lc]["abs_rare_change"].dropna().values for lc in lc_order]
    lc_labels_short = ["Medit.", "W-Eur\nmixed", "C-Eur\nagri-for.", "N/E-Eur\nforest"]
    bp = ax.boxplot(lc_data, labels=lc_labels_short, patch_artist=True,
                    boxprops=dict(facecolor="lightblue", alpha=0.7))
    colors_box = ["#e6550d", "#fdae6b", "#74c476", "#31a354"]
    for patch, color in zip(bp["boxes"], colors_box):
        patch.set_facecolor(color)
    ax.set_ylabel("|Rarefied richness change|")
    ax.set_title(f"A) Land-cover type\n(ANOVA F={f_stat:.2f}, p={p_anova:.3f})")
    ax.axhline(y=df["abs_rare_change"].median(), color="grey", ls="--", alpha=0.5, label="overall median")

    # Panel B: |rare_change| vs effort growth
    ax = axes[0, 1]
    ax.scatter(valid["effort_growth"], valid["abs_rare_change"], alpha=0.6, c="steelblue", edgecolors="k", s=40)
    if len(valid) > 2:
        z = np.polyfit(valid["effort_growth"], valid["abs_rare_change"], 1)
        xline = np.linspace(valid["effort_growth"].min(), valid["effort_growth"].max(), 50)
        ax.plot(xline, np.polyval(z, xline), "r-", lw=2)
    ax.set_xlabel("Effort growth rate [log(n_late/n_early)]")
    ax.set_ylabel("|Rarefied richness change|")
    ax.set_title(f"B) Effort growth\n(r={r_eff:.3f}, p={p_eff:.3f})")

    # Panel C: |rare_change| vs urban fraction
    ax = axes[1, 0]
    ax.scatter(valid_u["urban_fraction"], valid_u["abs_rare_change"], alpha=0.6, c="darkorange", edgecolors="k", s=40)
    if len(valid_u) > 2 and not np.isnan(r_urb):
        z = np.polyfit(valid_u["urban_fraction"], valid_u["abs_rare_change"], 1)
        xline = np.linspace(valid_u["urban_fraction"].min(), valid_u["urban_fraction"].max(), 50)
        ax.plot(xline, np.polyval(z, xline), "r-", lw=2)
    ax.set_xlabel("Urban species fraction")
    ax.set_ylabel("|Rarefied richness change|")
    ax.set_title(f"C) Urbanisation index\n(r={r_urb:.3f}, p={p_urb:.3f})")

    # Panel D: Geographic map coloured by |rare_change| with landcover outlines
    ax = axes[1, 1]
    sc = ax.scatter(df["lon"], df["lat"], c=df["abs_rare_change"],
                    cmap="YlOrRd", s=60, edgecolors="k", linewidths=0.5, zorder=5)
    plt.colorbar(sc, ax=ax, label="|Rarefied richness change|")
    ax.axhline(43, color="blue", ls="--", alpha=0.4, label="lat=43")
    ax.axhline(52, color="blue", ls="--", alpha=0.4, label="lat=52")
    ax.axvline(5, color="green", ls="--", alpha=0.4, label="lon=5")
    ax.axvline(15, color="green", ls="--", alpha=0.4, label="lon=15")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("D) Geographic distribution\n(dashed = land-cover boundaries)")
    ax.legend(fontsize=7, loc="lower left")

    fig.suptitle("Land-cover predictors of extreme bird richness change", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_FIG, dpi=150, bbox_inches="tight")
    print(f"Figure saved to {OUT_FIG}")

    # ── Results markdown ───────────────────────────────────────────────
    md = f"""# Land-cover predictors of extreme bird richness change

## Data
- Input: {INPUT_CSV}
- Reliable cells: {len(df)} (S1>=30, n_early>=50)
- Land-cover categories assigned by lat/lon (CORINE proxy)
- Urban species index: fraction of GBIF records from {len(URBAN_SPECIES)} known urban-adapted species

## Land-cover category distribution
{df['landcover'].value_counts().to_string()}

## 1. ANOVA: |rare_change| ~ land-cover type
- F-statistic: {f_stat:.3f}
- p-value: {p_anova:.4f}
- Kruskal-Wallis H: {h_stat:.3f}, p = {p_kw:.4f}

### Per-group summary (|rarefied richness change|)
{group_summary.to_string()}

**Interpretation**: {"Land-cover category significantly predicts extreme richness change magnitude." if p_anova < 0.05 else "Land-cover category does NOT significantly predict extreme richness change magnitude (p > 0.05)."}

## 2. Bivariate correlations with |rare_change|
| Predictor | r | p-value |
|-----------|---|---------|
| Effort growth | {r_eff:.3f} | {p_eff:.4f} |
| Urban fraction | {r_urb:.3f} | {p_urb:.4f} |

## 3. Multiple regression: |rare_change| ~ lat + lon + effort + S1 + urban_frac + landcover
```
{reg_summary}
```

## 4. Reduced regression: |rare_change| ~ effort_growth + urban_fraction
```
{simple_summary}
```

## 5. Correlation matrix
```
{corr_matrix.round(3).to_string()}
```

## Key findings
- Effort growth rate (citizen-science expansion) correlation with |rare_change|: r = {r_eff:.3f} (p = {p_eff:.4f})
- Urban species fraction correlation with |rare_change|: r = {r_urb:.3f} (p = {p_urb:.4f})
- ANOVA by land-cover type: F = {f_stat:.2f}, p = {p_anova:.4f}
- Full model R-squared: {r2:.3f} (adj. {adj_r2:.3f})

## Figure
![Land-cover predictors](fig_fat_tailed_landcover.png)

## Methods note
Land-cover categories are coarse lat/lon proxies (not pixel-level MODIS), suitable
for testing broad biogeographic patterns but not fine-scale land-use effects.
The urban species index uses GBIF occurrence counts for {len(URBAN_SPECIES)} known
European synanthropic species relative to total bird occurrences.
"""
    with open(OUT_MD, "w") as f:
        f.write(md)
    print(f"Results saved to {OUT_MD}")
    print("\nDone.")


if __name__ == "__main__":
    main()
