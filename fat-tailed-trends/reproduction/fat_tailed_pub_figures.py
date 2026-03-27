#!/usr/bin/env python3
"""
Generate publication-quality figures for the fat-tailed biodiversity trends manuscript.

Figure 1: Distribution fits (2x2: birds/plants × histogram/Q-Q)
Figure 2: Cross-taxon comparison (kurtosis, asymmetry, AIC)
Figure 3: Spatial predictors (map, urban fraction, S1)
"""

import csv
import os
import numpy as np
from scipy import stats
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

DATA_DIR = "/home/ayu/ecolab2/shared/data"
OUT_DIR = "/home/ayu/ecolab2/reports/figures"
os.makedirs(OUT_DIR, exist_ok=True)

# --- Load data ---
def load_filtered(path, min_s1=30, min_n_early=50):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                rc = float(r['rare_change'])
                s1 = float(r['s1'])
                ne = float(r['n_early'])
                if not np.isnan(rc) and s1 >= min_s1 and ne >= min_n_early:
                    rows.append(r)
            except:
                pass
    return rows

bird_rows = load_filtered(os.path.join(DATA_DIR, "fat_tailed_full_data.csv"))
plant_rows = load_filtered(os.path.join(DATA_DIR, "fat_tailed_plants_cache.csv"))

bird_data = np.array([float(r['rare_change']) for r in bird_rows])
plant_data = np.array([float(r['rare_change']) for r in plant_rows])
bird_lats = np.array([float(r['lat']) for r in bird_rows])
bird_lons = np.array([float(r['lon']) for r in bird_rows])

# --- Asymmetric Laplace ---
def al_logpdf(x, mu, b1, b2):
    x = np.asarray(x)
    return np.where(x < mu, -np.log(b1+b2) - np.abs(x-mu)/b1, -np.log(b1+b2) - np.abs(x-mu)/b2)

def fit_al(data):
    def neg_ll(p):
        mu, lb1, lb2 = p
        return -np.sum(al_logpdf(data, mu, np.exp(lb1), np.exp(lb2)))
    med = np.median(data)
    mad = max(np.median(np.abs(data - med)), 0.01)
    res = minimize(neg_ll, [med, np.log(mad), np.log(mad)], method='Nelder-Mead', options={'maxiter':5000})
    return res.x[0], np.exp(res.x[1]), np.exp(res.x[2])

# Fit distributions
bird_norm = stats.norm.fit(bird_data)
bird_lap = stats.laplace.fit(bird_data)
bird_t = stats.t.fit(bird_data)
bird_al_mu, bird_al_b1, bird_al_b2 = fit_al(bird_data)

plant_norm = stats.norm.fit(plant_data)
plant_lap = stats.laplace.fit(plant_data)
plant_t = stats.t.fit(plant_data)
plant_al_mu, plant_al_b1, plant_al_b2 = fit_al(plant_data)

# Color scheme
C_BIRD = '#2171b5'
C_PLANT = '#238b45'
C_NORMAL = '#e41a1c'
C_LAPLACE = '#ff7f00'
C_STUDENTT = '#984ea3'
C_ASYMLAP = '#000000'

# ======================================================================
# FIGURE 1: Distribution fits (2×2)
# ======================================================================
fig1, axes = plt.subplots(2, 2, figsize=(10, 8))

for idx, (data, label, color, norm_p, lap_p, t_p, al_mu, al_b1, al_b2) in enumerate([
    (bird_data, 'Birds (Aves)', C_BIRD, bird_norm, bird_lap, bird_t, bird_al_mu, bird_al_b1, bird_al_b2),
    (plant_data, 'Plants (Tracheophyta)', C_PLANT, plant_norm, plant_lap, plant_t, plant_al_mu, plant_al_b1, plant_al_b2),
]):
    # Histogram + fits
    ax = axes[idx, 0]
    ax.hist(data, bins='auto', density=True, alpha=0.4, color=color, edgecolor='white', zorder=2)
    x = np.linspace(data.min()-0.3, min(data.max()+0.3, 3.0), 300)
    ax.plot(x, stats.norm.pdf(x, *norm_p), color=C_NORMAL, lw=1.5, ls='-', label='Normal', zorder=3)
    ax.plot(x, stats.laplace.pdf(x, *lap_p), color=C_LAPLACE, lw=1.5, ls='--', label='Laplace', zorder=3)
    ax.plot(x, stats.t.pdf(x, *t_p), color=C_STUDENTT, lw=1.5, ls=':', label=f'Student-t', zorder=3)
    al_pdf = np.exp(al_logpdf(x, al_mu, al_b1, al_b2))
    ax.plot(x, al_pdf, color=C_ASYMLAP, lw=2.5, ls='-', label='Asym. Laplace', zorder=4)
    ax.set_xlabel('Proportional richness change', fontsize=9)
    ax.set_ylabel('Density', fontsize=9)
    kurt = stats.kurtosis(data, fisher=True)
    n = len(data)
    panel = 'a' if idx == 0 else 'c'
    ax.set_title(f'({panel}) {label} (n={n}, κ={kurt:.1f})', fontsize=10, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_xlim(data.min()-0.3, min(data.max()+0.3, 3.0))

    # Q-Q plot
    ax = axes[idx, 1]
    sorted_d = np.sort(data)
    theoretical = stats.norm.ppf((np.arange(1, n+1) - 0.5)/n)
    ax.scatter(theoretical, sorted_d, s=15, c=color, alpha=0.7, edgecolors='k', linewidths=0.3, zorder=3)
    lims = [min(theoretical.min(), sorted_d.min())-0.2, max(theoretical.max(), sorted_d.max())+0.2]
    ax.plot(lims, lims, 'r-', lw=1.5, alpha=0.5, zorder=2)
    ax.set_xlabel('Normal theoretical quantiles', fontsize=9)
    ax.set_ylabel('Sample quantiles', fontsize=9)
    panel = 'b' if idx == 0 else 'd'
    ax.set_title(f'({panel}) Q-Q plot vs Normal', fontsize=10, fontweight='bold')

plt.tight_layout()
fig1.savefig(os.path.join(OUT_DIR, "fig1_distribution_fits.png"), dpi=300, bbox_inches='tight')
fig1.savefig(os.path.join(OUT_DIR, "fig1_distribution_fits.pdf"), bbox_inches='tight')
print("Figure 1 saved")

# ======================================================================
# FIGURE 2: Cross-taxon comparison (3 panels)
# ======================================================================
fig2, axes = plt.subplots(1, 3, figsize=(12, 4))

# (a) Kurtosis comparison
ax = axes[0]
taxa = ['Birds', 'Plants']
kurts = [stats.kurtosis(bird_data, fisher=True), stats.kurtosis(plant_data, fisher=True)]
bars = ax.bar(taxa, kurts, color=[C_BIRD, C_PLANT], edgecolor='k', linewidth=0.8, width=0.5)
ax.axhline(0, color='grey', ls='--', alpha=0.5)
ax.set_ylabel('Excess kurtosis', fontsize=10)
ax.set_title('(a) Leptokurtosis', fontsize=11, fontweight='bold')
for bar, k in zip(bars, kurts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, f'{k:.1f}',
            ha='center', va='bottom', fontsize=10, fontweight='bold')

# (b) Asymmetry ratio
ax = axes[1]
ratios = [bird_al_b2/bird_al_b1, plant_al_b2/plant_al_b1]
bars = ax.bar(taxa, ratios, color=[C_BIRD, C_PLANT], edgecolor='k', linewidth=0.8, width=0.5)
ax.axhline(1, color='grey', ls='--', alpha=0.5, label='Symmetric')
ax.set_ylabel('Asymmetry ratio (b₂/b₁)', fontsize=10)
ax.set_title('(b) Right-tail heaviness', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)
for bar, r in zip(bars, ratios):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, f'{r:.1f}×',
            ha='center', va='bottom', fontsize=10, fontweight='bold')

# (c) AIC comparison
ax = axes[2]
# Bird AICs
b_ll_n = np.sum(stats.norm.logpdf(bird_data, *bird_norm))
b_ll_l = np.sum(stats.laplace.logpdf(bird_data, *bird_lap))
b_ll_t = np.sum(stats.t.logpdf(bird_data, *bird_t))
b_ll_al = np.sum(al_logpdf(bird_data, bird_al_mu, bird_al_b1, bird_al_b2))
b_aics = [4-2*b_ll_n, 4-2*b_ll_l, 6-2*b_ll_t, 6-2*b_ll_al]
# Plant AICs
p_ll_n = np.sum(stats.norm.logpdf(plant_data, *plant_norm))
p_ll_l = np.sum(stats.laplace.logpdf(plant_data, *plant_lap))
p_ll_t = np.sum(stats.t.logpdf(plant_data, *plant_t))
p_ll_al = np.sum(al_logpdf(plant_data, plant_al_mu, plant_al_b1, plant_al_b2))
p_aics = [4-2*p_ll_n, 4-2*p_ll_l, 6-2*p_ll_t, 6-2*p_ll_al]

models = ['Normal', 'Laplace', 'Student-t', 'Asym.\nLaplace']
x_pos = np.arange(len(models))
w = 0.35
ax.bar(x_pos - w/2, [a - min(b_aics) for a in b_aics], w, color=C_BIRD, edgecolor='k', linewidth=0.5, label='Birds')
ax.bar(x_pos + w/2, [a - min(p_aics) for a in p_aics], w, color=C_PLANT, edgecolor='k', linewidth=0.5, label='Plants')
ax.set_xticks(x_pos)
ax.set_xticklabels(models, fontsize=9)
ax.set_ylabel('ΔAIC (from best)', fontsize=10)
ax.set_title('(c) Model comparison', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)

plt.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "fig2_cross_taxon.png"), dpi=300, bbox_inches='tight')
fig2.savefig(os.path.join(OUT_DIR, "fig2_cross_taxon.pdf"), bbox_inches='tight')
print("Figure 2 saved")

# ======================================================================
# FIGURE 3: Spatial predictors (3 panels)
# ======================================================================
# Load enriched landcover data
lc_path = os.path.join(DATA_DIR, "fat_tailed_landcover_enriched.csv")
import pandas as pd
lc_df = pd.read_csv(lc_path)

fig3, axes = plt.subplots(1, 3, figsize=(13, 4.5))

# (a) Geographic map
ax = axes[0]
sc = ax.scatter(lc_df['lon'], lc_df['lat'], c=lc_df['abs_rare_change'],
                cmap='YlOrRd', s=50, edgecolors='k', linewidths=0.5, vmin=0, vmax=0.8)
plt.colorbar(sc, ax=ax, label='|Richness change|', shrink=0.8)
ax.set_xlabel('Longitude', fontsize=9)
ax.set_ylabel('Latitude', fontsize=9)
ax.set_title('(a) Geographic distribution', fontsize=11, fontweight='bold')

# (b) Urban fraction vs |change|
ax = axes[1]
valid = lc_df.dropna(subset=['urban_fraction', 'abs_rare_change'])
ax.scatter(valid['urban_fraction'], valid['abs_rare_change'], s=30, alpha=0.6,
           c=C_BIRD, edgecolors='k', linewidths=0.3)
r_u, p_u = stats.pearsonr(valid['urban_fraction'], valid['abs_rare_change'])
z = np.polyfit(valid['urban_fraction'], valid['abs_rare_change'], 1)
xline = np.linspace(valid['urban_fraction'].min(), valid['urban_fraction'].max(), 50)
ax.plot(xline, np.polyval(z, xline), 'r-', lw=2)
ax.set_xlabel('Urban species fraction', fontsize=9)
ax.set_ylabel('|Richness change|', fontsize=9)
ax.set_title(f'(b) Urbanisation (r={r_u:.2f}, p={p_u:.3f})', fontsize=11, fontweight='bold')

# (c) S1 vs |change|
ax = axes[2]
ax.scatter(lc_df['s1'], lc_df['abs_rare_change'], s=30, alpha=0.6,
           c=C_BIRD, edgecolors='k', linewidths=0.3)
r_s, p_s = stats.pearsonr(lc_df['s1'], lc_df['abs_rare_change'])
z = np.polyfit(lc_df['s1'], lc_df['abs_rare_change'], 1)
xline = np.linspace(lc_df['s1'].min(), lc_df['s1'].max(), 50)
ax.plot(xline, np.polyval(z, xline), 'r-', lw=2)
ax.set_xlabel('Baseline species richness (S₁)', fontsize=9)
ax.set_ylabel('|Richness change|', fontsize=9)
ax.set_title(f'(c) Baseline richness (r={r_s:.2f}, p={p_s:.3f})', fontsize=11, fontweight='bold')

plt.tight_layout()
fig3.savefig(os.path.join(OUT_DIR, "fig3_spatial_predictors.png"), dpi=300, bbox_inches='tight')
fig3.savefig(os.path.join(OUT_DIR, "fig3_spatial_predictors.pdf"), bbox_inches='tight')
print("Figure 3 saved")

print(f"\nAll figures saved to {OUT_DIR}")
