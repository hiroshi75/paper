#!/usr/bin/env python3
"""
Diagnostic Framework Simulation (Package A)
============================================
Evaluates the performance of placebo/negative-control diagnostics
for detecting observer-effort bias in citizen-science ecological data.

Model:
  X_i ~ N(0,1)                          # environmental covariate
  B_i = rho_B * X_i + eta_i             # true biological signal
  E_i = rho_E * X_i + zeta_i            # observer effort
  Y_main_i = alpha + beta_B*B_i + beta_E*E_i + eps_i   # main taxon count
  Y_diag_i = alpha + beta_E*E_i + eps_diag_i           # diagnostic taxon (no biology)

The diagnostic taxon lacks beta_B — any signal in Y_diag ~ X is purely effort artefact.
"""

import itertools
import os
import time
import warnings

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.api import OLS, add_constant

warnings.filterwarnings("ignore")

# ── Configuration ──────────────────────────────────────────────────────────
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

PARAM_GRID = {
    "rho_B": [0.0, 0.2, 0.5, 1.0],
    "rho_E": [0.0, 0.3, 0.6, 0.9],
    "beta_B": [0.0, 0.3, 0.5, 1.0],
    "beta_E": [0.5, 1.0, 2.0],
    "n_cells": [500, 5000],
    "noise_sd": [1.0, 2.0],
}

N_ITER = 200  # iterations per scenario
ALPHA_INTERCEPT = 0.0
SIGNIFICANCE = 0.05
RATIO_THRESHOLD = 0.5  # diagnostic/main ratio threshold for screening

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})


def run_single_iteration(rng, rho_B, rho_E, beta_B, beta_E, n_cells, noise_sd):
    """Run one Monte-Carlo iteration and return summary statistics."""
    X = rng.standard_normal(n_cells)
    eta = rng.standard_normal(n_cells)
    zeta = rng.standard_normal(n_cells)

    B = rho_B * X + eta
    E = rho_E * X + zeta

    eps_main = rng.standard_normal(n_cells) * noise_sd
    eps_diag = rng.standard_normal(n_cells) * noise_sd

    Y_main = ALPHA_INTERCEPT + beta_B * B + beta_E * E + eps_main
    Y_diag = ALPHA_INTERCEPT + beta_E * E + eps_diag

    # Naive regressions: Y ~ X
    X_design = add_constant(X)

    res_main = OLS(Y_main, X_design).fit()
    res_diag = OLS(Y_diag, X_design).fit()

    coef_main = res_main.params[1]
    p_main = res_main.pvalues[1]
    coef_diag = res_diag.params[1]
    p_diag = res_diag.pvalues[1]

    abs_main = np.abs(coef_main)
    abs_diag = np.abs(coef_diag)
    ratio = abs_diag / abs_main if abs_main > 1e-10 else np.nan

    return {
        "coef_main": coef_main,
        "p_main": p_main,
        "coef_diag": coef_diag,
        "p_diag": p_diag,
        "ratio": ratio,
        "main_sig": int(p_main < SIGNIFICANCE),
        "diag_sig": int(p_diag < SIGNIFICANCE),
    }


def run_scenario(params):
    """Run all iterations for one parameter combination."""
    rho_B, rho_E, beta_B, beta_E, n_cells, noise_sd = params
    rng = np.random.default_rng(seed=42)

    results = []
    for _ in range(N_ITER):
        res = run_single_iteration(rng, rho_B, rho_E, beta_B, beta_E, n_cells, noise_sd)
        results.append(res)

    ratios = np.array([r["ratio"] for r in results])
    ratios_clean = ratios[np.isfinite(ratios)]

    main_sig = np.array([r["main_sig"] for r in results])
    diag_sig = np.array([r["diag_sig"] for r in results])
    ratios_arr = ratios.copy()

    # Rates
    n_total = len(results)
    fpr = float(main_sig.mean()) if beta_B == 0 else np.nan
    tpr = float(main_sig.mean()) if beta_B > 0 else np.nan
    diag_detection = float(diag_sig.mean()) if rho_E > 0 else np.nan

    # Corrected FPR: main significant AND ratio < threshold
    if beta_B == 0:
        corrected_sig = np.array([
            r["main_sig"] and (r["ratio"] < RATIO_THRESHOLD if np.isfinite(r["ratio"]) else False)
            for r in results
        ])
        corrected_fpr = float(corrected_sig.mean())
    else:
        corrected_fpr = np.nan

    return {
        "rho_B": rho_B,
        "rho_E": rho_E,
        "beta_B": beta_B,
        "beta_E": beta_E,
        "n_cells": n_cells,
        "noise_sd": noise_sd,
        "fpr": fpr,
        "tpr": tpr,
        "diag_detection_rate": diag_detection,
        "ratio_mean": float(np.nanmean(ratios_clean)) if len(ratios_clean) > 0 else np.nan,
        "ratio_median": float(np.nanmedian(ratios_clean)) if len(ratios_clean) > 0 else np.nan,
        "ratio_p5": float(np.nanpercentile(ratios_clean, 5)) if len(ratios_clean) > 0 else np.nan,
        "ratio_p95": float(np.nanpercentile(ratios_clean, 95)) if len(ratios_clean) > 0 else np.nan,
        "corrected_fpr": corrected_fpr,
        "main_sig_rate": float(main_sig.mean()),
        "diag_sig_rate": float(diag_sig.mean()),
    }


def build_param_combos():
    """Generate all parameter combinations."""
    keys = list(PARAM_GRID.keys())
    vals = [PARAM_GRID[k] for k in keys]
    combos = list(itertools.product(*vals))
    return combos


def run_all_scenarios():
    """Run the full simulation grid."""
    combos = build_param_combos()
    total = len(combos)
    print(f"Total scenarios: {total}")
    print(f"Iterations per scenario: {N_ITER}")
    print(f"Estimated total regressions: {total * N_ITER * 2:,}")

    results = []
    t0 = time.time()
    for i, combo in enumerate(combos):
        if (i + 1) % 50 == 0 or i == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (total - i - 1) / rate if rate > 0 else 0
            print(f"  Scenario {i+1}/{total}  ({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")
        results.append(run_scenario(combo))

    elapsed = time.time() - t0
    print(f"Simulation complete in {elapsed:.1f}s")
    return pd.DataFrame(results)


# ── Plotting ───────────────────────────────────────────────────────────────

def fig_fpr_heatmap(df):
    """
    Figure 1: Heatmap of false positive rate (naive vs diagnostic-corrected)
    across rho_E x beta_E conditions.
    """
    # Filter: beta_B == 0 (false positive conditions), default n_cells and noise
    sub = df[(df["beta_B"] == 0) & (df["n_cells"] == 500) & (df["noise_sd"] == 1.0)].copy()
    if sub.empty:
        print("WARNING: No data for FPR heatmap")
        return

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.8), sharey=True)

    for ax, col, title in zip(
        axes,
        ["fpr", "corrected_fpr"],
        ["Naive analysis", "With diagnostic screening"],
    ):
        pivot = sub.pivot_table(index="rho_E", columns="beta_E", values=col, aggfunc="mean")
        pivot = pivot.sort_index(ascending=False)

        sns.heatmap(
            pivot,
            ax=ax,
            annot=True,
            fmt=".2f",
            cmap="RdYlGn_r",
            vmin=0,
            vmax=1.0,
            linewidths=0.5,
            linecolor="white",
            cbar_kws={"label": "False positive rate", "shrink": 0.8},
            annot_kws={"size": 8},
        )
        ax.set_title(title, fontweight="bold", pad=8)
        ax.set_xlabel(r"Effort contamination ($\beta_E$)")
        ax.set_ylabel(r"Effort–environment confounding ($\rho_E$)")

    # Add a dashed line at nominal alpha
    fig.suptitle("", y=1.02)
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, "fig_simulation_fpr.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved: {path}")


def fig_ratio_distribution(df):
    """
    Figure 2: Distribution of diagnostic/main coefficient ratio across conditions.
    """
    # Use a subset of interesting conditions
    sub = df[
        (df["n_cells"] == 500) & (df["noise_sd"] == 1.0) & (df["beta_B"].isin([0.0, 0.5]))
    ].copy()

    sub["label"] = sub.apply(
        lambda r: f"ρ_E={r.rho_E:.1f}, β_E={r.beta_E:.1f}", axis=1
    )
    sub["beta_B_label"] = sub["beta_B"].map({0.0: "No signal (β_B=0)", 0.5: "True signal (β_B=0.5)"})

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.0), sharey=False)

    for ax, (bb, grp) in zip(axes, sub.groupby("beta_B_label")):
        # Plot ratio_median vs rho_E, colored by beta_E
        for beta_e, g2 in grp.groupby("beta_E"):
            g2s = g2.sort_values("rho_E")
            ax.plot(g2s["rho_E"], g2s["ratio_median"], "o-", markersize=4,
                    label=f"β_E={beta_e:.1f}")
            ax.fill_between(g2s["rho_E"], g2s["ratio_p5"], g2s["ratio_p95"], alpha=0.15)

        ax.axhline(RATIO_THRESHOLD, color="red", ls="--", lw=0.8, label="Threshold (0.5)")
        ax.set_title(bb, fontweight="bold", pad=6)
        ax.set_xlabel(r"Effort–environment confounding ($\rho_E$)")
        ax.set_ylabel("Diagnostic / Main ratio")
        ax.legend(loc="best", frameon=True, framealpha=0.9)
        ax.set_ylim(bottom=0)

    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, "fig_simulation_ratio_dist.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved: {path}")


def fig_roc_curve(df):
    """
    Figure 3: ROC-style curve — true positive rate vs false positive rate
    with and without diagnostic screening.
    """
    # We need paired FPR/TPR across matching conditions
    # Strategy: for each (rho_E, beta_E, n_cells, noise_sd),
    # FPR comes from beta_B=0, TPR from beta_B>0
    # We vary rho_B to get different TPR levels

    n_cells_val = 500
    noise_val = 1.0

    base = df[(df["n_cells"] == n_cells_val) & (df["noise_sd"] == noise_val)].copy()

    fpr_df = base[base["beta_B"] == 0][
        ["rho_E", "beta_E", "rho_B", "fpr", "corrected_fpr"]
    ].copy()
    tpr_df = base[base["beta_B"] > 0][
        ["rho_E", "beta_E", "rho_B", "beta_B", "tpr", "main_sig_rate", "ratio_median"]
    ].copy()

    # For the "corrected TPR", a true positive passes if main is significant AND ratio < threshold
    # We need to re-derive this. Use main_sig_rate as naive TPR.
    # For corrected TPR: when there IS a real signal but also effort bias,
    # the ratio should be < 1 (diagnostic smaller than main), so true positives survive.
    # We approximate: corrected_tpr ≈ tpr * P(ratio < threshold)
    # But more accurately, we should compute it in the simulation. Let's add it.

    # Instead, re-run a focused simulation for ROC data
    rng = np.random.default_rng(42)
    roc_data = []

    rho_E_vals = [0.0, 0.3, 0.6, 0.9]
    beta_E_vals = [0.5, 1.0, 2.0]
    beta_B_vals = [0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5]
    rho_B_val = 0.5  # fixed moderate biological correlation

    for rho_E in rho_E_vals:
        for beta_E in beta_E_vals:
            for beta_B in beta_B_vals:
                naive_sig = 0
                corrected_sig = 0
                total = 0
                for _ in range(N_ITER):
                    res = run_single_iteration(
                        rng, rho_B_val, rho_E, beta_B, beta_E,
                        n_cells_val, noise_val
                    )
                    total += 1
                    if res["main_sig"]:
                        naive_sig += 1
                        if np.isfinite(res["ratio"]) and res["ratio"] < RATIO_THRESHOLD:
                            corrected_sig += 1

                roc_data.append({
                    "rho_E": rho_E,
                    "beta_E": beta_E,
                    "beta_B": beta_B,
                    "naive_rate": naive_sig / total,
                    "corrected_rate": corrected_sig / total,
                })

    roc_df = pd.DataFrame(roc_data)

    # Plot: for each (rho_E, beta_E) condition, plot the curve of
    # (FPR at beta_B=0, TPR at beta_B>0) as beta_B varies
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.8), sharey=True)

    for ax, beta_E in zip(axes, beta_E_vals):
        for rho_E in rho_E_vals:
            sub = roc_df[(roc_df["beta_E"] == beta_E) & (roc_df["rho_E"] == rho_E)].sort_values("beta_B")

            # FPR is the rate at beta_B=0
            fpr_naive = sub.loc[sub["beta_B"] == 0, "naive_rate"].values[0]
            fpr_corr = sub.loc[sub["beta_B"] == 0, "corrected_rate"].values[0]

            tpr_naive = sub.loc[sub["beta_B"] > 0, "naive_rate"].values
            tpr_corr = sub.loc[sub["beta_B"] > 0, "corrected_rate"].values
            bb_vals = sub.loc[sub["beta_B"] > 0, "beta_B"].values

            # Plot naive as dashed, corrected as solid
            ax.plot(
                [fpr_naive] * len(tpr_naive), tpr_naive,
                "x", markersize=4, alpha=0.5, color=f"C{rho_E_vals.index(rho_E)}"
            )
            ax.plot(
                [fpr_corr] * len(tpr_corr), tpr_corr,
                "o", markersize=3, alpha=0.7, color=f"C{rho_E_vals.index(rho_E)}",
                label=f"ρ_E={rho_E}" if beta_E == beta_E_vals[0] else ""
            )

            # Connect with arrows from naive to corrected
            for tn, tc, fn, fc in zip(tpr_naive, tpr_corr, [fpr_naive]*len(tpr_naive), [fpr_corr]*len(tpr_corr)):
                ax.annotate(
                    "", xy=(fc, tc), xytext=(fn, tn),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=0.4, alpha=0.3),
                )

        ax.plot([0, 1], [0, 1], "k--", lw=0.5, alpha=0.3)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(f"β_E = {beta_E}", fontweight="bold", pad=6)
        ax.set_xlabel("False positive rate")
        if beta_E == beta_E_vals[0]:
            ax.set_ylabel("True positive rate")
        ax.set_aspect("equal")

    # Build legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="x", color="gray", linestyle="None", markersize=5, label="Naive"),
        Line2D([0], [0], marker="o", color="gray", linestyle="None", markersize=4, label="With diagnostic"),
    ] + [
        Line2D([0], [0], marker="o", color=f"C{i}", linestyle="None", markersize=4, label=f"ρ_E={v}")
        for i, v in enumerate(rho_E_vals)
    ]
    axes[-1].legend(handles=legend_elements, loc="lower right", frameon=True, framealpha=0.9, fontsize=7)

    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, "fig_simulation_roc.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved: {path}")


# ── Main ───────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("Diagnostic Framework Simulation — Package A")
    print("=" * 60)

    df = run_all_scenarios()

    # Save CSV
    csv_path = os.path.join(OUTPUT_DIR, "simulation_results.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nSaved results: {csv_path}")
    print(f"Shape: {df.shape}")

    # Summary statistics
    print("\n── Key Results ──")
    fpr_rows = df[df["beta_B"] == 0].dropna(subset=["fpr"])
    if not fpr_rows.empty:
        print(f"Naive FPR range: {fpr_rows['fpr'].min():.3f} – {fpr_rows['fpr'].max():.3f}")
        print(f"Corrected FPR range: {fpr_rows['corrected_fpr'].min():.3f} – {fpr_rows['corrected_fpr'].max():.3f}")

    tpr_rows = df[df["beta_B"] > 0].dropna(subset=["tpr"])
    if not tpr_rows.empty:
        print(f"TPR range: {tpr_rows['tpr'].min():.3f} – {tpr_rows['tpr'].max():.3f}")

    diag_rows = df[df["rho_E"] > 0].dropna(subset=["diag_detection_rate"])
    if not diag_rows.empty:
        print(f"Diagnostic detection rate range: {diag_rows['diag_detection_rate'].min():.3f} – {diag_rows['diag_detection_rate'].max():.3f}")

    # Generate figures
    print("\n── Generating Figures ──")
    fig_fpr_heatmap(df)
    fig_ratio_distribution(df)
    fig_roc_curve(df)

    print("\nDone.")
