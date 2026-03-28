#!/usr/bin/env Rscript
# ============================================================================
# RPP Uncertainty Propagation via Monte Carlo Simulation
# Paper 7 reviewer response: Is the ~48x correction robust?
# ============================================================================
#
# The paper claims Ambrosia is 40.6x more abundant than Zea in raw pollen,
# but after RPP correction this reverses to ~0.86x (Zea becomes more abundant).
# The reviewer questions whether this dramatic correction is robust.
#
# Approach: Monte Carlo simulation with uniform RPP distributions from
# published literature ranges.
# ============================================================================

set.seed(42)
n_iter <- 10000

# --- Raw pollen percentages from the paper ---
Zea_raw <- 0.208      # mean pollen % across ENA sites
Ambrosia_raw <- 8.44   # mean pollen % across ENA sites

# --- RPP ranges from literature (relative to Poaceae = 1) ---
# Ambrosia: high pollen producer, well-documented RPP range 4-8
# Zea: very poor pollen disperser (large heavy grains), RPP 0.05-0.20
Ambrosia_RPP <- runif(n_iter, 4.0, 8.0)
Zea_RPP <- runif(n_iter, 0.05, 0.20)

# --- Corrected vegetation representation ---
# Divide raw pollen % by RPP to estimate actual vegetation abundance
Zea_corrected <- Zea_raw / Zea_RPP
Ambrosia_corrected <- Ambrosia_raw / Ambrosia_RPP

# --- Ratio analysis ---
ratio_amb_to_zea <- Ambrosia_corrected / Zea_corrected

cat("=== RPP Uncertainty Propagation Results ===\n\n")

cat("Input values:\n")
cat(sprintf("  Zea raw pollen: %.3f%%\n", Zea_raw))
cat(sprintf("  Ambrosia raw pollen: %.2f%%\n", Ambrosia_raw))
cat(sprintf("  Raw ratio (Ambrosia/Zea): %.1fx\n\n", Ambrosia_raw / Zea_raw))

cat("RPP distributions:\n")
cat("  Ambrosia RPP: Uniform(4.0, 8.0)\n")
cat("  Zea RPP: Uniform(0.05, 0.20)\n\n")

cat("Corrected Zea vegetation representation:\n")
cat(sprintf("  Median: %.3f%%\n", median(Zea_corrected)))
cat(sprintf("  95%% CI: %.3f%% - %.3f%%\n",
            quantile(Zea_corrected, 0.025), quantile(Zea_corrected, 0.975)))

cat("\nCorrected Ambrosia vegetation representation:\n")
cat(sprintf("  Median: %.3f%%\n", median(Ambrosia_corrected)))
cat(sprintf("  95%% CI: %.3f%% - %.3f%%\n",
            quantile(Ambrosia_corrected, 0.025), quantile(Ambrosia_corrected, 0.975)))

cat("\nCorrected ratio (Ambrosia/Zea):\n")
cat(sprintf("  Median: %.2fx\n", median(ratio_amb_to_zea)))
cat(sprintf("  Mean: %.2fx\n", mean(ratio_amb_to_zea)))
cat(sprintf("  95%% CI: %.2fx - %.2fx\n",
            quantile(ratio_amb_to_zea, 0.025), quantile(ratio_amb_to_zea, 0.975)))

cat(sprintf("\nPr(corrected Zea > corrected Ambrosia): %.1f%%\n",
            100 * mean(Zea_corrected > Ambrosia_corrected)))

# --- Sensitivity analysis ---
cat("\n=== Sensitivity: What RPPs would preserve the 'Ambrosia dominates' story? ===\n")
cat("For corrected ratio > 5:\n")
cat(sprintf("  Zea_RPP/Ambrosia_RPP must exceed %.4f\n", 5 * Zea_raw / Ambrosia_raw))
cat(sprintf("  Pr(ratio > 5): %.1f%%\n", 100 * mean(ratio_amb_to_zea > 5)))
cat(sprintf("  Pr(ratio > 10): %.1f%%\n", 100 * mean(ratio_amb_to_zea > 10)))
cat(sprintf("  Pr(ratio > 20): %.1f%%\n", 100 * mean(ratio_amb_to_zea > 20)))

cat("\nRequired Zea RPP for Ambrosia to remain 5x larger:\n")
target <- 5 * Zea_raw / Ambrosia_raw
cat(sprintf("  If Ambrosia RPP = 4.0: Zea RPP must be %.3f (vs max published 0.20)\n", target * 4.0))
cat(sprintf("  If Ambrosia RPP = 6.0: Zea RPP must be %.3f\n", target * 6.0))
cat(sprintf("  If Ambrosia RPP = 8.0: Zea RPP must be %.3f\n", target * 8.0))
cat("  All values far exceed published Zea RPP range (0.05-0.20)\n")
