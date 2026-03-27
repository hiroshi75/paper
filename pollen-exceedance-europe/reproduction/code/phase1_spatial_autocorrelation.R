#!/usr/bin/env Rscript
# Phase 1: Spatial Autocorrelation Analysis of Pollen Signal Onset Ages
# Tests Moran's I across 3 European regions to assess site independence
# Output: phase1_spatial_autocorrelation_results.md

library(spdep)
library(sf)

cat("=== Phase 1: Spatial Autocorrelation of Signal Onset Ages ===\n\n")

# ---- 1. Load and combine data ----

uk <- read.csv("/home/ayu/archeco/shared/pilot_signal_phase_uk.csv",
               stringsAsFactors = FALSE)
uk$region <- "UK"

scan <- read.csv("/home/ayu/archeco/shared/signal_phase_scandinavia.csv",
                 stringsAsFactors = FALSE)
scan$region <- "Scandinavia"

ce <- read.csv("/home/ayu/archeco/shared/signal_phase_central_europe.csv",
               stringsAsFactors = FALSE)
ce$region <- "Central_Europe"

# Standardize columns
common_cols <- c("datasetid", "sitename", "lat", "lon", "signal_onset_age", "region")
uk_sub <- uk[, common_cols]
scan_sub <- scan[, common_cols]
ce_sub <- ce[, common_cols]

combined <- rbind(uk_sub, scan_sub, ce_sub)
cat(sprintf("Raw combined: %d rows\n", nrow(combined)))

# Remove NAs in coordinates or signal_onset_age
combined <- combined[!is.na(combined$lat) & !is.na(combined$lon) &
                     !is.na(combined$signal_onset_age), ]
cat(sprintf("After NA removal: %d rows\n", nrow(combined)))

# Remove exact duplicate coordinates (keep first)
combined$coord_key <- paste(combined$lat, combined$lon)
dups <- duplicated(combined$coord_key)
cat(sprintf("Duplicate coordinates: %d\n", sum(dups)))
combined <- combined[!dups, ]
combined$coord_key <- NULL
cat(sprintf("After dedup: %d rows\n\n", nrow(combined)))

for (r in unique(combined$region)) {
  cat(sprintf("  %s: %d sites\n", r, sum(combined$region == r)))
}
cat("\n")

# ---- 2. Function: Moran's I analysis ----

run_moran <- function(data, label, threshold_km) {
  cat(sprintf("--- %s (n=%d, threshold=%dkm) ---\n", label, nrow(data), threshold_km))

  if (nrow(data) < 5) {
    cat("  Too few sites, skipping.\n\n")
    return(NULL)
  }

  coords <- cbind(data$lon, data$lat)

  # Distance-based neighbors
  # dnearneigh uses great-circle distance when longlat=TRUE, threshold in km
  nb <- dnearneigh(coords, d1 = 0, d2 = threshold_km, longlat = TRUE)

  # Check connectivity
  n_no_neighbors <- sum(card(nb) == 0)
  cat(sprintf("  Sites with no neighbors within %dkm: %d / %d\n",
              threshold_km, n_no_neighbors, nrow(data)))

  if (n_no_neighbors == nrow(data)) {
    cat("  No neighbors found at this threshold. Skipping.\n\n")
    return(NULL)
  }

  # For sites with no neighbors, we still include them (they contribute 0 to Moran's I)
  # But moran.test requires connected graph - use zero.policy=TRUE
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

  # Moran's I test
  mt <- moran.test(data$signal_onset_age, lw, zero.policy = TRUE)

  cat(sprintf("  Moran's I:   %.4f\n", mt$estimate["Moran I statistic"]))
  cat(sprintf("  Expected I:  %.4f\n", mt$estimate["Expectation"]))
  cat(sprintf("  Variance:    %.6f\n", mt$estimate["Variance"]))
  cat(sprintf("  Z-score:     %.4f\n", mt$statistic))
  cat(sprintf("  p-value:     %.6f\n", mt$p.value))
  cat(sprintf("  Significant: %s\n", ifelse(mt$p.value < 0.05, "YES", "NO")))

  # Effective sample size via spatial clustering
  eff_n <- NA
  ratio <- NA
  if (mt$p.value < 0.05) {
    # Cluster sites using hierarchical clustering on geographic distance
    dist_mat <- as.dist(spDists(coords, longlat = TRUE))
    hc <- hclust(dist_mat, method = "complete")
    clusters <- cutree(hc, h = threshold_km)
    eff_n <- length(unique(clusters))
    ratio <- eff_n / nrow(data)
    cat(sprintf("  Effective n: %d (%.1f%% of original %d)\n",
                eff_n, ratio * 100, nrow(data)))
    cat(sprintf("  CI widening factor: %.2fx\n", 1 / sqrt(ratio)))
  }
  cat("\n")

  return(list(
    label = label,
    n = nrow(data),
    threshold_km = threshold_km,
    moran_I = as.numeric(mt$estimate["Moran I statistic"]),
    expected_I = as.numeric(mt$estimate["Expectation"]),
    variance = as.numeric(mt$estimate["Variance"]),
    z_score = as.numeric(mt$statistic),
    p_value = mt$p.value,
    significant = mt$p.value < 0.05,
    n_no_neighbors = n_no_neighbors,
    eff_n = eff_n,
    ratio = ratio
  ))
}

# ---- 3. Run analysis ----

results <- list()

# Per-region at 50km
for (r in unique(combined$region)) {
  sub <- combined[combined$region == r, ]
  res <- run_moran(sub, paste0(r, " (50km)"), 50)
  if (!is.null(res)) results[[length(results) + 1]] <- res
}

# Combined at 50km
res <- run_moran(combined, "Combined (50km)", 50)
if (!is.null(res)) results[[length(results) + 1]] <- res

# Sensitivity: 100km threshold
cat("\n=== Sensitivity Analysis: 100km threshold ===\n\n")

for (r in unique(combined$region)) {
  sub <- combined[combined$region == r, ]
  res <- run_moran(sub, paste0(r, " (100km)"), 100)
  if (!is.null(res)) results[[length(results) + 1]] <- res
}

res <- run_moran(combined, "Combined (100km)", 100)
if (!is.null(res)) results[[length(results) + 1]] <- res

# ---- 4. Summary statistics ----

cat("\n=== Summary of signal_onset_age by region ===\n\n")
for (r in unique(combined$region)) {
  sub <- combined[combined$region == r, ]
  cat(sprintf("%s: mean=%.0f, median=%.0f, sd=%.0f, range=[%.0f, %.0f]\n",
              r, mean(sub$signal_onset_age), median(sub$signal_onset_age),
              sd(sub$signal_onset_age),
              min(sub$signal_onset_age), max(sub$signal_onset_age)))
}
all_ages <- combined$signal_onset_age
cat(sprintf("Combined: mean=%.0f, median=%.0f, sd=%.0f, range=[%.0f, %.0f]\n",
            mean(all_ages), median(all_ages), sd(all_ages),
            min(all_ages), max(all_ages)))

# ---- 5. Generate markdown report ----

md <- character()
md <- c(md, "# Phase 1: Spatial Autocorrelation of Signal Onset Ages")
md <- c(md, "")
md <- c(md, sprintf("**Date**: %s", Sys.Date()))
md <- c(md, sprintf("**Total sites analyzed**: %d (after removing NAs and coordinate duplicates)", nrow(combined)))
md <- c(md, "")

# Region counts
md <- c(md, "## Sample Sizes by Region")
md <- c(md, "")
md <- c(md, "| Region | n |")
md <- c(md, "|--------|---|")
for (r in unique(combined$region)) {
  md <- c(md, sprintf("| %s | %d |", r, sum(combined$region == r)))
}
md <- c(md, sprintf("| **Combined** | **%d** |", nrow(combined)))
md <- c(md, "")

# Moran's I table
md <- c(md, "## Moran's I Results")
md <- c(md, "")
md <- c(md, "| Analysis | n | Moran's I | Expected I | Z-score | p-value | Significant | Isolated sites |")
md <- c(md, "|----------|---|-----------|------------|---------|---------|-------------|----------------|")
for (res in results) {
  sig_str <- ifelse(res$significant, "**YES**", "no")
  md <- c(md, sprintf("| %s | %d | %.4f | %.4f | %.3f | %.4f | %s | %d |",
                       res$label, res$n, res$moran_I, res$expected_I,
                       res$z_score, res$p_value, sig_str, res$n_no_neighbors))
}
md <- c(md, "")

# Effective sample size table
sig_results <- Filter(function(x) x$significant, results)
if (length(sig_results) > 0) {
  md <- c(md, "## Effective Sample Sizes (where autocorrelation detected)")
  md <- c(md, "")
  md <- c(md, "| Analysis | Original n | Effective n | Ratio | CI widening |")
  md <- c(md, "|----------|-----------|-------------|-------|-------------|")
  for (res in sig_results) {
    ci_widen <- sprintf("%.2fx", 1 / sqrt(res$ratio))
    md <- c(md, sprintf("| %s | %d | %d | %.1f%% | %s |",
                         res$label, res$n, res$eff_n,
                         res$ratio * 100, ci_widen))
  }
  md <- c(md, "")
}

# Descriptive stats
md <- c(md, "## Descriptive Statistics: signal_onset_age (cal BP)")
md <- c(md, "")
md <- c(md, "| Region | Mean | Median | SD | Min | Max |")
md <- c(md, "|--------|------|--------|-----|-----|-----|")
for (r in unique(combined$region)) {
  sub <- combined[combined$region == r, ]
  md <- c(md, sprintf("| %s | %.0f | %.0f | %.0f | %.0f | %.0f |",
                       r, mean(sub$signal_onset_age), median(sub$signal_onset_age),
                       sd(sub$signal_onset_age),
                       min(sub$signal_onset_age), max(sub$signal_onset_age)))
}
md <- c(md, sprintf("| **Combined** | **%.0f** | **%.0f** | **%.0f** | **%.0f** | **%.0f** |",
                     mean(all_ages), median(all_ages), sd(all_ages),
                     min(all_ages), max(all_ages)))
md <- c(md, "")

# Interpretation
md <- c(md, "## Interpretation")
md <- c(md, "")

any_sig <- any(sapply(results, function(x) x$significant))
if (any_sig) {
  md <- c(md, "Significant spatial autocorrelation was detected, confirming the reviewer's concern")
  md <- c(md, "that nearby pollen sites are not fully independent. Sites within close proximity")
  md <- c(md, "share overlapping pollen source areas, respond to the same regional climate forcing,")
  md <- c(md, "and reflect activities of the same human populations.")
  md <- c(md, "")
  md <- c(md, "### Implications for confidence intervals")
  md <- c(md, "")
  md <- c(md, "The effective sample size is smaller than the nominal count. Confidence intervals")
  md <- c(md, "for regional mean onset ages should be widened by the CI widening factor shown above.")
  md <- c(md, "For example, if the CI widening factor is 1.5x, a 95% CI of +/-200 years becomes +/-300 years.")
  md <- c(md, "")
  md <- c(md, "### Recommended adjustments")
  md <- c(md, "")
  md <- c(md, "1. Report effective sample sizes alongside nominal counts")
  md <- c(md, "2. Use spatial bootstrap or block bootstrap for uncertainty estimation")
  md <- c(md, "3. Consider mixed-effects models with spatial random effects")
  md <- c(md, "4. The 50km threshold corresponds roughly to pollen source areas for arboreal taxa")
} else {
  md <- c(md, "No significant spatial autocorrelation was detected at the tested thresholds.")
  md <- c(md, "Sites can be treated as approximately independent for statistical inference.")
}
md <- c(md, "")

# Method note
md <- c(md, "## Methods")
md <- c(md, "")
md <- c(md, "- **Moran's I**: Computed using `spdep::moran.test()` with row-standardized weights")
md <- c(md, "- **Neighbor definition**: Distance-based (`dnearneigh`), great-circle distance, thresholds 50km and 100km")
md <- c(md, "- **Effective n**: Hierarchical clustering (complete linkage) on geographic distance matrix, cut at threshold distance")
md <- c(md, "- **CI widening factor**: 1 / sqrt(effective_n / original_n)")
md <- c(md, "- **Zero policy**: Sites with no neighbors within threshold are included (zero.policy=TRUE)")
md <- c(md, "")

writeLines(md, "/home/ayu/archeco/shared/phase1_spatial_autocorrelation_results.md")
cat("\n\nReport written to shared/phase1_spatial_autocorrelation_results.md\n")
