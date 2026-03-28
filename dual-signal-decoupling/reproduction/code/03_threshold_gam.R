#!/usr/bin/env Rscript
# =============================================================================
# Analysis 3: Threshold / Nonlinearity Analysis (GAM + Segmented Regression)
#
# Paper 3 (Dual Signal) — Additional Analysis
#
# 1. Pool all 200-year bins from all 3 European regions
# 2. Fit GAM: AP_change ~ s(Cerealia_level) + s(Plantago_level) + region
#    where AP_change = AP(t) - AP(t-1)
# 3. Find Cerealia level where GAM smooth crosses zero (AP starts declining)
# 4. Fit segmented regression: breakpoint in AP_change ~ Cerealia_level
# 5. Output: threshold Cerealia % at which AP responds, with CI
#
# Dependencies: neotoma2, dplyr, tidyr, mgcv, segmented
# Output: shared/analysis/output/threshold_gam_results.csv
#         shared/analysis/output/threshold_gam_summary.txt
# =============================================================================

library(neotoma2)
library(dplyr)
library(tidyr)
library(mgcv)

# segmented is optional — install if needed
if (!requireNamespace("segmented", quietly = TRUE)) {
  message("Package 'segmented' not installed. Segmented regression will be skipped.")
  message("Install with: install.packages('segmented')")
  HAS_SEGMENTED <- FALSE
} else {
  library(segmented)
  HAS_SEGMENTED <- TRUE
}

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BIN_WIDTH <- 200
MIN_SITES_PER_BIN <- 3
OUTPUT_DIR <- "output"

regions <- list(
  britain = list(
    name = "Britain",
    bbox = c(-6, 50, 2, 57),
    ag_start = 6000
  ),
  scandinavia = list(
    name = "Scandinavia",
    bbox = c(5, 55, 30, 70),
    ag_start = 5500
  ),
  alps = list(
    name = "Alps",
    bbox = c(6, 45, 16, 48),
    ag_start = 7500
  )
)

AP_GROUPS <- c("TRSH", "MANG", "PALM")
NAP_GROUPS <- c("UPHE", "SUCC")

CEREALIA_TAXA <- c("Cerealia", "Cerealia-type", "Cerealia undiff.",
                   "Triticum", "Hordeum", "Secale", "Avena",
                   "Triticum-type", "Hordeum-type", "Secale-type")

PLANTAGO_TAXA <- c("Plantago lanceolata", "Plantago lanceolata-type",
                   "Plantago")

bin_edges <- seq(0, 10000, by = BIN_WIDTH)
bin_mids  <- bin_edges[-length(bin_edges)] + BIN_WIDTH / 2

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

assign_bin <- function(age) {
  idx <- findInterval(age, bin_edges, rightmost.closed = TRUE)
  idx <- pmin(pmax(idx, 1), length(bin_mids))
  bin_mids[idx]
}


#' Fetch pollen data from Neotoma for a region.
#' Returns site-level binned data with AP, Cerealia, and Plantago.
fetch_and_bin_data <- function(bbox) {
  sites <- get_sites(loc = bbox)
  if (length(sites) == 0) return(data.frame())

  datasets <- get_datasets(sites, datasettype = "pollen")
  if (length(datasets) == 0) return(data.frame())

  dl <- get_downloads(datasets)

  all_results <- list()

  for (i in seq_along(dl)) {
    tryCatch({
      samples_df <- samples(dl[[i]])
      if (nrow(samples_df) == 0) next

      sample_results <- samples_df %>%
        mutate(
          is_ap       = ecologicalgroup %in% AP_GROUPS,
          is_nap      = ecologicalgroup %in% NAP_GROUPS,
          is_cerealia = variablename %in% CEREALIA_TAXA,
          is_plantago = variablename %in% PLANTAGO_TAXA
        ) %>%
        group_by(sampleid, age) %>%
        summarise(
          ap_sum       = sum(value[is_ap], na.rm = TRUE),
          nap_sum      = sum(value[is_nap], na.rm = TRUE),
          total_sum    = ap_sum + nap_sum,
          cerealia_ct  = sum(value[is_cerealia], na.rm = TRUE),
          plantago_ct  = sum(value[is_plantago], na.rm = TRUE),
          total_pollen = sum(value, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        filter(total_sum > 10, !is.na(age)) %>%
        mutate(
          ap_ratio      = ap_sum / total_sum,
          cerealia_pct  = (cerealia_ct / total_pollen) * 100,
          plantago_pct  = (plantago_ct / total_pollen) * 100,
          bin           = assign_bin(age)
        )

      if (nrow(sample_results) < 3) next

      sid <- unique(dl[[i]]@sites@sites$siteid)[1]
      sample_results$site_id <- sid
      all_results[[length(all_results) + 1]] <- sample_results

    }, error = function(e) {
      message("Skipping site ", i, ": ", e$message)
    })
  }

  if (length(all_results) == 0) return(data.frame())
  bind_rows(all_results)
}


#' Build regional composite with AP, Cerealia, and Plantago
build_composite_full <- function(region_df) {
  # Site-level binning first
  site_bins <- region_df %>%
    group_by(site_id, bin) %>%
    summarise(
      ap_ratio     = mean(ap_ratio, na.rm = TRUE),
      cerealia_pct = mean(cerealia_pct, na.rm = TRUE),
      plantago_pct = mean(plantago_pct, na.rm = TRUE),
      .groups = "drop"
    )

  # Regional composite
  composite <- site_bins %>%
    group_by(bin) %>%
    summarise(
      mean_ap       = mean(ap_ratio, na.rm = TRUE),
      mean_cerealia = mean(cerealia_pct, na.rm = TRUE),
      mean_plantago = mean(plantago_pct, na.rm = TRUE),
      n_sites       = n_distinct(site_id),
      .groups = "drop"
    ) %>%
    filter(n_sites >= MIN_SITES_PER_BIN) %>%
    arrange(desc(bin))   # oldest first

  return(composite)
}

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Step 1: Fetch and build composites for all regions
all_composites <- list()

for (reg_key in names(regions)) {
  reg <- regions[[reg_key]]
  cat("\n=== Fetching:", reg$name, "===\n")

  region_df <- fetch_and_bin_data(reg$bbox)
  if (nrow(region_df) == 0) {
    cat("  No data. Skipping.\n")
    next
  }

  composite <- build_composite_full(region_df)
  composite$region <- reg$name
  all_composites[[reg_key]] <- composite
  cat("  Composite bins:", nrow(composite), "\n")
}

# Step 2: Compute AP_change = AP(t) - AP(t-1) for each region
# Bins are ordered oldest-first (descending cal BP).
# "Forward in time" = moving down the rows.
# AP_change(t) = AP(t) - AP(t-1), where t is the younger bin.

pooled_data <- data.frame()

for (reg_key in names(all_composites)) {
  comp <- all_composites[[reg_key]]

  if (nrow(comp) < 2) next

  # comp is ordered oldest-first; compute change going forward in time
  comp <- comp %>%
    arrange(desc(bin)) %>%   # ensure oldest first
    mutate(
      ap_change = mean_ap - lag(mean_ap),   # positive = AP increase
      cerealia_level = mean_cerealia,
      plantago_level = mean_plantago
    ) %>%
    filter(!is.na(ap_change))   # drop first row (no lag)

  pooled_data <- bind_rows(pooled_data, comp)
}

cat("\n=== Pooled dataset ===\n")
cat("Total bins with AP_change:", nrow(pooled_data), "\n")
cat("Regions:", paste(unique(pooled_data$region), collapse = ", "), "\n")

# Save pooled data
pooled_file <- file.path(OUTPUT_DIR, "threshold_pooled_data.csv")
write.csv(pooled_data, pooled_file, row.names = FALSE)
cat("Pooled data saved to:", pooled_file, "\n")

# ---------------------------------------------------------------------------
# Step 3: Fit GAM
# AP_change ~ s(cerealia_level) + s(plantago_level) + region
# ---------------------------------------------------------------------------

cat("\n=== Fitting GAM ===\n")

# Ensure region is a factor
pooled_data$region <- as.factor(pooled_data$region)

# Check for sufficient unique values of predictors
n_unique_cer <- length(unique(pooled_data$cerealia_level))
n_unique_pla <- length(unique(pooled_data$plantago_level))
cat("Unique Cerealia levels:", n_unique_cer, "\n")
cat("Unique Plantago levels:", n_unique_pla, "\n")

# Set appropriate basis dimension (k) — cannot exceed unique values
k_cer <- min(10, n_unique_cer - 1)
k_pla <- min(10, n_unique_pla - 1)

# Fit GAM with adaptive smoothing
# If too few unique values, fall back to simpler model
if (k_cer >= 4 && k_pla >= 4) {
  gam_fit <- gam(
    ap_change ~ s(cerealia_level, k = k_cer, bs = "tp") +
                s(plantago_level, k = k_pla, bs = "tp") +
                region,
    data = pooled_data,
    method = "REML"
  )
  cat("\nFull GAM fitted successfully.\n")
} else if (k_cer >= 4) {
  # Plantago has too few unique values — use linear
  gam_fit <- gam(
    ap_change ~ s(cerealia_level, k = k_cer, bs = "tp") +
                plantago_level +
                region,
    data = pooled_data,
    method = "REML"
  )
  cat("\nGAM fitted with linear Plantago term (insufficient unique values for smooth).\n")
} else {
  # Both have too few — use linear terms
  gam_fit <- gam(
    ap_change ~ cerealia_level + plantago_level + region,
    data = pooled_data,
    method = "REML"
  )
  cat("\nLinear model fitted (insufficient unique values for smooth terms).\n")
}

cat("\nGAM Summary:\n")
print(summary(gam_fit))

# ---------------------------------------------------------------------------
# Step 4: Find the Cerealia threshold where the smooth term crosses zero
# ---------------------------------------------------------------------------

cat("\n=== Finding Cerealia threshold (GAM smooth zero-crossing) ===\n")

# Create prediction grid: vary Cerealia, hold Plantago and region at median/mode
cer_range <- seq(min(pooled_data$cerealia_level),
                 max(pooled_data$cerealia_level),
                 length.out = 500)

# Use median Plantago and most common region for prediction
median_plantago <- median(pooled_data$plantago_level)
mode_region <- names(sort(table(pooled_data$region), decreasing = TRUE))[1]

pred_df <- data.frame(
  cerealia_level = cer_range,
  plantago_level = median_plantago,
  region = factor(mode_region, levels = levels(pooled_data$region))
)

# Predict with SE for CI
pred <- predict(gam_fit, newdata = pred_df, se.fit = TRUE)
pred_df$fit <- pred$fit
pred_df$se  <- pred$se.fit
pred_df$lower <- pred$fit - 1.96 * pred$se.fit
pred_df$upper <- pred$fit + 1.96 * pred$se.fit

# Find zero-crossing of the fitted curve
# (transition from positive to negative AP_change)
zero_crossings <- which(diff(sign(pred_df$fit)) != 0)

gam_threshold <- NA_real_
gam_threshold_ci_lower <- NA_real_
gam_threshold_ci_upper <- NA_real_

if (length(zero_crossings) > 0) {
  # Take the first crossing where fit goes from positive to negative
  for (zc in zero_crossings) {
    if (pred_df$fit[zc] > 0 && pred_df$fit[zc + 1] <= 0) {
      # Linear interpolation for precise crossing
      x1 <- pred_df$cerealia_level[zc]
      x2 <- pred_df$cerealia_level[zc + 1]
      y1 <- pred_df$fit[zc]
      y2 <- pred_df$fit[zc + 1]
      gam_threshold <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)

      # Approximate CI from upper/lower bounds
      upper_cross <- which(diff(sign(pred_df$upper)) != 0)
      lower_cross <- which(diff(sign(pred_df$lower)) != 0)

      # Find relevant crossings near the main threshold
      if (length(upper_cross) > 0) {
        uc <- upper_cross[which.min(abs(pred_df$cerealia_level[upper_cross] -
                                        gam_threshold))]
        gam_threshold_ci_upper <- pred_df$cerealia_level[uc]
      }
      if (length(lower_cross) > 0) {
        lc <- lower_cross[which.min(abs(pred_df$cerealia_level[lower_cross] -
                                        gam_threshold))]
        gam_threshold_ci_lower <- pred_df$cerealia_level[lc]
      }

      break
    }
  }
}

if (!is.na(gam_threshold)) {
  cat("GAM Cerealia threshold (zero-crossing):", round(gam_threshold, 4), "%\n")
  cat("  Approximate 95% CI: [",
      round(min(gam_threshold_ci_lower, gam_threshold_ci_upper, na.rm = TRUE), 4),
      ",",
      round(max(gam_threshold_ci_lower, gam_threshold_ci_upper, na.rm = TRUE), 4),
      "] %\n")
} else {
  cat("No zero-crossing found in GAM smooth for Cerealia.\n")
  cat("This may indicate a monotonic relationship or insufficient data range.\n")
}

# Save prediction curve
pred_file <- file.path(OUTPUT_DIR, "threshold_gam_predictions.csv")
write.csv(pred_df, pred_file, row.names = FALSE)

# ---------------------------------------------------------------------------
# Step 5: Segmented regression
# ---------------------------------------------------------------------------

cat("\n=== Segmented Regression ===\n")

seg_threshold <- NA_real_
seg_ci_lower  <- NA_real_
seg_ci_upper  <- NA_real_

if (HAS_SEGMENTED) {
  tryCatch({
    # Fit base linear model
    lm_fit <- lm(ap_change ~ cerealia_level + region, data = pooled_data)

    # Fit segmented model with one breakpoint
    # Use median as initial guess for breakpoint
    seg_fit <- segmented(lm_fit,
                         seg.Z = ~ cerealia_level,
                         psi = median(pooled_data$cerealia_level))

    cat("\nSegmented regression summary:\n")
    print(summary(seg_fit))

    # Extract breakpoint and CI
    bp <- seg_fit$psi
    seg_threshold <- bp[1, "Est."]
    seg_ci_lower  <- confint(seg_fit)$cerealia_level[1, 1]
    seg_ci_upper  <- confint(seg_fit)$cerealia_level[1, 2]

    cat("\nSegmented regression breakpoint:", round(seg_threshold, 4), "%\n")
    cat("  95% CI: [", round(seg_ci_lower, 4), ",",
        round(seg_ci_upper, 4), "] %\n")

  }, error = function(e) {
    cat("Segmented regression failed:", e$message, "\n")
    cat("This may occur if the relationship is not piecewise linear.\n")
  })
} else {
  cat("Package 'segmented' not available. Skipping segmented regression.\n")
}

# ---------------------------------------------------------------------------
# Output summary
# ---------------------------------------------------------------------------

cat("\n\n========================================\n")
cat("THRESHOLD ANALYSIS RESULTS\n")
cat("========================================\n\n")

results <- data.frame(
  method = c("GAM smooth zero-crossing", "Segmented regression"),
  threshold_cerealia_pct = c(
    ifelse(is.na(gam_threshold), NA, round(gam_threshold, 4)),
    ifelse(is.na(seg_threshold), NA, round(seg_threshold, 4))
  ),
  ci_lower = c(
    ifelse(is.na(gam_threshold_ci_lower), NA,
           round(min(gam_threshold_ci_lower, gam_threshold_ci_upper,
                     na.rm = TRUE), 4)),
    ifelse(is.na(seg_ci_lower), NA, round(seg_ci_lower, 4))
  ),
  ci_upper = c(
    ifelse(is.na(gam_threshold_ci_upper), NA,
           round(max(gam_threshold_ci_lower, gam_threshold_ci_upper,
                     na.rm = TRUE), 4)),
    ifelse(is.na(seg_ci_upper), NA, round(seg_ci_upper, 4))
  ),
  stringsAsFactors = FALSE
)

print(results, row.names = FALSE)

# Save results
results_file <- file.path(OUTPUT_DIR, "threshold_gam_results.csv")
write.csv(results, results_file, row.names = FALSE)
cat("\nResults saved to:", results_file, "\n")

# Save text summary
summary_file <- file.path(OUTPUT_DIR, "threshold_gam_summary.txt")
sink(summary_file)
cat("THRESHOLD / NONLINEARITY ANALYSIS SUMMARY\n")
cat("==========================================\n\n")
cat("Data: Pooled 200-yr bins from Britain, Scandinavia, Alps\n")
cat("Total observations:", nrow(pooled_data), "\n\n")

cat("--- GAM ---\n")
cat("Model: AP_change ~ s(Cerealia) + s(Plantago) + region\n")
print(summary(gam_fit))
cat("\nCerealia threshold (zero-crossing):",
    ifelse(is.na(gam_threshold), "Not found",
           paste0(round(gam_threshold, 4), "%")), "\n")

if (HAS_SEGMENTED && !is.na(seg_threshold)) {
  cat("\n--- Segmented Regression ---\n")
  cat("Breakpoint:", round(seg_threshold, 4), "%\n")
  cat("95% CI: [", round(seg_ci_lower, 4), ",",
      round(seg_ci_upper, 4), "]\n")
}
sink()
cat("Summary saved to:", summary_file, "\n")

cat("\nAnalysis complete.\n")
