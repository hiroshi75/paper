#!/usr/bin/env Rscript
# =============================================================================
# Analysis 1: Formal Lag Estimation Between Agricultural Indicator Rise
#              and AP Decline, with Bootstrap Confidence Intervals
#
# Paper 3 (Dual Signal) — Additional Analysis
#
# For each of 3 European regions (Britain, Scandinavia, Alps):
#   1. Compute composite AP and Cerealia timeseries in 200-yr bins
#   2. Identify "indicator rise point": first bin where Cerealia > 2x
#      pre-agricultural baseline mean
#   3. Identify "AP decline point": first bin where AP < pre-agricultural
#      mean - 1 SD
#   4. Lag = AP decline point - indicator rise point
#   5. Bootstrap CI (resample sites within region, 1000 iterations)
#
# Dependencies: neotoma2, dplyr, tidyr, boot
# Output: shared/analysis/output/lag_estimation_results.csv
# =============================================================================

library(neotoma2)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BOOT_N <- 1000
BIN_WIDTH <- 200          # years
MIN_SITES_PER_BIN <- 3    # minimum contributing sites
OUTPUT_DIR <- "output"

# Region definitions: name, bounding box (W, S, E, N), agriculture start (cal BP)
regions <- list(
  britain = list(
    name = "Britain",
    bbox = c(-6, 50, 2, 57),
    ag_start = 6000   # Neolithic introduction
  ),
  scandinavia = list(
    name = "Scandinavia",
    bbox = c(5, 55, 30, 70),
    ag_start = 5500   # Funnel Beaker Culture
  ),
  alps = list(
    name = "Alps",
    bbox = c(6, 45, 16, 48),
    ag_start = 7500   # LBK
  )
)

# Neotoma ecological groups defining arboreal pollen
AP_GROUPS <- c("TRSH", "MANG", "PALM")
NAP_GROUPS <- c("UPHE", "SUCC")

# Agricultural indicator taxa (Cerealia-type and variants in Neotoma)
CEREALIA_TAXA <- c("Cerealia", "Cerealia-type", "Cerealia undiff.",
                   "Triticum", "Hordeum", "Secale", "Avena",
                   "Triticum-type", "Hordeum-type", "Secale-type")

# 200-year bin edges (0 to 10000 cal BP)
bin_edges <- seq(0, 10000, by = BIN_WIDTH)
bin_mids  <- bin_edges[-length(bin_edges)] + BIN_WIDTH / 2

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

#' Fetch pollen data from Neotoma for a given bounding box
#' Returns a list of per-site data frames with columns:
#'   site_id, age, ap_ratio, cerealia_pct
fetch_pollen_data <- function(bbox) {
  # Query Neotoma with bounding box
  # bbox format for neotoma2: c(xmin, ymin, xmax, ymax)
  sites <- get_sites(loc = bbox)

  if (length(sites) == 0) {
    warning("No sites found in bounding box: ", paste(bbox, collapse = ", "))
    return(list())
  }

  # Get datasets — filter to pollen
  datasets <- get_datasets(sites, datasettype = "pollen")

  if (length(datasets) == 0) {
    warning("No pollen datasets found")
    return(list())
  }

  # Download full sample data
  dl <- get_downloads(datasets)

  # Process each site
  site_data_list <- list()

  for (i in seq_along(dl)) {
    tryCatch({
      samples_df <- samples(dl[[i]])

      if (nrow(samples_df) == 0) next

      # Need age and ecological group / taxon info
      # Compute AP ratio and Cerealia % per sample
      sample_results <- samples_df %>%
        mutate(
          is_ap  = ecologicalgroup %in% AP_GROUPS,
          is_nap = ecologicalgroup %in% NAP_GROUPS,
          is_cerealia = variablename %in% CEREALIA_TAXA
        ) %>%
        group_by(sampleid, age) %>%
        summarise(
          ap_sum      = sum(value[is_ap], na.rm = TRUE),
          nap_sum     = sum(value[is_nap], na.rm = TRUE),
          total_sum   = ap_sum + nap_sum,
          cerealia_ct = sum(value[is_cerealia], na.rm = TRUE),
          total_pollen = sum(value, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        filter(total_sum > 10, !is.na(age)) %>%
        mutate(
          ap_ratio     = ap_sum / total_sum,
          cerealia_pct = (cerealia_ct / total_pollen) * 100
        )

      if (nrow(sample_results) < 3) next

      sid <- unique(dl[[i]]@sites@sites$siteid)[1]
      sample_results$site_id <- sid
      site_data_list[[length(site_data_list) + 1]] <- sample_results

    }, error = function(e) {
      message("Skipping site ", i, ": ", e$message)
    })
  }

  return(site_data_list)
}


#' Assign each sample to a 200-year bin based on its age (cal BP)
assign_bin <- function(age) {
  idx <- findInterval(age, bin_edges, rightmost.closed = TRUE)
  idx <- pmin(pmax(idx, 1), length(bin_mids))
  bin_mids[idx]
}


#' Build regional composites from a list of site data frames
#' Returns a data frame with columns: bin_mid, mean_ap, mean_cerealia, n_sites
build_composite <- function(site_data_list) {
  # Combine all sites
  all_samples <- bind_rows(site_data_list)

  # Assign bins
  all_samples$bin <- assign_bin(all_samples$age)

  # Average per site per bin first (site-level binning),
  # then average across sites
  site_bins <- all_samples %>%
    group_by(site_id, bin) %>%
    summarise(
      ap_ratio     = mean(ap_ratio, na.rm = TRUE),
      cerealia_pct = mean(cerealia_pct, na.rm = TRUE),
      .groups = "drop"
    )

  composite <- site_bins %>%
    group_by(bin) %>%
    summarise(
      mean_ap      = mean(ap_ratio, na.rm = TRUE),
      mean_cerealia = mean(cerealia_pct, na.rm = TRUE),
      n_sites      = n_distinct(site_id),
      .groups = "drop"
    ) %>%
    filter(n_sites >= MIN_SITES_PER_BIN) %>%
    arrange(desc(bin))   # oldest first (highest cal BP)

  return(composite)
}


#' Identify the indicator rise point:
#' First 200-yr bin where Cerealia > 2x pre-agricultural baseline mean.
#' Pre-agricultural baseline = all bins older than ag_start.
#' Returns the bin midpoint (cal BP), or NA if not found.
find_indicator_rise <- function(composite, ag_start) {
  baseline <- composite %>% filter(bin > ag_start)

  if (nrow(baseline) == 0) {
    warning("No pre-agricultural bins found for baseline")
    return(NA_real_)
  }

  baseline_mean <- mean(baseline$mean_cerealia, na.rm = TRUE)
  threshold <- 2 * baseline_mean

  # Search forward in time (decreasing cal BP) from ag_start
  post_ag <- composite %>%
    filter(bin <= ag_start) %>%
    arrange(desc(bin))   # oldest post-ag bin first

  rise_bin <- post_ag %>%
    filter(mean_cerealia > threshold) %>%
    slice(1)

  if (nrow(rise_bin) == 0) return(NA_real_)
  return(rise_bin$bin)
}


#' Identify the AP decline point:
#' First bin where AP < pre-agricultural mean - 1 SD.
#' Returns the bin midpoint (cal BP), or NA if not found.
find_ap_decline <- function(composite, ag_start) {
  baseline <- composite %>% filter(bin > ag_start)

  if (nrow(baseline) < 2) {
    warning("Insufficient pre-agricultural bins for SD calculation")
    return(NA_real_)
  }

  baseline_mean <- mean(baseline$mean_ap, na.rm = TRUE)
  baseline_sd   <- sd(baseline$mean_ap, na.rm = TRUE)
  threshold     <- baseline_mean - baseline_sd

  # Search forward in time (decreasing cal BP) from ag_start
  post_ag <- composite %>%
    filter(bin <= ag_start) %>%
    arrange(desc(bin))

  decline_bin <- post_ag %>%
    filter(mean_ap < threshold) %>%
    slice(1)

  if (nrow(decline_bin) == 0) return(NA_real_)
  return(decline_bin$bin)
}


#' Bootstrap lag estimation.
#' Resamples sites with replacement within a region, rebuilds composites,
#' re-identifies indicator rise and AP decline points, computes lag.
#' Returns a vector of lag values (length = n_boot).
bootstrap_lag <- function(site_data_list, ag_start, n_boot = BOOT_N) {
  n_sites <- length(site_data_list)
  lags <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    # Resample sites with replacement
    idx <- sample(seq_len(n_sites), size = n_sites, replace = TRUE)
    resampled <- site_data_list[idx]

    # Rebuild composite
    comp <- tryCatch(
      build_composite(resampled),
      error = function(e) NULL
    )

    if (is.null(comp) || nrow(comp) < 5) {
      lags[b] <- NA_real_
      next
    }

    rise <- find_indicator_rise(comp, ag_start)
    decline <- find_ap_decline(comp, ag_start)

    if (is.na(rise) || is.na(decline)) {
      lags[b] <- NA_real_
    } else {
      # Lag in years: rise is older (higher cal BP) than decline
      # Positive lag = AP decline comes after indicator rise
      lags[b] <- rise - decline
    }
  }

  return(lags)
}

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

results <- data.frame(
  region = character(),
  indicator_rise_calBP = numeric(),
  ap_decline_calBP = numeric(),
  lag_years = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  boot_valid_n = integer(),
  stringsAsFactors = FALSE
)

for (reg_key in names(regions)) {
  reg <- regions[[reg_key]]
  cat("\n=== Processing:", reg$name, "===\n")

  # Step 1: Fetch pollen data
  cat("Fetching pollen data from Neotoma...\n")
  site_data <- fetch_pollen_data(reg$bbox)
  cat("  Retrieved", length(site_data), "sites with usable data\n")

  if (length(site_data) < MIN_SITES_PER_BIN) {
    warning("Insufficient sites for ", reg$name, ". Skipping.")
    next
  }

  # Step 2: Build composite
  cat("Building 200-yr bin composite...\n")
  composite <- build_composite(site_data)
  cat("  Composite bins:", nrow(composite), "\n")

  # Step 3: Identify indicator rise and AP decline points
  rise_date <- find_indicator_rise(composite, reg$ag_start)
  decline_date <- find_ap_decline(composite, reg$ag_start)

  cat("  Indicator rise point:", rise_date, "cal BP\n")
  cat("  AP decline point:", decline_date, "cal BP\n")

  if (!is.na(rise_date) && !is.na(decline_date)) {
    point_lag <- rise_date - decline_date
    cat("  Point estimate lag:", point_lag, "years\n")
  } else {
    point_lag <- NA_real_
    cat("  Could not compute lag (missing rise or decline point)\n")
  }

  # Step 4: Bootstrap CI
  cat("Running bootstrap (n =", BOOT_N, ")...\n")
  boot_lags <- bootstrap_lag(site_data, reg$ag_start, BOOT_N)
  valid_lags <- boot_lags[!is.na(boot_lags)]
  cat("  Valid bootstrap iterations:", length(valid_lags), "/", BOOT_N, "\n")

  if (length(valid_lags) >= 50) {
    ci <- quantile(valid_lags, probs = c(0.025, 0.975))
    cat("  95% CI:", ci[1], "-", ci[2], "years\n")
  } else {
    ci <- c(NA_real_, NA_real_)
    cat("  Insufficient valid bootstraps for CI\n")
  }

  # Collect results
  results <- rbind(results, data.frame(
    region = reg$name,
    indicator_rise_calBP = ifelse(is.na(rise_date), NA, rise_date),
    ap_decline_calBP     = ifelse(is.na(decline_date), NA, decline_date),
    lag_years            = ifelse(is.na(point_lag), NA, point_lag),
    ci_lower             = ci[1],
    ci_upper             = ci[2],
    boot_valid_n         = length(valid_lags),
    stringsAsFactors = FALSE
  ))
}

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

cat("\n\n========================================\n")
cat("LAG ESTIMATION RESULTS\n")
cat("========================================\n\n")
print(results, row.names = FALSE)

out_file <- file.path(OUTPUT_DIR, "lag_estimation_results.csv")
write.csv(results, out_file, row.names = FALSE)
cat("\nResults saved to:", out_file, "\n")

# Also save the bootstrap distributions for diagnostic plots
cat("\nAnalysis complete.\n")
