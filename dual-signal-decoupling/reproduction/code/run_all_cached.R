#!/usr/bin/env Rscript
# =============================================================================
# Paper 3 (Dual Signal) — All Three Analyses Using Cached Data
#
# Reads from cached JSON files instead of querying Neotoma API.
# Cached data: shared/neotoma_{region}_agri_indicators.json
#
# Analysis 1: Lag estimation between Cerealia rise and AP decline
# Analysis 2: Site-level replication of dual signal pattern
# Analysis 3: GAM threshold analysis (AP change ~ Cerealia level)
# =============================================================================

library(dplyr)
library(tidyr)
library(jsonlite)
library(mgcv)

# segmented is optional
if (!requireNamespace("segmented", quietly = TRUE)) {
  message("Package 'segmented' not installed. Segmented regression will be skipped.")
  HAS_SEGMENTED <- FALSE
} else {
  library(segmented)
  HAS_SEGMENTED <- TRUE
}

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE_DIR <- "/home/ayu/archeco/shared"
OUTPUT_DIR <- "/home/ayu/archeco/shared/analysis/output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

BOOT_N <- 1000
BIN_WIDTH <- 200
MIN_SITES_PER_BIN <- 3
CEREALIA_THRESHOLD <- 0.01  # % threshold for "Cerealia present"

regions <- list(
  britain = list(
    name = "Britain",
    file = file.path(BASE_DIR, "neotoma_britain_agri_indicators.json"),
    ag_start = 6000
  ),
  scandinavia = list(
    name = "Scandinavia",
    file = file.path(BASE_DIR, "neotoma_scandinavia_agri_indicators.json"),
    ag_start = 5500
  ),
  alps = list(
    name = "Alps",
    file = file.path(BASE_DIR, "neotoma_alps_agri_indicators.json"),
    ag_start = 7500
  )
)

bin_edges <- seq(0, 10000, by = BIN_WIDTH)
bin_mids  <- bin_edges[-length(bin_edges)] + BIN_WIDTH / 2

# ---------------------------------------------------------------------------
# Data loading from cached JSON
# ---------------------------------------------------------------------------

#' Load cached JSON data and return a list of per-site data frames
#' Each data frame has: site_id, site_name, age, ap_ratio, cerealia_pct,
#'                       pastoral_pct, total
load_cached_data <- function(json_file) {
  cat("  Loading:", json_file, "\n")
  raw <- fromJSON(json_file, simplifyVector = FALSE)

  site_data_list <- list()

  for (site_entry in raw) {
    samples <- site_entry$samples
    if (length(samples) < 3) next

    df <- bind_rows(lapply(samples, function(s) {
      data.frame(
        age          = as.numeric(s$age),
        ap_ratio     = as.numeric(s$ap_ratio),
        cerealia_pct = as.numeric(s$cereal_pct),
        pastoral_pct = as.numeric(s$pastoral_pct),
        disturb_pct  = as.numeric(s$disturb_pct),
        total        = as.numeric(s$total),
        stringsAsFactors = FALSE
      )
    }))

    df <- df %>% filter(!is.na(age), total > 10)

    if (nrow(df) < 3) next

    df$site_id   <- site_entry$datasetid
    df$site_name <- site_entry$sitename

    site_data_list[[length(site_data_list) + 1]] <- df
  }

  return(site_data_list)
}

# ---------------------------------------------------------------------------
# Helper functions (shared across analyses)
# ---------------------------------------------------------------------------

assign_bin <- function(age) {
  idx <- findInterval(age, bin_edges, rightmost.closed = TRUE)
  idx <- pmin(pmax(idx, 1), length(bin_mids))
  bin_mids[idx]
}

build_composite <- function(site_data_list) {
  all_samples <- bind_rows(site_data_list)
  all_samples$bin <- assign_bin(all_samples$age)

  site_bins <- all_samples %>%
    group_by(site_id, bin) %>%
    summarise(
      ap_ratio     = mean(ap_ratio, na.rm = TRUE),
      cerealia_pct = mean(cerealia_pct, na.rm = TRUE),
      pastoral_pct = mean(pastoral_pct, na.rm = TRUE),
      .groups = "drop"
    )

  composite <- site_bins %>%
    group_by(bin) %>%
    summarise(
      mean_ap       = mean(ap_ratio, na.rm = TRUE),
      mean_cerealia = mean(cerealia_pct, na.rm = TRUE),
      mean_pastoral = mean(pastoral_pct, na.rm = TRUE),
      n_sites       = n_distinct(site_id),
      .groups = "drop"
    ) %>%
    filter(n_sites >= MIN_SITES_PER_BIN) %>%
    arrange(desc(bin))

  return(composite)
}

# ============================================================================
# ANALYSIS 1: LAG ESTIMATION
# ============================================================================

cat("\n################################################################\n")
cat("# ANALYSIS 1: LAG ESTIMATION\n")
cat("################################################################\n")

find_indicator_rise <- function(composite, ag_start) {
  baseline <- composite %>% filter(bin > ag_start)
  if (nrow(baseline) == 0) return(NA_real_)

  baseline_mean <- mean(baseline$mean_cerealia, na.rm = TRUE)
  threshold <- 2 * baseline_mean

  # If baseline is essentially zero, use absolute threshold
  if (threshold < 0.01) threshold <- 0.01

  post_ag <- composite %>%
    filter(bin <= ag_start) %>%
    arrange(desc(bin))

  rise_bin <- post_ag %>%
    filter(mean_cerealia > threshold) %>%
    slice(1)

  if (nrow(rise_bin) == 0) return(NA_real_)
  return(rise_bin$bin)
}

find_ap_decline <- function(composite, ag_start) {
  baseline <- composite %>% filter(bin > ag_start)
  if (nrow(baseline) < 2) return(NA_real_)

  baseline_mean <- mean(baseline$mean_ap, na.rm = TRUE)
  baseline_sd   <- sd(baseline$mean_ap, na.rm = TRUE)
  threshold     <- baseline_mean - baseline_sd

  post_ag <- composite %>%
    filter(bin <= ag_start) %>%
    arrange(desc(bin))

  decline_bin <- post_ag %>%
    filter(mean_ap < threshold) %>%
    slice(1)

  if (nrow(decline_bin) == 0) return(NA_real_)
  return(decline_bin$bin)
}

bootstrap_lag <- function(site_data_list, ag_start, n_boot = BOOT_N) {
  n_sites <- length(site_data_list)
  lags <- numeric(n_boot)

  set.seed(42)  # reproducibility

  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n_sites), size = n_sites, replace = TRUE)
    resampled <- site_data_list[idx]

    comp <- tryCatch(build_composite(resampled), error = function(e) NULL)

    if (is.null(comp) || nrow(comp) < 5) {
      lags[b] <- NA_real_
      next
    }

    rise <- find_indicator_rise(comp, ag_start)
    decline <- find_ap_decline(comp, ag_start)

    if (is.na(rise) || is.na(decline)) {
      lags[b] <- NA_real_
    } else {
      lags[b] <- rise - decline
    }
  }

  return(lags)
}

# --- Run Analysis 1 ---

lag_results <- data.frame(
  region = character(),
  indicator_rise_calBP = numeric(),
  ap_decline_calBP = numeric(),
  lag_years = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  boot_valid_n = integer(),
  n_sites = integer(),
  stringsAsFactors = FALSE
)

all_site_data <- list()  # cache for reuse in analyses 2 and 3

for (reg_key in names(regions)) {
  reg <- regions[[reg_key]]
  cat("\n=== Processing:", reg$name, "===\n")

  site_data <- load_cached_data(reg$file)
  cat("  Loaded", length(site_data), "sites\n")
  all_site_data[[reg_key]] <- site_data

  if (length(site_data) < MIN_SITES_PER_BIN) {
    warning("Insufficient sites for ", reg$name)
    next
  }

  composite <- build_composite(site_data)
  cat("  Composite bins:", nrow(composite), "\n")

  rise_date <- find_indicator_rise(composite, reg$ag_start)
  decline_date <- find_ap_decline(composite, reg$ag_start)

  cat("  Indicator rise point:", rise_date, "cal BP\n")
  cat("  AP decline point:", decline_date, "cal BP\n")

  if (!is.na(rise_date) && !is.na(decline_date)) {
    point_lag <- rise_date - decline_date
    cat("  Point estimate lag:", point_lag, "years\n")
  } else {
    point_lag <- NA_real_
    cat("  Could not compute lag\n")
  }

  cat("  Running bootstrap (n =", BOOT_N, ")...\n")
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

  lag_results <- rbind(lag_results, data.frame(
    region               = reg$name,
    indicator_rise_calBP = ifelse(is.na(rise_date), NA, rise_date),
    ap_decline_calBP     = ifelse(is.na(decline_date), NA, decline_date),
    lag_years            = ifelse(is.na(point_lag), NA, point_lag),
    ci_lower             = ci[1],
    ci_upper             = ci[2],
    boot_valid_n         = length(valid_lags),
    n_sites              = length(site_data),
    stringsAsFactors = FALSE
  ))
}

cat("\n\n========================================\n")
cat("LAG ESTIMATION RESULTS\n")
cat("========================================\n\n")
print(lag_results, row.names = FALSE)

lag_file <- file.path(OUTPUT_DIR, "lag_estimation_results.csv")
write.csv(lag_results, lag_file, row.names = FALSE)
cat("\nSaved to:", lag_file, "\n")


# ============================================================================
# ANALYSIS 2: SITE-LEVEL REPLICATION
# ============================================================================

cat("\n\n################################################################\n")
cat("# ANALYSIS 2: SITE-LEVEL REPLICATION\n")
cat("################################################################\n")

classify_site <- function(site_df) {
  site_df$bin <- assign_bin(site_df$age)

  site_bins <- site_df %>%
    group_by(bin) %>%
    summarise(
      ap_ratio     = mean(ap_ratio, na.rm = TRUE),
      cerealia_pct = mean(cerealia_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(bin))

  if (nrow(site_bins) < 4) {
    return(list(
      classification = "no clear pattern",
      cerealia_onset = NA_real_,
      ap_decline     = NA_real_,
      reason         = "insufficient bins"
    ))
  }

  # Find Cerealia onset (oldest bin going forward in time)
  cerealia_onset <- NA_real_
  for (j in seq_len(nrow(site_bins))) {
    if (site_bins$cerealia_pct[j] > CEREALIA_THRESHOLD) {
      cerealia_onset <- site_bins$bin[j]
      break
    }
  }

  # Find AP decline onset (2 consecutive decreases)
  ap_decline <- NA_real_
  for (j in 1:(nrow(site_bins) - 2)) {
    ap1 <- site_bins$ap_ratio[j]
    ap2 <- site_bins$ap_ratio[j + 1]
    ap3 <- site_bins$ap_ratio[j + 2]
    if (ap2 < ap1 && ap3 < ap2) {
      ap_decline <- site_bins$bin[j + 1]
      break
    }
  }

  # Classify
  if (is.na(cerealia_onset) && is.na(ap_decline)) {
    return(list(classification = "no clear pattern",
                cerealia_onset = NA_real_, ap_decline = NA_real_,
                reason = "neither event detected"))
  } else if (is.na(cerealia_onset)) {
    return(list(classification = "AP decline first",
                cerealia_onset = NA_real_, ap_decline = ap_decline,
                reason = "AP decline without Cerealia"))
  } else if (is.na(ap_decline)) {
    return(list(classification = "indicator first",
                cerealia_onset = cerealia_onset, ap_decline = NA_real_,
                reason = "Cerealia without AP decline"))
  } else if (cerealia_onset == ap_decline) {
    return(list(classification = "simultaneous",
                cerealia_onset = cerealia_onset, ap_decline = ap_decline,
                reason = "same bin"))
  } else if (cerealia_onset > ap_decline) {
    return(list(classification = "indicator first",
                cerealia_onset = cerealia_onset, ap_decline = ap_decline,
                reason = paste0("Cerealia at ", cerealia_onset, " BP, AP decline at ",
                               ap_decline, " BP")))
  } else {
    return(list(classification = "AP decline first",
                cerealia_onset = cerealia_onset, ap_decline = ap_decline,
                reason = paste0("AP decline at ", ap_decline, " BP, Cerealia at ",
                               cerealia_onset, " BP")))
  }
}

# --- Run Analysis 2 ---

site_results_all <- data.frame()
summary_all <- data.frame()

for (reg_key in names(regions)) {
  reg <- regions[[reg_key]]
  cat("\n=== Processing:", reg$name, "===\n")

  site_data <- all_site_data[[reg_key]]
  cat("  Sites:", length(site_data), "\n")

  if (length(site_data) == 0) next

  site_classifications <- list()

  for (i in seq_along(site_data)) {
    sdf <- site_data[[i]]
    sid <- unique(sdf$site_id)[1]
    sname <- unique(sdf$site_name)[1]
    result <- classify_site(sdf)

    site_classifications[[i]] <- data.frame(
      region         = reg$name,
      site_id        = sid,
      site_name      = sname,
      n_bins         = length(unique(assign_bin(sdf$age))),
      classification = result$classification,
      cerealia_onset = ifelse(is.na(result$cerealia_onset), NA, result$cerealia_onset),
      ap_decline_bin = ifelse(is.na(result$ap_decline), NA, result$ap_decline),
      reason         = result$reason,
      stringsAsFactors = FALSE
    )
  }

  site_results <- bind_rows(site_classifications)
  site_results_all <- bind_rows(site_results_all, site_results)

  n_total <- nrow(site_results)
  counts <- table(site_results$classification)

  get_count <- function(name) {
    if (name %in% names(counts)) as.integer(counts[name]) else 0L
  }

  summary_row <- data.frame(
    region             = reg$name,
    total_sites        = n_total,
    n_indicator_first  = get_count("indicator first"),
    n_AP_decline_first = get_count("AP decline first"),
    n_simultaneous     = get_count("simultaneous"),
    n_no_clear_pattern = get_count("no clear pattern"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      pct_indicator_first  = round(n_indicator_first / total_sites * 100, 1),
      pct_AP_decline_first = round(n_AP_decline_first / total_sites * 100, 1),
      pct_simultaneous     = round(n_simultaneous / total_sites * 100, 1),
      pct_no_clear_pattern = round(n_no_clear_pattern / total_sites * 100, 1)
    )

  summary_all <- bind_rows(summary_all, summary_row)

  cat("\n  Classification summary:\n")
  cat("    Indicator first: ", summary_row$n_indicator_first,
      " (", summary_row$pct_indicator_first, "%)\n", sep = "")
  cat("    AP decline first:", summary_row$n_AP_decline_first,
      " (", summary_row$pct_AP_decline_first, "%)\n", sep = "")
  cat("    Simultaneous:    ", summary_row$n_simultaneous,
      " (", summary_row$pct_simultaneous, "%)\n", sep = "")
  cat("    No clear pattern:", summary_row$n_no_clear_pattern,
      " (", summary_row$pct_no_clear_pattern, "%)\n", sep = "")
}

cat("\n\n========================================\n")
cat("SITE-LEVEL REPLICATION SUMMARY\n")
cat("========================================\n\n")
print(summary_all, row.names = FALSE)

site_file <- file.path(OUTPUT_DIR, "site_replication_results.csv")
write.csv(site_results_all, site_file, row.names = FALSE)
cat("\nSite-level results saved to:", site_file, "\n")

summary_file <- file.path(OUTPUT_DIR, "site_replication_summary.csv")
write.csv(summary_all, summary_file, row.names = FALSE)
cat("Summary saved to:", summary_file, "\n")


# ============================================================================
# ANALYSIS 3: THRESHOLD GAM
# ============================================================================

cat("\n\n################################################################\n")
cat("# ANALYSIS 3: THRESHOLD / NONLINEARITY (GAM)\n")
cat("################################################################\n")

# Build composites with pastoral indicator (proxy for Plantago)
all_composites <- list()

for (reg_key in names(regions)) {
  reg <- regions[[reg_key]]
  site_data <- all_site_data[[reg_key]]

  if (length(site_data) == 0) next

  composite <- build_composite(site_data)
  composite$region <- reg$name
  all_composites[[reg_key]] <- composite
  cat("  ", reg$name, "- composite bins:", nrow(composite), "\n")
}

# Compute AP_change for each region
pooled_data <- data.frame()

for (reg_key in names(all_composites)) {
  comp <- all_composites[[reg_key]]
  if (nrow(comp) < 2) next

  comp <- comp %>%
    arrange(desc(bin)) %>%
    mutate(
      ap_change      = mean_ap - lag(mean_ap),
      cerealia_level = mean_cerealia,
      plantago_level = mean_pastoral   # using pastoral_pct as Plantago proxy
    ) %>%
    filter(!is.na(ap_change))

  pooled_data <- bind_rows(pooled_data, comp)
}

cat("\nPooled dataset: ", nrow(pooled_data), " bins with AP_change\n")
cat("Regions:", paste(unique(pooled_data$region), collapse = ", "), "\n")

pooled_file <- file.path(OUTPUT_DIR, "threshold_pooled_data.csv")
write.csv(pooled_data, pooled_file, row.names = FALSE)

# --- Fit GAM ---

cat("\n=== Fitting GAM ===\n")

pooled_data$region <- as.factor(pooled_data$region)

n_unique_cer <- length(unique(pooled_data$cerealia_level))
n_unique_pla <- length(unique(pooled_data$plantago_level))
cat("Unique Cerealia levels:", n_unique_cer, "\n")
cat("Unique Plantago (pastoral) levels:", n_unique_pla, "\n")

k_cer <- min(10, n_unique_cer - 1)
k_pla <- min(10, n_unique_pla - 1)

if (k_cer >= 4 && k_pla >= 4) {
  gam_fit <- gam(
    ap_change ~ s(cerealia_level, k = k_cer, bs = "tp") +
                s(plantago_level, k = k_pla, bs = "tp") +
                region,
    data = pooled_data,
    method = "REML"
  )
  cat("Full GAM fitted.\n")
} else if (k_cer >= 4) {
  gam_fit <- gam(
    ap_change ~ s(cerealia_level, k = k_cer, bs = "tp") +
                plantago_level + region,
    data = pooled_data,
    method = "REML"
  )
  cat("GAM with linear Plantago term.\n")
} else {
  gam_fit <- gam(
    ap_change ~ cerealia_level + plantago_level + region,
    data = pooled_data,
    method = "REML"
  )
  cat("Linear model (insufficient unique values for smooths).\n")
}

cat("\nGAM Summary:\n")
print(summary(gam_fit))

# --- Find threshold ---

cat("\n=== Finding Cerealia threshold ===\n")

cer_range <- seq(min(pooled_data$cerealia_level),
                 max(pooled_data$cerealia_level),
                 length.out = 500)

median_plantago <- median(pooled_data$plantago_level)
mode_region <- names(sort(table(pooled_data$region), decreasing = TRUE))[1]

pred_df <- data.frame(
  cerealia_level = cer_range,
  plantago_level = median_plantago,
  region = factor(mode_region, levels = levels(pooled_data$region))
)

pred <- predict(gam_fit, newdata = pred_df, se.fit = TRUE)
pred_df$fit   <- pred$fit
pred_df$se    <- pred$se.fit
pred_df$lower <- pred$fit - 1.96 * pred$se.fit
pred_df$upper <- pred$fit + 1.96 * pred$se.fit

# Find zero-crossing
zero_crossings <- which(diff(sign(pred_df$fit)) != 0)

gam_threshold <- NA_real_
gam_threshold_ci_lower <- NA_real_
gam_threshold_ci_upper <- NA_real_

if (length(zero_crossings) > 0) {
  for (zc in zero_crossings) {
    if (pred_df$fit[zc] > 0 && pred_df$fit[zc + 1] <= 0) {
      x1 <- pred_df$cerealia_level[zc]
      x2 <- pred_df$cerealia_level[zc + 1]
      y1 <- pred_df$fit[zc]
      y2 <- pred_df$fit[zc + 1]
      gam_threshold <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)

      upper_cross <- which(diff(sign(pred_df$upper)) != 0)
      lower_cross <- which(diff(sign(pred_df$lower)) != 0)

      if (length(upper_cross) > 0) {
        uc <- upper_cross[which.min(abs(pred_df$cerealia_level[upper_cross] - gam_threshold))]
        gam_threshold_ci_upper <- pred_df$cerealia_level[uc]
      }
      if (length(lower_cross) > 0) {
        lc <- lower_cross[which.min(abs(pred_df$cerealia_level[lower_cross] - gam_threshold))]
        gam_threshold_ci_lower <- pred_df$cerealia_level[lc]
      }
      break
    }
  }
}

if (!is.na(gam_threshold)) {
  cat("GAM Cerealia threshold:", round(gam_threshold, 4), "%\n")
  cat("  Approx 95% CI: [",
      round(min(gam_threshold_ci_lower, gam_threshold_ci_upper, na.rm = TRUE), 4), ",",
      round(max(gam_threshold_ci_lower, gam_threshold_ci_upper, na.rm = TRUE), 4), "] %\n")
} else {
  cat("No zero-crossing found in GAM smooth.\n")
  cat("Checking if relationship is monotonically negative...\n")
  if (all(pred_df$fit < 0)) {
    cat("  All predicted AP_change values are negative.\n")
    cat("  This suggests Cerealia is associated with AP decline at all levels.\n")
  } else if (all(pred_df$fit > 0)) {
    cat("  All predicted AP_change values are positive.\n")
  } else {
    cat("  Mixed signs but no clean positive-to-negative transition.\n")
  }
}

pred_file <- file.path(OUTPUT_DIR, "threshold_gam_predictions.csv")
write.csv(pred_df, pred_file, row.names = FALSE)

# --- Segmented regression ---

cat("\n=== Segmented Regression ===\n")

seg_threshold <- NA_real_
seg_ci_lower  <- NA_real_
seg_ci_upper  <- NA_real_

if (HAS_SEGMENTED) {
  tryCatch({
    lm_fit <- lm(ap_change ~ cerealia_level + region, data = pooled_data)

    seg_fit <- segmented(lm_fit,
                         seg.Z = ~ cerealia_level,
                         psi = median(pooled_data$cerealia_level))

    cat("\nSegmented regression summary:\n")
    print(summary(seg_fit))

    bp <- seg_fit$psi
    seg_threshold <- bp[1, "Est."]
    seg_se <- bp[1, "St.Err"]
    seg_ci_lower <- seg_threshold - 1.96 * seg_se
    seg_ci_upper <- seg_threshold + 1.96 * seg_se

    cat("\nBreakpoint:", round(seg_threshold, 4), "%\n")
    cat("  95% CI: [", round(seg_ci_lower, 4), ",",
        round(seg_ci_upper, 4), "] %\n")

  }, error = function(e) {
    cat("Segmented regression failed:", e$message, "\n")
  })
} else {
  cat("Package 'segmented' not available. Skipping.\n")
}

# --- Output ---

cat("\n\n========================================\n")
cat("THRESHOLD ANALYSIS RESULTS\n")
cat("========================================\n\n")

threshold_results <- data.frame(
  method = c("GAM smooth zero-crossing", "Segmented regression"),
  threshold_cerealia_pct = c(
    ifelse(is.na(gam_threshold), NA, round(gam_threshold, 4)),
    ifelse(is.na(seg_threshold), NA, round(seg_threshold, 4))
  ),
  ci_lower = c(
    ifelse(is.na(gam_threshold_ci_lower), NA,
           round(min(gam_threshold_ci_lower, gam_threshold_ci_upper, na.rm = TRUE), 4)),
    ifelse(is.na(seg_ci_lower), NA, round(seg_ci_lower, 4))
  ),
  ci_upper = c(
    ifelse(is.na(gam_threshold_ci_upper), NA,
           round(max(gam_threshold_ci_lower, gam_threshold_ci_upper, na.rm = TRUE), 4)),
    ifelse(is.na(seg_ci_upper), NA, round(seg_ci_upper, 4))
  ),
  stringsAsFactors = FALSE
)

print(threshold_results, row.names = FALSE)

results_file <- file.path(OUTPUT_DIR, "threshold_gam_results.csv")
write.csv(threshold_results, results_file, row.names = FALSE)
cat("\nResults saved to:", results_file, "\n")

# Text summary
summary_file <- file.path(OUTPUT_DIR, "threshold_gam_summary.txt")
sink(summary_file)
cat("THRESHOLD / NONLINEARITY ANALYSIS SUMMARY\n")
cat("==========================================\n\n")
cat("Data: Pooled 200-yr bins from Britain, Scandinavia, Alps\n")
cat("Note: Plantago proxy = pastoral_pct from cached data\n")
cat("Total observations:", nrow(pooled_data), "\n\n")
cat("--- GAM ---\n")
cat("Model: AP_change ~ s(Cerealia) + s(Pastoral) + region\n")
print(summary(gam_fit))
cat("\nCerealia threshold (zero-crossing):",
    ifelse(is.na(gam_threshold), "Not found",
           paste0(round(gam_threshold, 4), "%")), "\n")
if (HAS_SEGMENTED && !is.na(seg_threshold)) {
  cat("\n--- Segmented Regression ---\n")
  cat("Breakpoint:", round(seg_threshold, 4), "%\n")
  cat("95% CI: [", round(seg_ci_lower, 4), ",", round(seg_ci_upper, 4), "]\n")
}
sink()
cat("Summary saved to:", summary_file, "\n")

cat("\n\n################################################################\n")
cat("# ALL ANALYSES COMPLETE\n")
cat("################################################################\n")
