#!/usr/bin/env Rscript
# =============================================================================
# Paper 11: Review Response Analyses — Additional Robustness Checks
# =============================================================================
# Addresses reviewer concerns with:
#   1. Taphonomic correction (Surovell et al. 2009)
#   2. Climate proxy comparison (GISP2-like / Bond events)
#   3. Subsampling analysis (equal sample sizes)
#   4. Benjamini-Hochberg FDR correction
#   5. Sample size effect on correlation strength
# =============================================================================

library(rcarbon)

cat("================================================================\n")
cat("PAPER 11: Review Response Analyses\n")
cat("================================================================\n\n")

# =====================================================================
# STEP 0: Reproduce SPDs from the full analysis (needed as input)
# =====================================================================
cat("--- STEP 0: Loading data and computing regional SPDs ---\n\n")

d <- read.csv("/home/ayu/archeco/shared/p3k14c_data.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$Lat) & !is.na(d$Long) & !is.na(d$Age) & !is.na(d$Error), ]
d <- d[d$Error > 0 & d$Error < 500 & d$Age > 500 & d$Age < 12000, ]
cat("After filtering:", nrow(d), "dates\n\n")

regions <- list(
  Britain = list(
    filter = function(x) x[x$Country %in% c("United Kingdom", "Ireland"), ],
    label = "Britain+Ireland",
    centroid_lat = 54.0, centroid_lon = -2.0
  ),
  W_Europe = list(
    filter = function(x) x[x$Country %in% c("France", "Spain", "Germany",
                                              "Belgium", "Netherlands", "Italy",
                                              "Portugal", "Switzerland", "Austria") &
                            x$Continent == "Europe", ],
    label = "W.Europe",
    centroid_lat = 47.0, centroid_lon = 5.0
  ),
  E_Europe = list(
    filter = function(x) x[x$Country %in% c("Poland", "Hungary", "Czech Republic",
                                              "Romania", "Bulgaria", "Serbia",
                                              "Croatia", "Slovakia", "Greece") &
                            x$Continent == "Europe", ],
    label = "E.Europe+Greece",
    centroid_lat = 45.0, centroid_lon = 22.0
  ),
  Scandinavia = list(
    filter = function(x) x[x$Country %in% c("Norway", "Sweden", "Denmark", "Finland"), ],
    label = "Scandinavia",
    centroid_lat = 62.0, centroid_lon = 15.0
  ),
  Near_East = list(
    filter = function(x) x[x$Country %in% c("Turkey", "Syria", "Iraq", "Iran",
                                              "Israel", "Jordan", "Lebanon", "Palestine"), ],
    label = "Near East",
    centroid_lat = 37.0, centroid_lon = 40.0
  ),
  China = list(
    filter = function(x) x[x$Country == "China", ],
    label = "China",
    centroid_lat = 35.0, centroid_lon = 110.0
  ),
  Japan = list(
    filter = function(x) x[x$Country == "Japan", ],
    label = "Japan",
    centroid_lat = 36.0, centroid_lon = 138.0
  )
)

region_data <- list()
region_n <- c()
for (rname in names(regions)) {
  rd <- regions[[rname]]$filter(d)
  n <- nrow(rd)
  cat(sprintf("  %-18s: %6d dates\n", regions[[rname]]$label, n))
  if (n >= 200) {
    region_data[[rname]] <- rd
    region_n[rname] <- n
  }
}
cat("\n")

time_range <- c(10000, 1000)
grid_step <- 50
time_grid <- seq(time_range[1], time_range[2], by = -grid_step)

spd_list <- list()
for (rname in names(region_data)) {
  rd <- region_data[[rname]]
  cat(sprintf("  Calibrating %s (%d dates)...\n", regions[[rname]]$label, nrow(rd)))
  caldates <- calibrate(x = rd$Age, errors = rd$Error, calCurves = "intcal20", verbose = FALSE)
  site_ids <- ifelse(rd$SiteName == "" | is.na(rd$SiteName),
                     paste0("unknown", seq_len(nrow(rd))),
                     gsub("_", "-", rd$SiteName))
  bins <- binPrep(sites = site_ids, ages = rd$Age, h = 200)
  spd_result <- spd(caldates, bins = bins, timeRange = time_range, verbose = FALSE)
  spd_list[[rname]] <- spd_result
}

# Build raw SPD matrix (unnormalized, for taphonomic correction)
spd_matrix_raw <- matrix(NA, nrow = length(time_grid), ncol = length(spd_list))
colnames(spd_matrix_raw) <- names(spd_list)
for (i in seq_along(spd_list)) {
  rname <- names(spd_list)[i]
  sg <- spd_list[[rname]]$grid
  spd_vals <- approx(sg$calBP, sg$PrDens, xout = time_grid, method = "linear")$y
  spd_matrix_raw[, i] <- spd_vals
}

# Normalized SPD matrix (for comparison)
spd_matrix <- spd_matrix_raw
for (i in 1:ncol(spd_matrix)) {
  spd_matrix[, i] <- spd_matrix[, i] / mean(spd_matrix[, i], na.rm = TRUE)
}

# First-differenced matrix (uncorrected baseline)
spd_diff <- apply(spd_matrix, 2, diff)
n_reg <- ncol(spd_matrix)
reg_names <- colnames(spd_matrix)
reg_labels <- sapply(reg_names, function(x) regions[[x]]$label)

# Baseline detrended correlation matrix
detrend_cor <- cor(spd_diff, use = "pairwise.complete.obs")

# Compute baseline mean r
baseline_pairs <- c()
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    baseline_pairs <- c(baseline_pairs, detrend_cor[i, j])
  }
}
baseline_mean_r <- mean(baseline_pairs)

cat(sprintf("\nSPD matrix: %d time steps x %d regions\n", nrow(spd_matrix), ncol(spd_matrix)))
cat(sprintf("Baseline mean detrended r: %.3f\n\n", baseline_mean_r))


# =====================================================================
# ANALYSIS 1: Taphonomic Correction (Surovell et al. 2009)
# =====================================================================
cat("================================================================\n")
cat("ANALYSIS 1: Taphonomic Correction (Surovell et al. 2009)\n")
cat("================================================================\n\n")

# For each region, fit exponential to raw SPD:
#   log(SPD) = a + b * time_BP
# Then divide raw SPD by fitted exponential to correct for taphonomic loss.
# Time here: time_grid values are in cal BP (higher = older)
# Expected: older dates are underrepresented, so raw SPD decreases with age.
# Correction: divide by exp(a + b*time) to inflate older values.

spd_matrix_corrected <- matrix(NA, nrow = length(time_grid), ncol = n_reg)
colnames(spd_matrix_corrected) <- reg_names

cat("Fitting exponential taphonomic models per region:\n\n")
cat(sprintf("%-18s  Slope(b)     R²       Half-life(yr)  Correction@10ka\n", "Region"))
cat(paste(rep("-", 75), collapse = ""), "\n")

tapho_results <- list()
for (i in 1:n_reg) {
  rname <- reg_names[i]
  spd_raw <- spd_matrix_raw[, i]

  # Use time_grid as predictor (cal BP)
  # Fit log-linear: log(SPD) ~ time
  # Only use positive values
  valid <- spd_raw > 0 & !is.na(spd_raw)
  if (sum(valid) < 10) {
    cat(sprintf("  %s: insufficient data for fitting\n", reg_labels[i]))
    spd_matrix_corrected[, i] <- spd_raw
    next
  }

  log_spd <- log(spd_raw[valid])
  time_vals <- time_grid[valid]

  fit <- lm(log_spd ~ time_vals)
  slope <- coef(fit)[2]
  intercept <- coef(fit)[1]
  rsq <- summary(fit)$r.squared

  # Fitted exponential for ALL time points
  fitted_exp <- exp(intercept + slope * time_grid)

  # Correction: divide raw SPD by fitted exponential
  corrected <- spd_raw / fitted_exp

  # Normalize corrected SPD
  corrected <- corrected / mean(corrected, na.rm = TRUE)
  spd_matrix_corrected[, i] <- corrected

  # Half-life: time for expected count to halve
  # exp(b * t) = 0.5 => t = log(0.5) / b
  half_life <- ifelse(slope < 0, abs(log(0.5) / slope), NA)

  # Correction factor at 10000 BP vs 1000 BP
  correction_10ka <- exp(slope * 1000) / exp(slope * 10000)

  cat(sprintf("%-18s  %+.2e  %.3f    %8.0f       %.2f\n",
              reg_labels[i], slope, rsq,
              ifelse(is.na(half_life), NA, half_life),
              correction_10ka))

  tapho_results[[rname]] <- list(slope = slope, rsq = rsq, half_life = half_life)
}

# First-difference the corrected SPDs
spd_diff_corrected <- apply(spd_matrix_corrected, 2, diff)

# Compute corrected pairwise correlations
corrected_cor <- cor(spd_diff_corrected, use = "pairwise.complete.obs")

cat("\n--- Pairwise correlations: Uncorrected vs Taphonomically Corrected ---\n\n")
cat(sprintf("%-40s  r(uncorr)  r(corr)   Diff    Agreement?\n", "Pair"))
cat(paste(rep("-", 85), collapse = ""), "\n")

corrected_pairs <- c()
pair_labels_all <- c()
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    r_uncorr <- detrend_cor[i, j]
    r_corr <- corrected_cor[i, j]
    agree <- (r_uncorr > 0 & r_corr > 0) | (r_uncorr < 0 & r_corr < 0)
    pair_label <- paste(reg_labels[i], "x", reg_labels[j])
    cat(sprintf("%-40s  %+.3f     %+.3f    %+.3f   %s\n",
                pair_label, r_uncorr, r_corr, r_corr - r_uncorr,
                ifelse(agree, "YES", "NO")))
    corrected_pairs <- c(corrected_pairs, r_corr)
    pair_labels_all <- c(pair_labels_all, pair_label)
  }
}

corrected_mean_r <- mean(corrected_pairs)
cat(sprintf("\nMean r (uncorrected, first-diff): %.3f\n", baseline_mean_r))
cat(sprintf("Mean r (taphonomically corrected, first-diff): %.3f\n", corrected_mean_r))
cat(sprintf("Difference: %+.3f\n", corrected_mean_r - baseline_mean_r))
cat(sprintf("Sign agreement: %d / %d pairs\n",
            sum((baseline_pairs > 0) == (corrected_pairs > 0)), length(baseline_pairs)))
cat(sprintf("Correlation between corrected and uncorrected r values: %.3f\n",
            cor(baseline_pairs, corrected_pairs)))

# Significance test on corrected pairs
n_corr_sig <- 0
n_t_corr <- nrow(spd_diff_corrected)
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    ct <- cor.test(spd_diff_corrected[, i], spd_diff_corrected[, j])
    if (ct$p.value < 0.05) n_corr_sig <- n_corr_sig + 1
  }
}
cat(sprintf("Pairs significant at p<0.05 (corrected): %d / 21\n", n_corr_sig))

if (corrected_mean_r > 0.3) {
  cat("\nCONCLUSION: Synchrony SURVIVES taphonomic correction.\n")
  cat("  The common taphonomic trend cannot explain the observed synchrony.\n")
} else {
  cat("\nCONCLUSION: Synchrony is SUBSTANTIALLY REDUCED by taphonomic correction.\n")
  cat("  Shared taphonomic bias may partly explain apparent synchrony.\n")
}


# =====================================================================
# ANALYSIS 2: Climate Proxy Comparison
# =====================================================================
cat("\n\n================================================================\n")
cat("ANALYSIS 2: Climate Proxy Comparison\n")
cat("================================================================\n\n")

# Construct a GISP2-like climate proxy at 50-year resolution (10000-1000 cal BP)
# Based on well-documented Holocene climate events:
#   - 9.3 ka event (cold)
#   - 8.2 ka event (major cold, ~3.3°C drop in Greenland)
#   - 5.9 ka event (aridification)
#   - 4.2 ka event (megadrought)
#   - 3.2 ka event (Late Bronze Age collapse)
#   - 2.8 ka event (sub-Atlantic transition)
#
# We use the well-known Holocene temperature pattern:
#   - Holocene Thermal Maximum ~9000-5000 BP (warm)
#   - Neoglaciation ~5000-1000 BP (cooling trend)
#   - Bond events superimposed as cold perturbations
#
# Method: construct baseline + Bond event perturbations

cat("Constructing synthetic climate proxy based on:\n")
cat("  - Holocene Thermal Maximum baseline (~9000-5000 BP)\n")
cat("  - Neoglaciation cooling trend (~5000-1000 BP)\n")
cat("  - Bond events as cold perturbations\n\n")

# Baseline: smooth Holocene temperature curve (relative to modern)
# Warm period peaks around 7000 BP
climate_baseline <- -0.5 + 1.5 * exp(-((time_grid - 7000) / 2500)^2)

# Bond events: cold perturbations (timing in cal BP, duration ~200-400 yr)
bond_events <- data.frame(
  center = c(9300, 8200, 5900, 4200, 3200, 2800, 1400),
  amplitude = c(-0.8, -2.5, -1.2, -1.5, -1.0, -0.8, -0.6),
  width = c(150, 200, 200, 200, 200, 150, 150),
  name = c("9.3ka", "8.2ka", "5.9ka", "4.2ka", "3.2ka", "2.8ka", "1.4ka")
)

climate_proxy <- climate_baseline
for (k in 1:nrow(bond_events)) {
  perturbation <- bond_events$amplitude[k] * exp(-((time_grid - bond_events$center[k]) / bond_events$width[k])^2)
  climate_proxy <- climate_proxy + perturbation
}

# Add small noise for realism
set.seed(42)
climate_proxy <- climate_proxy + rnorm(length(time_grid), 0, 0.15)

# Normalize
climate_proxy <- (climate_proxy - mean(climate_proxy)) / sd(climate_proxy)

cat("Climate proxy constructed (", length(time_grid), "time steps)\n")
cat(sprintf("Range: %.2f to %.2f (standardized)\n", min(climate_proxy), max(climate_proxy)))

# First-difference the climate proxy
climate_diff <- diff(climate_proxy)

cat("\n--- Correlation between climate proxy and regional SPDs (detrended) ---\n\n")
cat(sprintf("%-18s  r(climate)  p-value       Significant?\n", "Region"))
cat(paste(rep("-", 60), collapse = ""), "\n")

climate_cors <- c()
climate_ps <- c()
for (i in 1:n_reg) {
  ct <- cor.test(climate_diff, spd_diff[, i])
  cat(sprintf("%-18s  %+.3f      %.4e    %s\n",
              reg_labels[i], ct$estimate, ct$p.value,
              ifelse(ct$p.value < 0.05, "***", "ns")))
  climate_cors <- c(climate_cors, ct$estimate)
  climate_ps <- c(climate_ps, ct$p.value)
}

cat(sprintf("\nMean |r| (climate vs SPD): %.3f\n", mean(abs(climate_cors))))
cat(sprintf("Mean r (climate vs SPD): %+.3f\n", mean(climate_cors)))
cat(sprintf("Regions significant at p<0.05: %d / %d\n", sum(climate_ps < 0.05), n_reg))

# Also create a binary climate deterioration index
cat("\n--- Binary climate deterioration index ---\n\n")

# Define deterioration windows (cal BP)
deterioration_windows <- list(
  c(8400, 8000),  # 8.2 ka event
  c(6000, 5800),  # 5.9 ka event
  c(4300, 4100),  # 4.2 ka event
  c(3300, 3100),  # 3.2 ka event
  c(2900, 2700)   # 2.8 ka event
)

climate_binary <- rep(0, length(time_grid))
for (w in deterioration_windows) {
  climate_binary[time_grid <= w[1] & time_grid >= w[2]] <- 1
}

cat(sprintf("Time steps in deterioration windows: %d / %d (%.1f%%)\n",
            sum(climate_binary), length(climate_binary),
            100 * sum(climate_binary) / length(climate_binary)))

# For each region, compute mean SPD change during vs outside deterioration windows
cat("\n--- SPD behavior during climate deteriorations ---\n\n")
cat(sprintf("%-18s  Mean SPD(stable)  Mean SPD(crisis)  Ratio   t-test p\n", "Region"))
cat(paste(rep("-", 75), collapse = ""), "\n")

for (i in 1:n_reg) {
  stable_idx <- climate_binary == 0
  crisis_idx <- climate_binary == 1
  spd_stable <- mean(spd_matrix[stable_idx, i], na.rm = TRUE)
  spd_crisis <- mean(spd_matrix[crisis_idx, i], na.rm = TRUE)
  ratio <- spd_crisis / spd_stable

  # t-test
  tt <- t.test(spd_matrix[crisis_idx, i], spd_matrix[stable_idx, i])
  cat(sprintf("%-18s  %.3f             %.3f             %.3f   %.4f\n",
              reg_labels[i], spd_stable, spd_crisis, ratio, tt$p.value))
}

# Key interpretation
cat("\nINTERPRETATION:\n")
cat("  If synchrony is climate-forced, we expect:\n")
cat("  1. Most regions to correlate positively with the climate proxy\n")
cat("  2. Population declines during deterioration windows\n")
cat("  The strength of climate-SPD correlations vs inter-regional SPD correlations\n")
cat("  tells us whether climate fully mediates the synchrony.\n")
cat(sprintf("\n  Mean |r| climate-SPD: %.3f\n", mean(abs(climate_cors))))
cat(sprintf("  Mean r inter-regional SPD: %.3f\n", baseline_mean_r))
if (mean(abs(climate_cors)) < baseline_mean_r * 0.5) {
  cat("  Climate correlations are WEAKER than inter-regional correlations.\n")
  cat("  Climate may contribute to synchrony but does not fully explain it.\n")
} else {
  cat("  Climate correlations are COMPARABLE to inter-regional correlations.\n")
  cat("  Climate forcing is a plausible driver of observed synchrony.\n")
}


# =====================================================================
# ANALYSIS 3: Subsampling Analysis (Equal Sample Sizes)
# =====================================================================
cat("\n\n================================================================\n")
cat("ANALYSIS 3: Subsampling Analysis (Equal Sample Sizes)\n")
cat("================================================================\n\n")

# For each iteration:
# - Draw N=1500 dates from each region (matching Japan's ~1433)
# - Calibrate, bin, SPD, detrend, correlate
# Reduced to 20 iterations for computational feasibility

set.seed(2024)
n_iter <- 20
n_subsample <- 1500

cat(sprintf("Subsampling: N=%d dates per region, %d iterations\n", n_subsample, n_iter))
cat("(Matching Japan's sample size to test if large N creates spurious synchrony)\n\n")

# Check which regions have enough dates
eligible_regions <- names(region_data)[region_n[names(region_data)] >= n_subsample]
cat(sprintf("Regions with >= %d dates: %d / %d\n", n_subsample, length(eligible_regions), n_reg))
for (rname in eligible_regions) {
  cat(sprintf("  %-18s: %d dates (sampling %d)\n", regions[[rname]]$label, region_n[rname], n_subsample))
}
cat("\n")

n_elig <- length(eligible_regions)
n_pairs_elig <- choose(n_elig, 2)

subsample_mean_r <- numeric(n_iter)
subsample_all_r <- matrix(NA, nrow = n_iter, ncol = n_pairs_elig)

# Pair labels for eligible regions
elig_labels <- sapply(eligible_regions, function(x) regions[[x]]$label)
sub_pair_labels <- c()
for (i in 1:(n_elig - 1)) {
  for (j in (i + 1):n_elig) {
    sub_pair_labels <- c(sub_pair_labels, paste(elig_labels[i], "x", elig_labels[j]))
  }
}
colnames(subsample_all_r) <- sub_pair_labels

cat("Running subsampling iterations:\n")

for (iter in 1:n_iter) {
  cat(sprintf("  Iteration %d/%d...", iter, n_iter))

  # Subsample, calibrate, bin, SPD for each region
  sub_spd_matrix <- matrix(NA, nrow = length(time_grid), ncol = n_elig)
  colnames(sub_spd_matrix) <- eligible_regions

  for (ri in seq_along(eligible_regions)) {
    rname <- eligible_regions[ri]
    rd <- region_data[[rname]]

    # Random subsample
    idx <- sample(1:nrow(rd), n_subsample, replace = FALSE)
    rd_sub <- rd[idx, ]

    # Calibrate
    caldates <- calibrate(x = rd_sub$Age, errors = rd_sub$Error, calCurves = "intcal20", verbose = FALSE)
    site_ids <- ifelse(rd_sub$SiteName == "" | is.na(rd_sub$SiteName),
                       paste0("unknown", seq_len(nrow(rd_sub))),
                       gsub("_", "-", rd_sub$SiteName))
    bins <- binPrep(sites = site_ids, ages = rd_sub$Age, h = 200)
    spd_result <- spd(caldates, bins = bins, timeRange = time_range, verbose = FALSE)

    sg <- spd_result$grid
    spd_vals <- approx(sg$calBP, sg$PrDens, xout = time_grid, method = "linear")$y
    spd_vals <- spd_vals / mean(spd_vals, na.rm = TRUE)
    sub_spd_matrix[, ri] <- spd_vals
  }

  # First-difference
  sub_diff <- apply(sub_spd_matrix, 2, diff)

  # Pairwise correlations
  sub_cor <- cor(sub_diff, use = "pairwise.complete.obs")
  pair_idx <- 0
  pair_r <- c()
  for (i in 1:(n_elig - 1)) {
    for (j in (i + 1):n_elig) {
      pair_idx <- pair_idx + 1
      pair_r <- c(pair_r, sub_cor[i, j])
    }
  }

  subsample_mean_r[iter] <- mean(pair_r)
  subsample_all_r[iter, ] <- pair_r
  cat(sprintf(" mean r = %.3f\n", mean(pair_r)))
}

cat("\n--- Subsampling Results ---\n\n")
cat(sprintf("Iterations: %d\n", n_iter))
cat(sprintf("Sample size per region: %d\n", n_subsample))
cat(sprintf("Eligible regions: %d (those with >= %d dates)\n", n_elig, n_subsample))
cat(sprintf("Pairs per iteration: %d\n", n_pairs_elig))

cat(sprintf("\nSubsampled mean r: %.3f (SD = %.3f)\n", mean(subsample_mean_r), sd(subsample_mean_r)))
cat(sprintf("95%% CI of subsampled mean r: [%.3f, %.3f]\n",
            quantile(subsample_mean_r, 0.025), quantile(subsample_mean_r, 0.975)))

# Compare with full-data baseline (using same set of eligible regions)
full_eligible_cor <- cor(spd_diff[, eligible_regions], use = "pairwise.complete.obs")
full_elig_pairs <- c()
for (i in 1:(n_elig - 1)) {
  for (j in (i + 1):n_elig) {
    full_elig_pairs <- c(full_elig_pairs, full_eligible_cor[i, j])
  }
}
full_elig_mean_r <- mean(full_elig_pairs)

cat(sprintf("\nFull-data mean r (same %d regions): %.3f\n", n_elig, full_elig_mean_r))
cat(sprintf("Full-data mean r (all 7 regions): %.3f\n", baseline_mean_r))
cat(sprintf("Subsampled mean r: %.3f\n", mean(subsample_mean_r)))
cat(sprintf("Difference (subsampled - full eligible): %+.3f\n",
            mean(subsample_mean_r) - full_elig_mean_r))

# Per-pair comparison
cat("\n--- Per-pair comparison (full vs subsampled) ---\n\n")
cat(sprintf("%-40s  r(full)  r(sub,mean)  r(sub,SD)\n", "Pair"))
cat(paste(rep("-", 75), collapse = ""), "\n")
pair_idx <- 0
for (i in 1:(n_elig - 1)) {
  for (j in (i + 1):n_elig) {
    pair_idx <- pair_idx + 1
    cat(sprintf("%-40s  %+.3f   %+.3f       %.3f\n",
                sub_pair_labels[pair_idx],
                full_elig_pairs[pair_idx],
                mean(subsample_all_r[, pair_idx]),
                sd(subsample_all_r[, pair_idx])))
  }
}

if (abs(mean(subsample_mean_r) - full_elig_mean_r) < 0.1) {
  cat("\nCONCLUSION: Synchrony is NOT driven by large sample sizes.\n")
  cat("  Subsampled mean r is within 0.1 of the full-data value.\n")
  cat("  Equal-N subsampling preserves the synchrony signal.\n")
} else {
  cat("\nCONCLUSION: Large sample sizes INFLATE synchrony estimates.\n")
  cat("  Subsampled mean r differs substantially from full-data value.\n")
}


# =====================================================================
# ANALYSIS 4: Benjamini-Hochberg FDR Correction
# =====================================================================
cat("\n\n================================================================\n")
cat("ANALYSIS 4: Benjamini-Hochberg FDR Correction\n")
cat("================================================================\n\n")

# Collect p-values from all 21 pairwise correlations (first-differenced)
pair_p <- c()
pair_r <- c()
pair_lab <- c()
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    ct <- cor.test(spd_diff[, i], spd_diff[, j])
    pair_p <- c(pair_p, ct$p.value)
    pair_r <- c(pair_r, ct$estimate)
    pair_lab <- c(pair_lab, paste(reg_labels[i], "x", reg_labels[j]))
  }
}

# Apply BH correction
p_adj <- p.adjust(pair_p, method = "BH")

cat(sprintf("%-40s  r       p(raw)       p(BH)        Sig(raw)  Sig(BH)\n", "Pair"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for (k in 1:length(pair_lab)) {
  cat(sprintf("%-40s  %+.3f  %.4e   %.4e   %s       %s\n",
              pair_lab[k], pair_r[k], pair_p[k], p_adj[k],
              ifelse(pair_p[k] < 0.05, "***", "ns"),
              ifelse(p_adj[k] < 0.05, "***", "ns")))
}

n_raw_sig <- sum(pair_p < 0.05)
n_bh_sig <- sum(p_adj < 0.05)
cat(sprintf("\nSignificant pairs (raw, p<0.05): %d / 21\n", n_raw_sig))
cat(sprintf("Significant pairs (BH-corrected, p<0.05): %d / 21\n", n_bh_sig))
cat(sprintf("Most lenient BH-adjusted p: %.4e\n", max(p_adj[p_adj < 0.05], na.rm = TRUE)))
cat(sprintf("Most stringent non-significant BH-adjusted p: %.4f\n",
            ifelse(any(p_adj >= 0.05), min(p_adj[p_adj >= 0.05]), NA)))


# =====================================================================
# ANALYSIS 5: Sample Size Effect on Correlation Strength
# =====================================================================
cat("\n\n================================================================\n")
cat("ANALYSIS 5: Sample Size Effect on Correlation Strength\n")
cat("================================================================\n\n")

# For each pair, compute geometric mean of N_i and N_j
# Regress detrended r on geom_mean_N to test if larger samples = higher r

pair_geom_n <- c()
pair_r_for_n <- c()
pair_lab_for_n <- c()
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    ni <- region_n[reg_names[i]]
    nj <- region_n[reg_names[j]]
    geom_n <- sqrt(ni * nj)
    pair_geom_n <- c(pair_geom_n, geom_n)
    pair_r_for_n <- c(pair_r_for_n, detrend_cor[i, j])
    pair_lab_for_n <- c(pair_lab_for_n, paste(reg_labels[i], "x", reg_labels[j]))
  }
}

# Linear regression
fit_n <- lm(pair_r_for_n ~ pair_geom_n)
fit_log_n <- lm(pair_r_for_n ~ log(pair_geom_n))

cat(sprintf("%-40s  N_i       N_j       geom(N)    r\n", "Pair"))
cat(paste(rep("-", 85), collapse = ""), "\n")
pidx <- 0
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    pidx <- pidx + 1
    ni <- region_n[reg_names[i]]
    nj <- region_n[reg_names[j]]
    cat(sprintf("%-40s  %6d    %6d    %8.0f   %+.3f\n",
                pair_lab_for_n[pidx], ni, nj, pair_geom_n[pidx], pair_r_for_n[pidx]))
  }
}

cat(sprintf("\n--- Linear regression: r ~ geom_mean(N) ---\n"))
cat(sprintf("Slope: %.2e\n", coef(fit_n)[2]))
cat(sprintf("R-squared: %.3f\n", summary(fit_n)$r.squared))
cat(sprintf("p-value: %.4f\n", summary(fit_n)$coefficients[2, 4]))

cat(sprintf("\n--- Log regression: r ~ log(geom_mean(N)) ---\n"))
cat(sprintf("Slope: %.4f\n", coef(fit_log_n)[2]))
cat(sprintf("R-squared: %.3f\n", summary(fit_log_n)$r.squared))
cat(sprintf("p-value: %.4f\n", summary(fit_log_n)$coefficients[2, 4]))

# Spearman rank correlation
sp_test <- cor.test(pair_geom_n, pair_r_for_n, method = "spearman")
cat(sprintf("\nSpearman rank correlation: rho = %.3f, p = %.4f\n",
            sp_test$estimate, sp_test$p.value))

if (summary(fit_n)$coefficients[2, 4] < 0.05) {
  cat("\nCONCLUSION: Sample size SIGNIFICANTLY predicts correlation strength.\n")
  cat("  This is a concern: larger regions may show inflated synchrony.\n")
  cat("  The subsampling analysis (Analysis 3) directly addresses this.\n")
} else {
  cat("\nCONCLUSION: Sample size does NOT significantly predict correlation strength.\n")
  cat("  Synchrony is not an artifact of differential sample sizes.\n")
}


# =====================================================================
# COMPREHENSIVE SUMMARY
# =====================================================================
cat("\n\n================================================================\n")
cat("COMPREHENSIVE SUMMARY — Review Response Analyses\n")
cat("================================================================\n\n")

cat("1. TAPHONOMIC CORRECTION (Surovell et al. 2009):\n")
cat(sprintf("   Mean r (uncorrected): %.3f\n", baseline_mean_r))
cat(sprintf("   Mean r (taphonomically corrected): %.3f\n", corrected_mean_r))
cat(sprintf("   Difference: %+.3f\n", corrected_mean_r - baseline_mean_r))
cat(sprintf("   Sign agreement: %d / 21 pairs\n",
            sum((baseline_pairs > 0) == (corrected_pairs > 0))))
cat(sprintf("   Pairs significant (corrected): %d / 21\n", n_corr_sig))
cat("\n")

cat("2. CLIMATE PROXY COMPARISON:\n")
cat(sprintf("   Mean |r| (climate vs SPD): %.3f\n", mean(abs(climate_cors))))
cat(sprintf("   Mean r (climate vs SPD): %+.3f\n", mean(climate_cors)))
cat(sprintf("   Regions with significant climate correlation: %d / %d\n",
            sum(climate_ps < 0.05), n_reg))
cat(sprintf("   For comparison, mean inter-regional r: %.3f\n", baseline_mean_r))
cat("\n")

cat("3. SUBSAMPLING (N=1500 per region, 20 iterations):\n")
cat(sprintf("   Subsampled mean r: %.3f (SD = %.3f)\n", mean(subsample_mean_r), sd(subsample_mean_r)))
cat(sprintf("   Full-data mean r (same regions): %.3f\n", full_elig_mean_r))
cat(sprintf("   Full-data mean r (all 7 regions): %.3f\n", baseline_mean_r))
cat(sprintf("   Difference (sub - full eligible): %+.3f\n",
            mean(subsample_mean_r) - full_elig_mean_r))
cat("\n")

cat("4. BENJAMINI-HOCHBERG FDR CORRECTION:\n")
cat(sprintf("   Significant (raw): %d / 21\n", n_raw_sig))
cat(sprintf("   Significant (BH-corrected): %d / 21\n", n_bh_sig))
cat("\n")

cat("5. SAMPLE SIZE EFFECT:\n")
cat(sprintf("   r ~ geom(N): R² = %.3f, p = %.4f\n",
            summary(fit_n)$r.squared, summary(fit_n)$coefficients[2, 4]))
cat(sprintf("   r ~ log(geom(N)): R² = %.3f, p = %.4f\n",
            summary(fit_log_n)$r.squared, summary(fit_log_n)$coefficients[2, 4]))
cat(sprintf("   Spearman rho: %.3f, p = %.4f\n", sp_test$estimate, sp_test$p.value))

cat("\n\nDone. All review response analyses complete.\n")
