#!/usr/bin/env Rscript
# =============================================================================
# Paper 11: Full Analysis — Global Demographic Synchrony
# =============================================================================
# Builds on pilot results. Runs:
#   1. Distance-decay of synchrony (Mantel test)
#   2. Block bootstrap confidence intervals (21 pairs)
#   3. Effective degrees of freedom correction (Bartlett)
#   4. Alternative detrending (loess residuals, span=0.3)
#   5. Sub-regional analysis within W.Europe (France, Iberia, Central)
# =============================================================================

library(rcarbon)

cat("================================================================\n")
cat("PAPER 11: Full Analysis — Global Demographic Synchrony\n")
cat("================================================================\n\n")

# =====================================================================
# STEP 0: Reproduce pilot SPDs (needed as input for all analyses)
# =====================================================================
cat("--- STEP 0: Loading data and computing regional SPDs ---\n\n")

d <- read.csv("/home/ayu/archeco/shared/p3k14c_data.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$Lat) & !is.na(d$Long) & !is.na(d$Age) & !is.na(d$Error), ]
d <- d[d$Error > 0 & d$Error < 500 & d$Age > 500 & d$Age < 12000, ]
cat("After filtering:", nrow(d), "dates\n\n")

# Region definitions (same as pilot)
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

# Extract region data
region_data <- list()
for (rname in names(regions)) {
  rd <- regions[[rname]]$filter(d)
  n <- nrow(rd)
  cat(sprintf("  %-18s: %6d dates\n", regions[[rname]]$label, n))
  if (n >= 200) region_data[[rname]] <- rd
}
cat("\n")

# Calibrate and compute SPDs
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

# Build SPD matrix
spd_matrix <- matrix(NA, nrow = length(time_grid), ncol = length(spd_list))
colnames(spd_matrix) <- names(spd_list)
for (i in seq_along(spd_list)) {
  rname <- names(spd_list)[i]
  sg <- spd_list[[rname]]$grid
  spd_vals <- approx(sg$calBP, sg$PrDens, xout = time_grid, method = "linear")$y
  spd_vals <- spd_vals / mean(spd_vals, na.rm = TRUE)
  spd_matrix[, i] <- spd_vals
}

# First-differenced matrix
spd_diff <- apply(spd_matrix, 2, diff)
n_reg <- ncol(spd_matrix)
reg_names <- colnames(spd_matrix)
reg_labels <- sapply(reg_names, function(x) regions[[x]]$label)

cat("\nSPD matrix:", nrow(spd_matrix), "time steps x", ncol(spd_matrix), "regions\n")
cat("Detrended matrix:", nrow(spd_diff), "time steps\n\n")

# Store all results for output
results <- list()

# =====================================================================
# ANALYSIS 1: Distance-decay of synchrony (Mantel test)
# =====================================================================
cat("================================================================\n")
cat("ANALYSIS 1: Distance-Decay of Synchrony (Mantel Test)\n")
cat("================================================================\n\n")

# Haversine distance between centroids (km)
haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371  # Earth radius in km
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat / 2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R * c
}

# Build distance matrix
dist_matrix <- matrix(0, nrow = n_reg, ncol = n_reg)
rownames(dist_matrix) <- colnames(dist_matrix) <- reg_names
for (i in 1:n_reg) {
  for (j in 1:n_reg) {
    dist_matrix[i, j] <- haversine(
      regions[[reg_names[i]]]$centroid_lat, regions[[reg_names[i]]]$centroid_lon,
      regions[[reg_names[j]]]$centroid_lat, regions[[reg_names[j]]]$centroid_lon
    )
  }
}

# Detrended correlation matrix
detrend_cor <- cor(spd_diff, use = "pairwise.complete.obs")

# Extract upper triangles
dist_vec <- dist_matrix[upper.tri(dist_matrix)]
cor_vec <- detrend_cor[upper.tri(detrend_cor)]

# Pair labels for output
pair_labels <- c()
pair_dist <- c()
pair_cor <- c()
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    pair_labels <- c(pair_labels, paste(reg_labels[i], "x", reg_labels[j]))
    pair_dist <- c(pair_dist, dist_matrix[i, j])
    pair_cor <- c(pair_cor, detrend_cor[i, j])
  }
}

cat("Distance vs Detrended Correlation:\n")
cat(sprintf("%-40s  Dist(km)  Detrended r\n", "Pair"))
cat(paste(rep("-", 65), collapse = ""), "\n")
ord <- order(pair_dist)
for (k in ord) {
  cat(sprintf("%-40s  %7.0f   %+.3f\n", pair_labels[k], pair_dist[k], pair_cor[k]))
}

# Linear regression: correlation ~ distance
fit_dist <- lm(cor_vec ~ dist_vec)
cat(sprintf("\nLinear fit: r = %.4f + %.6f * distance\n", coef(fit_dist)[1], coef(fit_dist)[2]))
cat(sprintf("R-squared: %.3f\n", summary(fit_dist)$r.squared))
cat(sprintf("Slope p-value: %.4f\n", summary(fit_dist)$coefficients[2, 4]))

# Mantel test (permutation-based)
set.seed(42)
n_perm <- 9999
observed_mantel <- cor(dist_vec, -cor_vec)  # negative because we expect r to decrease with distance
cat(sprintf("\nObserved Mantel r (distance vs -correlation): %.4f\n", observed_mantel))

# Permutation
mantel_null <- numeric(n_perm)
n_items <- n_reg
for (p in 1:n_perm) {
  perm <- sample(1:n_items)
  perm_dist <- dist_matrix[perm, perm]
  perm_dist_vec <- perm_dist[upper.tri(perm_dist)]
  mantel_null[p] <- cor(perm_dist_vec, -cor_vec)
}
mantel_p <- (sum(mantel_null >= observed_mantel) + 1) / (n_perm + 1)
cat(sprintf("Mantel test p-value (one-sided): %.4f  (n_perm=%d)\n", mantel_p, n_perm))
cat(sprintf("Null distribution: mean=%.4f, sd=%.4f\n", mean(mantel_null), sd(mantel_null)))

if (mantel_p < 0.05) {
  cat("RESULT: Significant distance-decay. Synchrony decreases with distance.\n")
  cat("        Suggests diffusion / proximity effects contribute to synchrony.\n")
} else {
  cat("RESULT: No significant distance-decay. Synchrony is distance-independent.\n")
  cat("        Consistent with global climate forcing as driver.\n")
}

results$mantel_r <- observed_mantel
results$mantel_p <- mantel_p
results$dist_slope <- coef(fit_dist)[2]
results$dist_rsq <- summary(fit_dist)$r.squared

# =====================================================================
# ANALYSIS 2: Block Bootstrap Confidence Intervals
# =====================================================================
cat("\n================================================================\n")
cat("ANALYSIS 2: Block Bootstrap Confidence Intervals\n")
cat("================================================================\n\n")

set.seed(123)
n_boot <- 1000
block_size <- 5  # 5 time steps = 250 years
n_t <- nrow(spd_diff)

# Generate block bootstrap indices
boot_block_indices <- function(n, block_size) {
  n_blocks <- ceiling(n / block_size)
  start_points <- sample(1:(n - block_size + 1), n_blocks, replace = TRUE)
  idx <- unlist(lapply(start_points, function(s) s:(s + block_size - 1)))
  idx[1:n]
}

# Store bootstrap results
boot_results <- matrix(NA, nrow = n_boot, ncol = choose(n_reg, 2))
pair_idx <- 0
pair_names <- c()
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    pair_idx <- pair_idx + 1
    pair_names <- c(pair_names, paste(reg_labels[i], "x", reg_labels[j]))
  }
}
colnames(boot_results) <- pair_names

cat("Running", n_boot, "block bootstrap iterations (block size =", block_size, ")...\n")

for (b in 1:n_boot) {
  boot_idx <- boot_block_indices(n_t, block_size)
  boot_data <- spd_diff[boot_idx, ]
  pair_idx <- 0
  for (i in 1:(n_reg - 1)) {
    for (j in (i + 1):n_reg) {
      pair_idx <- pair_idx + 1
      boot_results[b, pair_idx] <- cor(boot_data[, i], boot_data[, j], use = "pairwise.complete.obs")
    }
  }
}

cat("Done.\n\n")

# Compute CIs and display
cat(sprintf("%-40s  Obs.r    95%% CI           Width  Sig?\n", "Pair"))
cat(paste(rep("-", 80), collapse = ""), "\n")

boot_ci_table <- data.frame(
  Pair = character(0), Obs_r = numeric(0),
  CI_lo = numeric(0), CI_hi = numeric(0),
  Significant = logical(0), stringsAsFactors = FALSE
)

pair_idx <- 0
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    pair_idx <- pair_idx + 1
    obs_r <- detrend_cor[i, j]
    ci <- quantile(boot_results[, pair_idx], probs = c(0.025, 0.975), na.rm = TRUE)
    sig <- ci[1] > 0  # CI doesn't cross zero
    width <- ci[2] - ci[1]
    cat(sprintf("%-40s  %+.3f   [%+.3f, %+.3f]   %.3f  %s\n",
                pair_names[pair_idx], obs_r, ci[1], ci[2], width,
                ifelse(sig, "YES", "NO")))
    boot_ci_table <- rbind(boot_ci_table, data.frame(
      Pair = pair_names[pair_idx], Obs_r = obs_r,
      CI_lo = ci[1], CI_hi = ci[2], Significant = sig,
      stringsAsFactors = FALSE
    ))
  }
}

n_boot_sig <- sum(boot_ci_table$Significant)
cat(sprintf("\nPairs with CI excluding zero: %d / %d\n", n_boot_sig, nrow(boot_ci_table)))
cat(sprintf("Mean CI width: %.3f\n", mean(boot_ci_table$CI_hi - boot_ci_table$CI_lo)))

results$boot_n_sig <- n_boot_sig
results$boot_ci_table <- boot_ci_table

# =====================================================================
# ANALYSIS 3: Effective Degrees of Freedom (Bartlett correction)
# =====================================================================
cat("\n================================================================\n")
cat("ANALYSIS 3: Effective Degrees of Freedom (Bartlett Correction)\n")
cat("================================================================\n\n")

# Compute lag-1 autocorrelation for each detrended series
lag1_ac <- numeric(n_reg)
names(lag1_ac) <- reg_names
for (i in 1:n_reg) {
  ac <- acf(spd_diff[, i], lag.max = 1, plot = FALSE, na.action = na.pass)
  lag1_ac[i] <- ac$acf[2]  # lag-1
}

cat("Lag-1 autocorrelations of detrended series:\n")
for (i in 1:n_reg) {
  cat(sprintf("  %-18s: rho1 = %.3f\n", reg_labels[i], lag1_ac[i]))
}
cat("\n")

# Bartlett formula: N_eff = N * (1 - rho1*rho2) / (1 + rho1*rho2)
N <- nrow(spd_diff)

cat(sprintf("%-40s  N_eff   t-stat   p(naive)   p(corrected)  Sig?\n", "Pair"))
cat(paste(rep("-", 95), collapse = ""), "\n")

bartlett_table <- data.frame(
  Pair = character(0), Obs_r = numeric(0),
  N_eff = numeric(0), p_naive = numeric(0),
  p_corrected = numeric(0), Significant = logical(0),
  stringsAsFactors = FALSE
)

for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    r <- detrend_cor[i, j]
    rho1 <- lag1_ac[i]
    rho2 <- lag1_ac[j]

    # Bartlett effective sample size
    N_eff <- N * (1 - rho1 * rho2) / (1 + rho1 * rho2)
    N_eff <- max(N_eff, 5)  # floor at 5

    # t-statistic with corrected df
    t_stat <- r * sqrt((N_eff - 2) / (1 - r^2))
    p_corrected <- 2 * (1 - pt(abs(t_stat), df = N_eff - 2))

    # Naive p-value (using full N)
    t_naive <- r * sqrt((N - 2) / (1 - r^2))
    p_naive <- 2 * (1 - pt(abs(t_naive), df = N - 2))

    sig <- p_corrected < 0.05
    pair_label <- paste(reg_labels[i], "x", reg_labels[j])

    cat(sprintf("%-40s  %5.1f   %6.2f   %.4e   %.4e   %s\n",
                pair_label, N_eff, t_stat, p_naive, p_corrected,
                ifelse(sig, "***", "ns")))

    bartlett_table <- rbind(bartlett_table, data.frame(
      Pair = pair_label, Obs_r = r, N_eff = N_eff,
      p_naive = p_naive, p_corrected = p_corrected,
      Significant = sig, stringsAsFactors = FALSE
    ))
  }
}

n_bartlett_sig <- sum(bartlett_table$Significant)
cat(sprintf("\nPairs significant after Bartlett correction: %d / %d\n",
            n_bartlett_sig, nrow(bartlett_table)))
cat(sprintf("Mean N_eff: %.1f (out of N=%d, ratio=%.1f%%)\n",
            mean(bartlett_table$N_eff), N, 100 * mean(bartlett_table$N_eff) / N))

results$bartlett_n_sig <- n_bartlett_sig
results$bartlett_table <- bartlett_table
results$mean_N_eff <- mean(bartlett_table$N_eff)

# =====================================================================
# ANALYSIS 4: Alternative Detrending — Loess Residuals (span=0.3)
# =====================================================================
cat("\n================================================================\n")
cat("ANALYSIS 4: Alternative Detrending — Loess Residuals (span=0.3)\n")
cat("================================================================\n\n")

spd_loess <- matrix(NA, nrow = nrow(spd_matrix), ncol = n_reg)
colnames(spd_loess) <- reg_names

for (i in 1:n_reg) {
  y <- spd_matrix[, i]
  x <- 1:length(y)
  lo <- loess(y ~ x, span = 0.3)
  spd_loess[, i] <- residuals(lo)
}

# Compute correlations on loess residuals
loess_cor <- cor(spd_loess, use = "pairwise.complete.obs")

cat("Comparison: First-Differenced vs Loess Residuals\n\n")
cat(sprintf("%-40s  r(diff)   r(loess)  Diff   Agreement?\n", "Pair"))
cat(paste(rep("-", 80), collapse = ""), "\n")

loess_table <- data.frame(
  Pair = character(0), r_diff = numeric(0),
  r_loess = numeric(0), stringsAsFactors = FALSE
)

for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    r_diff <- detrend_cor[i, j]
    r_lo <- loess_cor[i, j]
    agree <- (r_diff > 0 & r_lo > 0) | (r_diff < 0 & r_lo < 0)
    pair_label <- paste(reg_labels[i], "x", reg_labels[j])
    cat(sprintf("%-40s  %+.3f    %+.3f    %+.3f  %s\n",
                pair_label, r_diff, r_lo, r_lo - r_diff,
                ifelse(agree, "YES", "NO")))
    loess_table <- rbind(loess_table, data.frame(
      Pair = pair_label, r_diff = r_diff, r_loess = r_lo,
      stringsAsFactors = FALSE
    ))
  }
}

# Also compute p-values for loess
n_loess_sig <- 0
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    ct <- cor.test(spd_loess[, i], spd_loess[, j])
    if (ct$p.value < 0.05) n_loess_sig <- n_loess_sig + 1
  }
}

method_cor <- cor(loess_table$r_diff, loess_table$r_loess)
cat(sprintf("\nCorrelation between methods: r = %.3f\n", method_cor))
cat(sprintf("Mean r (first-diff): %.3f\n", mean(loess_table$r_diff)))
cat(sprintf("Mean r (loess):      %.3f\n", mean(loess_table$r_loess)))
cat(sprintf("Sign agreement: %d / %d pairs\n",
            sum((loess_table$r_diff > 0) == (loess_table$r_loess > 0)),
            nrow(loess_table)))
cat(sprintf("Pairs significant at p<0.05 (loess): %d / %d\n",
            n_loess_sig, nrow(loess_table)))

results$method_correlation <- method_cor
results$loess_mean_r <- mean(loess_table$r_loess)
results$loess_n_sig <- n_loess_sig

# =====================================================================
# ANALYSIS 5: Sub-regional Analysis within W.Europe
# =====================================================================
cat("\n================================================================\n")
cat("ANALYSIS 5: Sub-regional Analysis within W.Europe\n")
cat("================================================================\n\n")

# Define sub-regions of W.Europe
sub_regions <- list(
  France = list(
    filter = function(x) x[x$Country == "France", ],
    label = "France",
    centroid_lat = 46.5, centroid_lon = 2.5
  ),
  Iberia = list(
    filter = function(x) x[x$Country %in% c("Spain", "Portugal"), ],
    label = "Iberia",
    centroid_lat = 40.0, centroid_lon = -4.0
  ),
  Central = list(
    filter = function(x) x[x$Country %in% c("Germany", "Switzerland", "Austria"), ],
    label = "Central (DE/CH/AT)",
    centroid_lat = 48.5, centroid_lon = 11.0
  )
)

sub_data <- list()
for (rname in names(sub_regions)) {
  rd <- sub_regions[[rname]]$filter(d)
  n <- nrow(rd)
  cat(sprintf("  %-25s: %6d dates\n", sub_regions[[rname]]$label, n))
  if (n >= 100) sub_data[[rname]] <- rd
}
cat("\n")

# Compute SPDs for sub-regions
sub_spd_matrix <- matrix(NA, nrow = length(time_grid), ncol = length(sub_data))
colnames(sub_spd_matrix) <- names(sub_data)

for (i in seq_along(sub_data)) {
  rname <- names(sub_data)[i]
  rd <- sub_data[[rname]]
  cat(sprintf("  Calibrating %s (%d dates)...\n", sub_regions[[rname]]$label, nrow(rd)))
  caldates <- calibrate(x = rd$Age, errors = rd$Error, calCurves = "intcal20", verbose = FALSE)
  site_ids <- ifelse(rd$SiteName == "" | is.na(rd$SiteName),
                     paste0("unknown", seq_len(nrow(rd))),
                     gsub("_", "-", rd$SiteName))
  bins <- binPrep(sites = site_ids, ages = rd$Age, h = 200)
  spd_result <- spd(caldates, bins = bins, timeRange = time_range, verbose = FALSE)
  sg <- spd_result$grid
  spd_vals <- approx(sg$calBP, sg$PrDens, xout = time_grid, method = "linear")$y
  spd_vals <- spd_vals / mean(spd_vals, na.rm = TRUE)
  sub_spd_matrix[, i] <- spd_vals
}

sub_diff <- apply(sub_spd_matrix, 2, diff)
sub_labels <- sapply(names(sub_data), function(x) sub_regions[[x]]$label)
n_sub <- ncol(sub_spd_matrix)

cat("\n--- Intra-W.Europe correlations (detrended) ---\n\n")

intra_cors <- c()
for (i in 1:(n_sub - 1)) {
  for (j in (i + 1):n_sub) {
    r <- cor(sub_diff[, i], sub_diff[, j], use = "pairwise.complete.obs")
    ct <- cor.test(sub_diff[, i], sub_diff[, j])
    cat(sprintf("  %-30s x %-30s: r = %+.3f  p = %.4e\n",
                sub_labels[i], sub_labels[j], r, ct$p.value))
    intra_cors <- c(intra_cors, r)
  }
}
cat(sprintf("\nMean intra-W.Europe r: %.3f\n", mean(intra_cors)))

# Now compute cross-continental correlations (sub-regions vs non-European regions)
cat("\n--- Sub-regions vs Non-European regions (detrended) ---\n\n")

non_euro <- c("Near_East", "China", "Japan")
cross_cors <- c()

for (sub_name in names(sub_data)) {
  for (inter_name in non_euro) {
    if (inter_name %in% reg_names) {
      inter_idx <- which(reg_names == inter_name)
      r <- cor(sub_diff[, sub_name], spd_diff[, inter_idx], use = "pairwise.complete.obs")
      ct <- cor.test(sub_diff[, sub_name], spd_diff[, inter_idx])
      cat(sprintf("  %-20s x %-18s: r = %+.3f  p = %.4e\n",
                  sub_regions[[sub_name]]$label, reg_labels[inter_idx], r, ct$p.value))
      cross_cors <- c(cross_cors, r)
    }
  }
}
cat(sprintf("\nMean cross-continental r: %.3f\n", mean(cross_cors)))

# Test: is intra-European > inter-continental?
cat(sprintf("\nIntra-W.Europe mean r:    %.3f (n=%d pairs)\n", mean(intra_cors), length(intra_cors)))
cat(sprintf("Cross-continental mean r: %.3f (n=%d pairs)\n", mean(cross_cors), length(cross_cors)))
cat(sprintf("Difference: %.3f\n", mean(intra_cors) - mean(cross_cors)))

# Wilcoxon test (small sample, non-parametric)
if (length(intra_cors) >= 2 & length(cross_cors) >= 2) {
  wt <- wilcox.test(intra_cors, cross_cors, alternative = "greater")
  cat(sprintf("Wilcoxon test (intra > cross): W=%.0f, p=%.4f\n", wt$statistic, wt$p.value))
} else {
  cat("Insufficient pairs for Wilcoxon test; comparing numerically only.\n")
}

# Also compare with inter-continental pairs from main analysis
euro_regions <- c("Britain", "W_Europe", "E_Europe", "Scandinavia")
inter_continental_cors <- c()
intra_continental_cors <- c()

for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    r <- detrend_cor[i, j]
    ri <- reg_names[i]
    rj <- reg_names[j]
    both_euro <- (ri %in% euro_regions) & (rj %in% euro_regions)
    if (both_euro) {
      intra_continental_cors <- c(intra_continental_cors, r)
    } else {
      inter_continental_cors <- c(inter_continental_cors, r)
    }
  }
}

cat(sprintf("\n--- Full comparison (7-region main analysis) ---\n"))
cat(sprintf("Intra-European mean r:     %.3f (n=%d pairs)\n",
            mean(intra_continental_cors), length(intra_continental_cors)))
cat(sprintf("Inter-continental mean r:  %.3f (n=%d pairs)\n",
            mean(inter_continental_cors), length(inter_continental_cors)))
if (length(intra_continental_cors) >= 2 & length(inter_continental_cors) >= 2) {
  wt2 <- wilcox.test(intra_continental_cors, inter_continental_cors, alternative = "greater")
  cat(sprintf("Wilcoxon test (intra > inter): W=%.0f, p=%.4f\n", wt2$statistic, wt2$p.value))
}

results$intra_euro_mean <- mean(intra_cors)
results$cross_cont_mean <- mean(cross_cors)
results$intra_7reg_mean <- mean(intra_continental_cors)
results$inter_7reg_mean <- mean(inter_continental_cors)

# =====================================================================
# SUMMARY
# =====================================================================
cat("\n================================================================\n")
cat("COMPREHENSIVE SUMMARY\n")
cat("================================================================\n\n")

cat("1. DISTANCE-DECAY (Mantel test):\n")
cat(sprintf("   Mantel r = %.3f, p = %.4f\n", results$mantel_r, results$mantel_p))
cat(sprintf("   Linear fit R² = %.3f, slope = %.2e\n", results$dist_rsq, results$dist_slope))
cat("\n")

cat("2. BLOCK BOOTSTRAP CIs (1000 iterations, block=5):\n")
cat(sprintf("   Pairs with CI excluding zero: %d / 21\n", results$boot_n_sig))
cat(sprintf("   Mean CI width: %.3f\n", mean(boot_ci_table$CI_hi - boot_ci_table$CI_lo)))
cat("\n")

cat("3. BARTLETT CORRECTION:\n")
cat(sprintf("   Pairs significant after correction: %d / 21\n", results$bartlett_n_sig))
cat(sprintf("   Mean N_eff: %.1f (%.1f%% of N=%d)\n",
            results$mean_N_eff, 100 * results$mean_N_eff / N, N))
cat("\n")

cat("4. ALTERNATIVE DETRENDING (loess span=0.3):\n")
cat(sprintf("   Method correlation: r = %.3f\n", results$method_correlation))
cat(sprintf("   Mean r (first-diff): %.3f, Mean r (loess): %.3f\n",
            mean(loess_table$r_diff), results$loess_mean_r))
cat(sprintf("   Pairs significant (loess): %d / 21\n", results$loess_n_sig))
cat("\n")

cat("5. SUB-REGIONAL (W.Europe breakdown):\n")
cat(sprintf("   Intra-W.Europe mean r: %.3f\n", results$intra_euro_mean))
cat(sprintf("   Cross-continental mean r: %.3f\n", results$cross_cont_mean))
cat(sprintf("   Intra-European (7-reg): %.3f vs Inter-continental: %.3f\n",
            results$intra_7reg_mean, results$inter_7reg_mean))

cat("\n\nDone. All analyses complete.\n")
