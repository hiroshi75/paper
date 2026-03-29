#!/usr/bin/env Rscript
# =============================================================================
# Pilot: Demographic Synchrony Networks from p3k14c Radiocarbon Dates
# =============================================================================
# Question: Do population collapses (SPD drops) occur simultaneously across
# distant regions? Synchrony → shared driver (climate). Asynchrony → local causes.
#
# Method:
#   1. Load p3k14c, define sub-regions (W.Europe, E.Europe, Near East,
#      E.Asia-China, E.Asia-Japan, Britain)
#   2. Calibrate dates and compute regional SPDs with site binning
#   3. Compute raw cross-correlations between regional SPD time series
#   4. Apply first-differencing (detrending) to remove shared long-term trends
#   5. Re-compute cross-correlations on detrended series
#   6. Identify synchronous collapse episodes, especially around known events
#      (8.2ka, 4.2ka, 3.2ka/Bronze Age Collapse, 2.2ka)
# =============================================================================

library(rcarbon)

cat("========================================\n")
cat("PILOT: Archaeological Demographic Synchrony\n")
cat("========================================\n\n")

# --- 1. Load and filter data ---
d <- read.csv("/home/ayu/archeco/shared/p3k14c_data.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$Lat) & !is.na(d$Long) & !is.na(d$Age) & !is.na(d$Error), ]
d <- d[d$Error > 0 & d$Error < 500 & d$Age > 500 & d$Age < 12000, ]

cat("After filtering: ", nrow(d), " dates\n\n")

# --- 2. Define sub-regions ---
# Using bounding boxes for clear geographic separation
regions <- list(
  Britain = list(
    filter = function(x) x[x$Country %in% c("United Kingdom", "Ireland"), ],
    label = "Britain+Ireland"
  ),
  W_Europe = list(
    filter = function(x) x[x$Country %in% c("France", "Spain", "Germany",
                                              "Belgium", "Netherlands", "Italy",
                                              "Portugal", "Switzerland", "Austria") &
                            x$Continent == "Europe", ],
    label = "W.Europe"
  ),
  E_Europe = list(
    filter = function(x) x[x$Country %in% c("Poland", "Hungary", "Czech Republic",
                                              "Romania", "Bulgaria", "Serbia",
                                              "Croatia", "Slovakia", "Greece") &
                            x$Continent == "Europe", ],
    label = "E.Europe+Greece"
  ),
  Scandinavia = list(
    filter = function(x) x[x$Country %in% c("Norway", "Sweden", "Denmark", "Finland"), ],
    label = "Scandinavia"
  ),
  Near_East = list(
    filter = function(x) x[x$Country %in% c("Turkey", "Syria", "Iraq", "Iran",
                                              "Israel", "Jordan", "Lebanon", "Palestine"), ],
    label = "Near East"
  ),
  China = list(
    filter = function(x) x[x$Country == "China", ],
    label = "China"
  ),
  Japan = list(
    filter = function(x) x[x$Country == "Japan", ],
    label = "Japan"
  )
)

# Extract and count
region_data <- list()
for (rname in names(regions)) {
  rd <- regions[[rname]]$filter(d)
  n <- if (is.null(rd)) 0 else nrow(rd)
  cat(sprintf("%-18s: %6d dates\n", regions[[rname]]$label, n))
  if (n >= 200) {
    region_data[[rname]] <- rd
  }
}
cat("\n")

# --- 3. Calibrate and compute SPDs per region ---
cat("Calibrating dates and computing SPDs...\n")
cat("(This may take several minutes)\n\n")

# Common time range: 10000-1000 cal BP (8000 BC - AD 950)
time_range <- c(10000, 1000)

spd_list <- list()
for (rname in names(region_data)) {
  rd <- region_data[[rname]]
  cat(sprintf("  Calibrating %s (%d dates)...\n", regions[[rname]]$label, nrow(rd)))

  # Calibrate
  caldates <- calibrate(
    x = rd$Age,
    errors = rd$Error,
    calCurves = "intcal20",
    verbose = FALSE
  )

  # Site binning: use SiteName as site identifier, 200-year bins
  # This reduces overrepresentation of heavily-dated sites
  if (any(rd$SiteName != "" & !is.na(rd$SiteName))) {
    site_ids <- ifelse(rd$SiteName == "" | is.na(rd$SiteName),
                       paste0("unknown", seq_len(nrow(rd))),
                       gsub("_", "-", rd$SiteName))
    bins <- binPrep(sites = site_ids, ages = rd$Age, h = 200)
  } else {
    bins <- NA
  }

  # Compute SPD
  spd_result <- spd(
    caldates,
    bins = bins,
    timeRange = time_range,
    verbose = FALSE
  )

  spd_list[[rname]] <- spd_result
  cat(sprintf("    Done. Binned into %d bins.\n",
              length(unique(bins[!is.na(bins)]))))
}

cat("\nAll SPDs computed.\n\n")

# --- 4. Extract SPD values into a common time grid ---
# Use 50-year resolution for analysis
grid_step <- 50
time_grid <- seq(time_range[1], time_range[2], by = -grid_step)

spd_matrix <- matrix(NA, nrow = length(time_grid), ncol = length(spd_list))
colnames(spd_matrix) <- names(spd_list)

for (i in seq_along(spd_list)) {
  rname <- names(spd_list)[i]
  sg <- spd_list[[rname]]$grid
  # Interpolate to common grid
  spd_vals <- approx(sg$calBP, sg$PrDens, xout = time_grid, method = "linear")$y
  # Normalize each region to mean=1 for comparability
  spd_vals <- spd_vals / mean(spd_vals, na.rm = TRUE)
  spd_matrix[, i] <- spd_vals
}

cat("SPD matrix: ", nrow(spd_matrix), " time steps x ", ncol(spd_matrix), " regions\n\n")

# --- 5. Raw cross-correlations ---
cat("=== RAW CROSS-CORRELATIONS (lag 0) ===\n")
n_reg <- ncol(spd_matrix)
raw_cor <- cor(spd_matrix, use = "pairwise.complete.obs")
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    cat(sprintf("  %-18s x %-18s: r = %+.3f\n",
                regions[[colnames(spd_matrix)[i]]]$label,
                regions[[colnames(spd_matrix)[j]]]$label,
                raw_cor[i, j]))
  }
}

# --- 6. DETREND: First-differencing ---
cat("\n=== DETRENDED (first-difference) CROSS-CORRELATIONS ===\n")
spd_diff <- apply(spd_matrix, 2, diff)

detrend_cor <- cor(spd_diff, use = "pairwise.complete.obs")
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    # Significance test
    ct <- cor.test(spd_diff[, i], spd_diff[, j])
    sig <- ifelse(ct$p.value < 0.001, "***",
           ifelse(ct$p.value < 0.01, "**",
           ifelse(ct$p.value < 0.05, "*", "ns")))
    cat(sprintf("  %-18s x %-18s: r = %+.3f  p = %.4f  %s\n",
                regions[[colnames(spd_matrix)[i]]]$label,
                regions[[colnames(spd_matrix)[j]]]$label,
                detrend_cor[i, j], ct$p.value, sig))
  }
}

# --- 7. Focus on climate events ---
cat("\n=== REGIONAL SPD VALUES AT KNOWN CLIMATE EVENTS ===\n")
events <- data.frame(
  name = c("8.2ka event", "5.9ka aridification", "4.2ka event",
           "3.2ka Late Bronze Age", "2.2ka Roman Warm/end"),
  calBP = c(8200, 5900, 4200, 3200, 2200),
  stringsAsFactors = FALSE
)

# For each event, look at SPD in a window
window <- 200  # +/- 200 years

cat("\nZ-scores at event windows (negative = below-average population):\n\n")
cat(sprintf("%-25s", "Event"))
for (rname in colnames(spd_matrix)) {
  cat(sprintf("%-14s", regions[[rname]]$label))
}
cat("\n")
cat(paste(rep("-", 25 + 14 * n_reg), collapse = ""), "\n")

for (e in 1:nrow(events)) {
  # Find time indices in window
  event_idx <- which(time_grid >= (events$calBP[e] - window) &
                     time_grid <= (events$calBP[e] + window))

  cat(sprintf("%-25s", paste0(events$name[e], " (", events$calBP[e], ")")))

  for (rname in colnames(spd_matrix)) {
    # Z-score: how does event window compare to overall mean
    event_mean <- mean(spd_matrix[event_idx, rname], na.rm = TRUE)
    overall_mean <- mean(spd_matrix[, rname], na.rm = TRUE)
    overall_sd <- sd(spd_matrix[, rname], na.rm = TRUE)
    z <- (event_mean - overall_mean) / overall_sd

    marker <- ifelse(z < -1, " <<", ifelse(z > 1, " >>", ""))
    cat(sprintf("%-14s", paste0(sprintf("%+.2f", z), marker)))
  }
  cat("\n")
}

cat("\n<< = notably below average (potential collapse)\n")
cat(">> = notably above average (expansion)\n")

# --- 8. Rolling window synchrony ---
cat("\n=== ROLLING SYNCHRONY (500-year windows, first-differenced) ===\n")
window_size <- 10  # 10 steps of 50 years = 500 years
n_windows <- nrow(spd_diff) - window_size + 1

# Compute mean pairwise correlation in each window
rolling_sync <- data.frame(
  calBP = numeric(n_windows),
  mean_r = numeric(n_windows),
  n_sig_pairs = numeric(n_windows)
)

for (w in 1:n_windows) {
  idx <- w:(w + window_size - 1)
  window_data <- spd_diff[idx, ]

  # Pairwise correlations
  pairs_r <- c()
  for (i in 1:(n_reg - 1)) {
    for (j in (i + 1):n_reg) {
      r <- cor(window_data[, i], window_data[, j], use = "pairwise.complete.obs")
      if (!is.na(r)) pairs_r <- c(pairs_r, r)
    }
  }

  rolling_sync$calBP[w] <- time_grid[w + window_size %/% 2]
  rolling_sync$mean_r[w] <- mean(pairs_r, na.rm = TRUE)
  rolling_sync$n_sig_pairs[w] <- sum(pairs_r > 0.5, na.rm = TRUE)
}

# Find peaks of synchrony
high_sync <- rolling_sync[rolling_sync$mean_r > quantile(rolling_sync$mean_r, 0.9, na.rm=TRUE), ]
cat("\nTop 10% synchrony episodes (mean pairwise r):\n")
for (i in 1:min(15, nrow(high_sync))) {
  cat(sprintf("  %d cal BP: mean r = %.3f (%d strongly correlated pairs)\n",
              high_sync$calBP[i], high_sync$mean_r[i], high_sync$n_sig_pairs[i]))
}

low_sync <- rolling_sync[rolling_sync$mean_r < quantile(rolling_sync$mean_r, 0.1, na.rm=TRUE), ]
cat("\nBottom 10% (asynchrony episodes):\n")
for (i in 1:min(10, nrow(low_sync))) {
  cat(sprintf("  %d cal BP: mean r = %.3f\n",
              low_sync$calBP[i], low_sync$mean_r[i]))
}

# --- 9. Generate plot ---
cat("\nGenerating plots...\n")
pdf("/home/ayu/archeco/shared/pilot_archaeological_network_plots.pdf",
    width = 12, height = 16)

par(mfrow = c(4, 1), mar = c(4, 4, 3, 1))

# Panel 1: Raw regional SPDs
colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
            "#ff7f00", "#a65628", "#f781bf")
plot(time_grid, spd_matrix[, 1], type = "l", col = colors[1],
     xlim = rev(range(time_grid)), ylim = c(0, max(spd_matrix, na.rm = TRUE) * 1.1),
     xlab = "cal BP", ylab = "Normalized SPD",
     main = "Regional SPDs (normalized, site-binned)")
for (i in 2:ncol(spd_matrix)) {
  lines(time_grid, spd_matrix[, i], col = colors[i])
}
legend("topleft",
       legend = sapply(colnames(spd_matrix), function(x) regions[[x]]$label),
       col = colors[1:ncol(spd_matrix)], lty = 1, cex = 0.8, ncol = 2)
# Mark events
abline(v = events$calBP, lty = 2, col = "gray50")
text(events$calBP, max(spd_matrix, na.rm = TRUE) * 1.05, events$name,
     srt = 90, adj = 0, cex = 0.7)

# Panel 2: Detrended (first-differenced) SPDs
plot(time_grid[-1], spd_diff[, 1], type = "l", col = colors[1],
     xlim = rev(range(time_grid[-1])),
     ylim = range(spd_diff, na.rm = TRUE),
     xlab = "cal BP", ylab = "First difference",
     main = "Detrended (first-differenced) regional SPDs")
for (i in 2:ncol(spd_diff)) {
  lines(time_grid[-1], spd_diff[, i], col = colors[i])
}
abline(h = 0, lty = 2)
abline(v = events$calBP, lty = 2, col = "gray50")

# Panel 3: Rolling synchrony
plot(rolling_sync$calBP, rolling_sync$mean_r, type = "l", col = "black", lwd = 2,
     xlim = rev(range(rolling_sync$calBP)),
     xlab = "cal BP", ylab = "Mean pairwise r",
     main = "Rolling inter-regional synchrony (500-yr windows, detrended)")
abline(h = 0, lty = 2)
abline(h = mean(rolling_sync$mean_r, na.rm = TRUE), lty = 3, col = "blue")
abline(v = events$calBP, lty = 2, col = "gray50")
text(events$calBP, max(rolling_sync$mean_r, na.rm = TRUE) * 0.95, events$name,
     srt = 90, adj = 0, cex = 0.7)

# Panel 4: Correlation heatmap (raw vs detrended)
par(mfrow = c(1, 2), mar = c(5, 6, 3, 2))
labels <- sapply(colnames(spd_matrix), function(x) regions[[x]]$label)

# Raw correlations
image(1:n_reg, 1:n_reg, raw_cor, col = hcl.colors(20, "RdBu", rev = TRUE),
      zlim = c(-1, 1), axes = FALSE,
      xlab = "", ylab = "", main = "Raw correlations")
axis(1, at = 1:n_reg, labels = labels, las = 2, cex.axis = 0.8)
axis(2, at = 1:n_reg, labels = labels, las = 1, cex.axis = 0.8)
for (i in 1:n_reg) for (j in 1:n_reg) {
  text(i, j, sprintf("%.2f", raw_cor[i, j]), cex = 0.7)
}

# Detrended correlations
image(1:n_reg, 1:n_reg, detrend_cor, col = hcl.colors(20, "RdBu", rev = TRUE),
      zlim = c(-1, 1), axes = FALSE,
      xlab = "", ylab = "", main = "Detrended correlations")
axis(1, at = 1:n_reg, labels = labels, las = 2, cex.axis = 0.8)
axis(2, at = 1:n_reg, labels = labels, las = 1, cex.axis = 0.8)
for (i in 1:n_reg) for (j in 1:n_reg) {
  text(i, j, sprintf("%.2f", detrend_cor[i, j]), cex = 0.7)
}

dev.off()
cat("Plots saved to: shared/pilot_archaeological_network_plots.pdf\n")

# --- 10. Summary statistics for report ---
cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

cat("Raw correlation range: ", sprintf("%.3f to %.3f\n",
    min(raw_cor[upper.tri(raw_cor)]), max(raw_cor[upper.tri(raw_cor)])))
cat("Detrended correlation range: ", sprintf("%.3f to %.3f\n",
    min(detrend_cor[upper.tri(detrend_cor)]), max(detrend_cor[upper.tri(detrend_cor)])))

cat("\nMean raw pairwise r: ", sprintf("%.3f\n", mean(raw_cor[upper.tri(raw_cor)])))
cat("Mean detrended pairwise r: ", sprintf("%.3f\n", mean(detrend_cor[upper.tri(detrend_cor)])))

# How many pairs survive detrending with r > 0.2?
n_pairs <- sum(upper.tri(detrend_cor))
n_positive <- sum(detrend_cor[upper.tri(detrend_cor)] > 0.2)
cat(sprintf("\nPairs with detrended r > 0.2: %d / %d (%.0f%%)\n",
            n_positive, n_pairs, 100 * n_positive / n_pairs))

sig_pairs <- 0
for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    ct <- cor.test(spd_diff[, i], spd_diff[, j])
    if (ct$p.value < 0.05) sig_pairs <- sig_pairs + 1
  }
}
cat(sprintf("Pairs significant at p<0.05: %d / %d\n", sig_pairs, n_pairs))

cat("\nDone.\n")
