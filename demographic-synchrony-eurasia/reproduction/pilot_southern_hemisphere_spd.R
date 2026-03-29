#!/usr/bin/env Rscript
# ==============================================================================
# Pilot: Southern Hemisphere SPD Synchrony as Calibration Artifact Test
# Purpose: Test whether SH regions (SHCal20) show synchrony with NH regions
#          (IntCal20). If yes, calibration curve cannot explain Paper 11 synchrony.
# ==============================================================================

library(rcarbon)

cat("=== Loading p3k14c data ===\n")
d <- read.csv("/home/ayu/archeco/shared/p3k14c_data.csv", stringsAsFactors = FALSE)

# Filter: valid coords, age 500-12000, error > 0 and < 500
d <- d[!is.na(d$Lat) & !is.na(d$Long) & !is.na(d$Age) & !is.na(d$Error), ]
d <- d[d$Age >= 500 & d$Age <= 12000 & d$Error > 0 & d$Error < 500, ]

cat("Total dates after filtering:", nrow(d), "\n\n")

# ==============================================================================
# Define regions
# ==============================================================================

# --- Northern Hemisphere (IntCal20) ---
nh_regions <- list(
  "Britain" = list(
    filter = function(x) x$Country %in% c("United Kingdom", "Ireland"),
    cal = "intcal20"
  ),
  "W.Europe" = list(
    filter = function(x) x$Country %in% c("France", "Spain", "Germany", "Belgium",
                                            "Netherlands", "Italy", "Portugal",
                                            "Switzerland", "Austria"),
    cal = "intcal20"
  ),
  "Scandinavia" = list(
    filter = function(x) x$Country %in% c("Norway", "Sweden", "Denmark", "Finland"),
    cal = "intcal20"
  )
)

# --- Southern Hemisphere (SHCal20) ---
sh_regions <- list(
  "S.America" = list(
    filter = function(x) x$Continent == "South America" & x$Lat < 0,
    cal = "shcal20"
  ),
  "Sub-Saharan.Africa" = list(
    filter = function(x) x$Continent == "Africa" & x$Lat < 0,
    cal = "shcal20"
  ),
  "Australia" = list(
    filter = function(x) x$Continent == "Australia" & x$Lat < 0,
    cal = "shcal20"
  )
)

all_regions <- c(nh_regions, sh_regions)

# ==============================================================================
# Count dates per region
# ==============================================================================
cat("=== Date counts per region ===\n")
region_data <- list()
for (name in names(all_regions)) {
  sub <- d[all_regions[[name]]$filter(d), ]
  region_data[[name]] <- sub
  cat(sprintf("  %-20s: %5d dates, %4d sites\n", name,
              nrow(sub), length(unique(sub$SiteName))))
}

# ==============================================================================
# Compute SPDs with site binning
# ==============================================================================
cat("\n=== Computing SPDs (site binning, 200yr bins) ===\n")

time_range <- c(10000, 1000)  # cal BP
spd_list <- list()

for (name in names(all_regions)) {
  sub <- region_data[[name]]
  if (nrow(sub) < 100) {
    cat(sprintf("  SKIP %s (only %d dates)\n", name, nrow(sub)))
    next
  }

  cal_curve <- all_regions[[name]]$cal
  cat(sprintf("  Calibrating %s (%d dates, %s)...\n", name, nrow(sub), cal_curve))

  # Calibrate
  cal <- calibrate(sub$Age, sub$Error, calCurves = cal_curve, verbose = FALSE)

  # Site binning (200yr bins) — use SiteName since SiteID is mostly empty
  # Replace underscores (binPrep doesn't allow them)
  site_labels <- gsub("_", "-", sub$SiteName)
  bins <- binPrep(sites = site_labels, ages = sub$Age, h = 200)

  # SPD
  spd_result <- spd(cal, bins = bins, timeRange = time_range, verbose = FALSE)

  # Extract and interpolate to 50-year resolution
  spd_df <- data.frame(
    calBP = spd_result$grid$calBP,
    PrDens = spd_result$grid$PrDens
  )

  # Interpolate to 50-year steps
  target_years <- seq(time_range[1], time_range[2], by = -50)
  interp <- approx(spd_df$calBP, spd_df$PrDens, xout = target_years)

  spd_vec <- interp$y
  # Normalize to mean = 1
  spd_vec <- spd_vec / mean(spd_vec, na.rm = TRUE)

  spd_list[[name]] <- data.frame(calBP = target_years, spd = spd_vec)
  cat(sprintf("  Done: %s (binned groups: %d)\n", name, length(unique(bins))))
}

cat(sprintf("\n  Successfully computed %d SPDs\n", length(spd_list)))

# ==============================================================================
# First-difference detrending
# ==============================================================================
cat("\n=== Detrending (first-differencing) ===\n")

detrended <- list()
for (name in names(spd_list)) {
  s <- spd_list[[name]]$spd
  detrended[[name]] <- diff(s)  # length N-1
}

n_steps <- length(detrended[[1]])
cat(sprintf("  Detrended series length: %d steps\n", n_steps))

# ==============================================================================
# Pairwise correlations
# ==============================================================================
cat("\n=== Pairwise correlations (detrended) ===\n\n")

region_names <- names(detrended)
n_reg <- length(region_names)

# Classification of pairs
results <- data.frame(
  pair = character(),
  r = numeric(),
  p = numeric(),
  type = character(),
  stringsAsFactors = FALSE
)

nh_names <- intersect(names(nh_regions), region_names)
sh_names <- intersect(names(sh_regions), region_names)

for (i in 1:(n_reg - 1)) {
  for (j in (i + 1):n_reg) {
    r1 <- region_names[i]
    r2 <- region_names[j]

    ct <- cor.test(detrended[[r1]], detrended[[r2]])

    # Classify pair
    if (r1 %in% nh_names & r2 %in% nh_names) {
      ptype <- "NH x NH (IntCal20 x IntCal20)"
    } else if (r1 %in% sh_names & r2 %in% sh_names) {
      ptype <- "SH x SH (SHCal20 x SHCal20)"
    } else {
      ptype <- "NH x SH (IntCal20 x SHCal20)"
    }

    results <- rbind(results, data.frame(
      pair = paste(r1, "x", r2),
      r = round(ct$estimate, 3),
      p = ct$p.value,
      type = ptype,
      stringsAsFactors = FALSE
    ))
  }
}

# Print results grouped by type
for (tp in unique(results$type)) {
  cat(sprintf("--- %s ---\n", tp))
  sub <- results[results$type == tp, ]
  for (k in 1:nrow(sub)) {
    sig <- ifelse(sub$p[k] < 0.05, "*", " ")
    sig <- ifelse(sub$p[k] < 0.01, "**", sig)
    sig <- ifelse(sub$p[k] < 0.001, "***", sig)
    cat(sprintf("  %-35s r = %+.3f  p = %.4f %s\n",
                sub$pair[k], sub$r[k], sub$p[k], sig))
  }
  cat(sprintf("  Mean r = %.3f\n\n", mean(sub$r)))
}

# ==============================================================================
# Summary statistics
# ==============================================================================
cat("=== SUMMARY ===\n")
nh_nh <- results[results$type == "NH x NH (IntCal20 x IntCal20)", ]
nh_sh <- results[results$type == "NH x SH (IntCal20 x SHCal20)", ]
sh_sh <- results[results$type == "SH x SH (SHCal20 x SHCal20)", ]

cat(sprintf("  NH x NH mean r: %.3f (n=%d pairs, Paper 11 full: 0.458)\n",
            mean(nh_nh$r), nrow(nh_nh)))
cat(sprintf("  NH x SH mean r: %.3f (n=%d pairs) <- KEY TEST\n",
            mean(nh_sh$r), nrow(nh_sh)))
cat(sprintf("  SH x SH mean r: %.3f (n=%d pairs)\n",
            mean(sh_sh$r), nrow(sh_sh)))

cat(sprintf("\n  NH x SH positive pairs: %d/%d (%.0f%%)\n",
            sum(nh_sh$r > 0), nrow(nh_sh),
            100 * sum(nh_sh$r > 0) / nrow(nh_sh)))
cat(sprintf("  NH x SH significant (p<0.05): %d/%d (%.0f%%)\n",
            sum(nh_sh$p < 0.05), nrow(nh_sh),
            100 * sum(nh_sh$p < 0.05) / nrow(nh_sh)))

# Welch t-test comparing NH×NH vs NH×SH correlations
if (nrow(nh_nh) >= 2 & nrow(nh_sh) >= 2) {
  tt <- t.test(nh_nh$r, nh_sh$r)
  cat(sprintf("\n  Welch t-test NH×NH vs NH×SH: t=%.2f, p=%.4f\n",
              tt$statistic, tt$p.value))
  cat(sprintf("  -> %s\n",
              ifelse(tt$p.value > 0.05,
                     "NOT significantly different = calibration artifact REJECTED",
                     "Significantly different = more investigation needed")))
}

# ==============================================================================
# Interpretation
# ==============================================================================
cat("\n=== INTERPRETATION ===\n")
if (mean(nh_sh$r) > 0.1 & sum(nh_sh$p < 0.05) > nrow(nh_sh) / 2) {
  cat("  STRONG EVIDENCE against calibration artifact hypothesis.\n")
  cat("  NH-SH synchrony exists despite DIFFERENT calibration curves.\n")
  cat("  This means Paper 11's synchrony cannot be explained by IntCal20 alone.\n")
} else if (mean(nh_sh$r) > 0) {
  cat("  MODERATE EVIDENCE: NH-SH correlations are positive but weaker.\n")
  cat("  Some real synchrony exists, but calibration may amplify it.\n")
} else {
  cat("  WEAK/NO EVIDENCE: NH-SH correlations near zero.\n")
  cat("  Cannot rule out calibration artifact for Paper 11 findings.\n")
}

# ==============================================================================
# Save full results table
# ==============================================================================
write.csv(results, "/home/ayu/archeco/shared/pilot_sh_correlation_results.csv",
          row.names = FALSE)
cat("\n  Results saved to shared/pilot_sh_correlation_results.csv\n")

cat("\n=== PILOT COMPLETE ===\n")
