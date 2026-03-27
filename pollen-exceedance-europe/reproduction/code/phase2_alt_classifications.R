#!/usr/bin/env Rscript
# ============================================================================
# Phase 2: Alternative Classifications to Break Circularity
# ============================================================================
# Problem: Original H1 = pastoral/arable taxa exceed 2SD threshold.
#   Finding H1 has higher pastoral BC% is tautological.
# Solution: Test 3 independent classification methods.
#
# Alt1: Indicator taxa PRESENCE (binary: Plantago lanceolata OR Cerealia >0%)
# Alt2: Temporal clustering (onset within ±1500yr of regional Neolithic arrival)
# Alt3: Geographic position relative to Neolithic frontier (CE only, needs lat/lon)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(irr)  # for kappa
})

# Try to load lme4 for LMM; fall back to simple pooled test if unavailable
has_lme4 <- requireNamespace("lme4", quietly = TRUE)
if (has_lme4) library(lme4)

cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("Phase 2: Alternative Classifications for Circularity Breaking\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# ============================================================================
# 1. LOAD AND HARMONIZE DATA
# ============================================================================

# --- UK decomposition ---
uk_dec <- read.csv("shared/signal_decomposition_uk_corrected.csv",
                   stringsAsFactors = FALSE) %>%
  mutate(region = "UK",
         pastoral_arable_pct = pastoral_arable_pct,
         total_bc = total_bc,
         original_H1 = classification == "H1_anthropogenic")

# --- Scandinavia decomposition ---
scand_dec <- read.csv("shared/signal_decomposition_scandinavia.csv",
                      stringsAsFactors = FALSE) %>%
  mutate(region = "Scandinavia",
         pastoral_arable_pct = bc_anthro_pct,
         total_bc = bc_total,
         original_H1 = classification == "H1_anthropogenic")

# --- Central Europe decomposition ---
ce_dec <- read.csv("shared/signal_decomposition_central_europe.csv",
                   stringsAsFactors = FALSE) %>%
  mutate(region = "Central_Europe",
         original_H1 = classification == "H1_anthropogenic")

# Combine all with harmonized columns
all_dec <- bind_rows(
  uk_dec %>% select(datasetid, sitename, region, signal_onset_age,
                    pastoral_arable_pct, total_bc, original_H1,
                    classification, any_of(c("lat", "lon"))),
  scand_dec %>% select(datasetid, sitename, region, signal_onset_age,
                       pastoral_arable_pct, total_bc, original_H1,
                       classification),
  ce_dec %>% select(datasetid, sitename, region, signal_onset_age,
                    pastoral_arable_pct, total_bc, original_H1,
                    classification, any_of(c("lat", "lon")))
)

cat("Combined dataset: n =", nrow(all_dec), "sites\n")
cat("  UK:", sum(all_dec$region == "UK"),
    "  Scandinavia:", sum(all_dec$region == "Scandinavia"),
    "  CE:", sum(all_dec$region == "Central_Europe"), "\n")
cat("  Original H1:", sum(all_dec$original_H1),
    "  H3:", sum(!all_dec$original_H1), "\n\n")

# --- Indicator data ---
uk_ind <- read.csv("shared/signal_phase_indicator_test.csv",
                   stringsAsFactors = FALSE)
scand_ce_ind <- read.csv("shared/indicator_test_scand_ce.csv",
                         stringsAsFactors = FALSE)

# ============================================================================
# 2. ALTERNATIVE 1: Indicator Taxa PRESENCE (binary)
# ============================================================================
# H1_alt1 = site where Plantago lanceolata OR Cerealia-type is PRESENT (>0%)
# at any post-7000BP sample. "Present" = first_appearance_age is not NA
# and the taxon appears somewhere in the record.

cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("ALTERNATIVE 1: Indicator Taxa Presence (binary)\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# UK: Check Plantago_lanceolata or Cerealia presence
uk_presence <- uk_ind %>%
  filter(indicator %in% c("Plantago_lanceolata", "Cerealia")) %>%
  group_by(datasetid, sitename) %>%
  summarise(
    has_plantago = any(indicator == "Plantago_lanceolata" & !is.na(first_appearance_age)),
    has_cerealia = any(indicator == "Cerealia" & !is.na(first_appearance_age)),
    .groups = "drop"
  ) %>%
  mutate(alt1_H1 = has_plantago | has_cerealia)

# Scandinavia/CE: Check Plantago or Cerealia presence
scand_ce_presence <- scand_ce_ind %>%
  filter(indicator %in% c("Plantago", "Cerealia")) %>%
  group_by(datasetid, sitename, region) %>%
  summarise(
    has_plantago = any(indicator == "Plantago" & !is.na(first_appearance_age)),
    has_cerealia = any(indicator == "Cerealia" & !is.na(first_appearance_age)),
    .groups = "drop"
  ) %>%
  mutate(alt1_H1 = has_plantago | has_cerealia)

# Merge with decomposition data — deduplicate by datasetid first
alt1_lookup <- bind_rows(
  uk_presence %>% select(datasetid, alt1_H1),
  scand_ce_presence %>% select(datasetid, alt1_H1)
) %>%
  distinct(datasetid, .keep_all = TRUE)

all_dec <- all_dec %>%
  left_join(alt1_lookup, by = "datasetid")

# Sites without indicator data: use NA
cat("Alt1 classification:\n")
cat("  H1_alt1 (indicator present):", sum(all_dec$alt1_H1, na.rm=TRUE), "\n")
cat("  H0_alt1 (no indicator):", sum(!all_dec$alt1_H1, na.rm=TRUE), "\n")
cat("  NA (no indicator data):", sum(is.na(all_dec$alt1_H1)), "\n\n")

# ============================================================================
# 3. ALTERNATIVE 2: Temporal Clustering
# ============================================================================
# H1_alt2 = signal onset within ±1500yr of regional Neolithic arrival
# Regional Neolithic arrival estimates (cal BP):
#   UK/Ireland: ~6000 BP (Early Neolithic)
#   Scandinavia: ~5900 BP (southern) to ~4500 BP (northern) — use ~5200 mean
#   Central Europe: ~7500 BP (LBK) to ~6000 BP — use ~6500 mean

cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("ALTERNATIVE 2: Temporal Clustering (±1500yr of Neolithic)\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

neolithic_arrival <- c(UK = 6000, Scandinavia = 5200, Central_Europe = 6500)

all_dec <- all_dec %>%
  mutate(
    expected_neolithic = neolithic_arrival[region],
    onset_neo_diff = abs(signal_onset_age - expected_neolithic),
    alt2_H1 = onset_neo_diff <= 1500
  )

cat("Regional Neolithic arrivals used:\n")
for (r in names(neolithic_arrival)) {
  n_h1 <- sum(all_dec$region == r & all_dec$alt2_H1, na.rm=TRUE)
  n_total <- sum(all_dec$region == r)
  cat(sprintf("  %s: %d BP (H1_alt2: %d/%d)\n", r, neolithic_arrival[r], n_h1, n_total))
}
cat("\nAlt2 classification:\n")
cat("  H1_alt2 (near Neolithic):", sum(all_dec$alt2_H1, na.rm=TRUE), "\n")
cat("  H0_alt2 (far from Neolithic):", sum(!all_dec$alt2_H1, na.rm=TRUE), "\n\n")

# ============================================================================
# 4. ALTERNATIVE 3: Geographic Position (CE + UK with coords only)
# ============================================================================
# Expected Neolithic age = f(latitude, longitude) based on SE→NW gradient
# Simple linear model: expected_age = a + b*lat + c*lon
# Using well-known dates: LBK in Hungary ~7500 BP, Britain ~6000 BP, Scandinavia ~5000 BP

cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("ALTERNATIVE 3: Geographic Position (Neolithic frontier model)\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# Build simple Neolithic frontier model from anchor points
# These are approximate Neolithic arrival dates at key locations
anchor_points <- data.frame(
  lat = c(47.0, 52.0, 44.0, 48.5, 51.5, 55.0, 57.0, 60.0, 43.5, 54.0),
  lon = c(20.0, 12.0,  7.0,  2.0, -1.0, -6.0, -4.0, 10.0, -1.0, -6.5),
  neo_age = c(7500, 7000, 7000, 6500, 6000, 5800, 5500, 5000, 7200, 5800)
)

neo_model <- lm(neo_age ~ lat + lon, data = anchor_points)
cat("Neolithic frontier linear model:\n")
cat(sprintf("  expected_age = %.0f + %.1f*lat + %.1f*lon\n",
            coef(neo_model)[1], coef(neo_model)[2], coef(neo_model)[3]))
cat(sprintf("  R² = %.3f\n\n", summary(neo_model)$r.squared))

# Apply only to sites with lat/lon
sites_with_coords <- all_dec %>% filter(!is.na(lat) & !is.na(lon))
cat("Sites with coordinates:", nrow(sites_with_coords), "out of", nrow(all_dec), "\n")

if (nrow(sites_with_coords) > 0) {
  sites_with_coords <- sites_with_coords %>%
    mutate(
      expected_neo_geo = predict(neo_model, newdata = data.frame(lat = lat, lon = lon)),
      onset_neo_geo_diff = abs(signal_onset_age - expected_neo_geo),
      alt3_H1 = onset_neo_geo_diff <= 1000
    )

  cat("Alt3 classification (±1000yr of expected local Neolithic):\n")
  cat("  H1_alt3:", sum(sites_with_coords$alt3_H1, na.rm=TRUE), "\n")
  cat("  H0_alt3:", sum(!sites_with_coords$alt3_H1, na.rm=TRUE), "\n\n")

  # Merge back — deduplicate to avoid many-to-many
  alt3_lookup <- sites_with_coords %>%
    select(datasetid, region, alt3_H1, expected_neo_geo, onset_neo_geo_diff) %>%
    distinct(datasetid, region, .keep_all = TRUE)
  all_dec <- all_dec %>%
    left_join(alt3_lookup, by = c("datasetid", "region"))
} else {
  all_dec$alt3_H1 <- NA
  cat("  No sites with coordinates available for Alt3.\n\n")
}

# ============================================================================
# 5. STATISTICAL TESTS
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("STATISTICAL TESTS\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

results_list <- list()

run_comparison <- function(data, class_var, class_label) {
  d <- data %>% filter(!is.na(!!sym(class_var)))

  if (sum(d[[class_var]]) < 3 || sum(!d[[class_var]]) < 3) {
    cat(sprintf("\n--- %s: SKIPPED (insufficient group sizes: H1=%d, H0=%d) ---\n",
                class_label, sum(d[[class_var]]), sum(!d[[class_var]])))
    return(NULL)
  }

  h1_pastoral <- d$pastoral_arable_pct[d[[class_var]]]
  h0_pastoral <- d$pastoral_arable_pct[!d[[class_var]]]
  h1_bc <- d$total_bc[d[[class_var]]]
  h0_bc <- d$total_bc[!d[[class_var]]]

  cat(sprintf("\n--- %s ---\n", class_label))
  cat(sprintf("  H1: n=%d, H0: n=%d\n", length(h1_pastoral), length(h0_pastoral)))

  # Pastoral arable pct comparison
  cat("\n  [Pastoral/Arable %]\n")
  cat(sprintf("    H1 median: %.2f%%  (IQR: %.2f–%.2f)\n",
              median(h1_pastoral), quantile(h1_pastoral, 0.25), quantile(h1_pastoral, 0.75)))
  cat(sprintf("    H0 median: %.2f%%  (IQR: %.2f–%.2f)\n",
              median(h0_pastoral), quantile(h0_pastoral, 0.25), quantile(h0_pastoral, 0.75)))

  mw_pastoral <- wilcox.test(h1_pastoral, h0_pastoral)
  r_pastoral <- abs(qnorm(mw_pastoral$p.value / 2)) / sqrt(nrow(d))
  cat(sprintf("    Mann-Whitney U: W=%.0f, p=%.4f, r=%.3f\n",
              mw_pastoral$statistic, mw_pastoral$p.value, r_pastoral))

  # Total BC comparison
  cat("\n  [Total BC]\n")
  cat(sprintf("    H1 median: %.3f  (IQR: %.3f–%.3f)\n",
              median(h1_bc), quantile(h1_bc, 0.25), quantile(h1_bc, 0.75)))
  cat(sprintf("    H0 median: %.3f  (IQR: %.3f–%.3f)\n",
              median(h0_bc), quantile(h0_bc, 0.25), quantile(h0_bc, 0.75)))

  mw_bc <- wilcox.test(h1_bc, h0_bc)
  r_bc <- abs(qnorm(mw_bc$p.value / 2)) / sqrt(nrow(d))
  cat(sprintf("    Mann-Whitney U: W=%.0f, p=%.4f, r=%.3f\n",
              mw_bc$statistic, mw_bc$p.value, r_bc))

  # LMM with region as random effect (if lme4 available and >1 region)
  lmm_p_pastoral <- NA
  lmm_p_bc <- NA
  if (has_lme4 && length(unique(d$region)) > 1) {
    tryCatch({
      d$class_binary <- as.numeric(d[[class_var]])

      m_pastoral <- lmer(pastoral_arable_pct ~ class_binary + (1|region), data = d)
      coef_mat <- summary(m_pastoral)$coefficients
      t_val <- coef_mat[2, "t value"]
      lmm_p_pastoral <- 2 * pt(abs(t_val), df = nrow(d) - 3, lower.tail = FALSE)
      cat(sprintf("\n  [LMM: pastoral_pct ~ class + (1|region)]\n"))
      cat(sprintf("    Fixed effect (class): est=%.3f, t=%.3f, p=%.4f\n",
                  fixef(m_pastoral)[2], t_val, lmm_p_pastoral))

      m_bc <- lmer(total_bc ~ class_binary + (1|region), data = d)
      coef_mat_bc <- summary(m_bc)$coefficients
      t_val_bc <- coef_mat_bc[2, "t value"]
      lmm_p_bc <- 2 * pt(abs(t_val_bc), df = nrow(d) - 3, lower.tail = FALSE)
      cat(sprintf("  [LMM: total_bc ~ class + (1|region)]\n"))
      cat(sprintf("    Fixed effect (class): est=%.3f, t=%.3f, p=%.4f\n",
                  fixef(m_bc)[2], t_val_bc, lmm_p_bc))
    }, error = function(e) {
      cat(sprintf("    LMM failed: %s\n", e$message))
    })
  }

  # Agreement with original classification (Cohen's kappa)
  d_kappa <- d %>% filter(!is.na(original_H1))
  if (nrow(d_kappa) > 0 && length(unique(d_kappa[[class_var]])) == 2 &&
      length(unique(d_kappa$original_H1)) == 2) {
    kappa_val <- kappa2(data.frame(
      original = as.numeric(d_kappa$original_H1),
      alternative = as.numeric(d_kappa[[class_var]])
    ))
    cat(sprintf("\n  [Agreement with original classification]\n"))
    cat(sprintf("    Cohen's kappa = %.3f (p = %.4f)\n", kappa_val$value, kappa_val$p.value))

    # Confusion matrix
    cm <- table(Original_H1 = d_kappa$original_H1, Alt_H1 = d_kappa[[class_var]])
    cat("    Confusion matrix:\n")
    print(cm)
  }

  return(list(
    label = class_label,
    n_h1 = length(h1_pastoral),
    n_h0 = length(h0_pastoral),
    pastoral_h1_median = median(h1_pastoral),
    pastoral_h0_median = median(h0_pastoral),
    pastoral_mw_p = mw_pastoral$p.value,
    pastoral_r = r_pastoral,
    bc_h1_median = median(h1_bc),
    bc_h0_median = median(h0_bc),
    bc_mw_p = mw_bc$p.value,
    bc_r = r_bc,
    lmm_p_pastoral = lmm_p_pastoral,
    lmm_p_bc = lmm_p_bc
  ))
}

# --- Original classification ---
r0 <- run_comparison(all_dec, "original_H1", "ORIGINAL (2SD exceedance)")
results_list[["original"]] <- r0

# --- Alt 1: Indicator presence ---
r1 <- run_comparison(all_dec, "alt1_H1", "ALT1: Indicator presence (Plantago/Cerealia)")
results_list[["alt1"]] <- r1

# --- Alt 2: Temporal clustering ---
r2 <- run_comparison(all_dec, "alt2_H1", "ALT2: Temporal clustering (±1500yr of Neolithic)")
results_list[["alt2"]] <- r2

# --- Alt 3: Geographic position (subset with coords) ---
if ("alt3_H1" %in% names(all_dec)) {
  r3 <- run_comparison(all_dec, "alt3_H1", "ALT3: Geographic position (±1000yr of local Neolithic)")
  results_list[["alt3"]] <- r3
}

# ============================================================================
# 6. BY-REGION BREAKDOWN
# ============================================================================

cat("\n\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("BY-REGION BREAKDOWN\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")

for (reg in c("UK", "Scandinavia", "Central_Europe")) {
  cat(sprintf("\n\n### %s ###\n", reg))
  d_reg <- all_dec %>% filter(region == reg)

  for (cv in c("original_H1", "alt1_H1", "alt2_H1")) {
    label <- switch(cv,
      "original_H1" = "Original",
      "alt1_H1" = "Alt1 (presence)",
      "alt2_H1" = "Alt2 (temporal)")

    d_sub <- d_reg %>% filter(!is.na(!!sym(cv)))
    n1 <- sum(d_sub[[cv]])
    n0 <- sum(!d_sub[[cv]])

    if (n1 < 2 || n0 < 2) {
      cat(sprintf("  %s: SKIPPED (n1=%d, n0=%d)\n", label, n1, n0))
      next
    }

    h1_p <- d_sub$pastoral_arable_pct[d_sub[[cv]]]
    h0_p <- d_sub$pastoral_arable_pct[!d_sub[[cv]]]
    mw <- wilcox.test(h1_p, h0_p)
    cat(sprintf("  %s: H1 median=%.2f%% vs H0 median=%.2f%%, MW p=%.4f (n=%d vs %d)\n",
                label, median(h1_p), median(h0_p), mw$p.value, n1, n0))
  }

  # Alt3 for sites with coords
  if ("alt3_H1" %in% names(all_dec)) {
    d_sub <- d_reg %>% filter(!is.na(alt3_H1))
    n1 <- sum(d_sub$alt3_H1, na.rm=TRUE)
    n0 <- sum(!d_sub$alt3_H1, na.rm=TRUE)
    if (n1 >= 2 && n0 >= 2) {
      h1_p <- d_sub$pastoral_arable_pct[d_sub$alt3_H1]
      h0_p <- d_sub$pastoral_arable_pct[!d_sub$alt3_H1]
      mw <- wilcox.test(h1_p, h0_p)
      cat(sprintf("  Alt3 (geographic): H1 median=%.2f%% vs H0 median=%.2f%%, MW p=%.4f (n=%d vs %d)\n",
                  median(h1_p), median(h0_p), mw$p.value, n1, n0))
    } else {
      cat(sprintf("  Alt3 (geographic): SKIPPED (n1=%d, n0=%d)\n", n1, n0))
    }
  }
}

# ============================================================================
# 7. SUMMARY TABLE
# ============================================================================

cat("\n\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("SUMMARY COMPARISON\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

cat(sprintf("%-45s %5s %5s %8s %8s %8s %8s\n",
            "Classification", "n_H1", "n_H0",
            "past_p", "past_r", "bc_p", "bc_r"))
cat(paste(rep("-", 90), collapse=""), "\n")

for (key in names(results_list)) {
  r <- results_list[[key]]
  if (!is.null(r)) {
    cat(sprintf("%-45s %5d %5d %8.4f %8.3f %8.4f %8.3f\n",
                r$label, r$n_h1, r$n_h0,
                r$pastoral_mw_p, r$pastoral_r,
                r$bc_mw_p, r$bc_r))
  }
}

# ============================================================================
# 8. WRITE MARKDOWN REPORT
# ============================================================================

sink("shared/phase2_alt_classifications.md")
cat("# Phase 2: Alternative Classifications — Breaking the Circularity\n\n")
cat("## Problem\n\n")
cat("Original classification: H1 = pastoral/arable taxa exceed 2SD threshold.\n")
cat("Finding that H1 sites have higher pastoral BC% is potentially tautological.\n\n")

cat("## Methods\n\n")
cat("Three alternative classifications tested, each independent of pastoral BC%:\n\n")
cat("1. **Alt1 — Indicator presence**: H1 if Plantago lanceolata OR Cerealia-type is PRESENT (>0%)\n")
cat("   - Binary presence vs. quantitative 2SD exceedance = different criterion\n")
cat("   - Uses only the two most unambiguous anthropogenic indicators\n\n")
cat("2. **Alt2 — Temporal clustering**: H1 if signal onset within ±1500yr of regional Neolithic arrival\n")
cat("   - Fully independent of composition — based on timing alone\n")
cat(sprintf("   - Neolithic arrivals: UK=%d, Scandinavia=%d, CE=%d cal BP\n\n",
            neolithic_arrival["UK"], neolithic_arrival["Scandinavia"], neolithic_arrival["Central_Europe"]))
cat("3. **Alt3 — Geographic position**: H1 if onset within ±1000yr of expected local Neolithic\n")
cat(sprintf("   - Linear frontier model: expected_age = %.0f + %.1f×lat + %.1f×lon (R²=%.3f)\n",
            coef(neo_model)[1], coef(neo_model)[2], coef(neo_model)[3],
            summary(neo_model)$r.squared))
cat("   - Only applicable to sites with coordinates\n\n")

cat("## Results\n\n")
cat("### Summary Table\n\n")
cat("| Classification | n_H1 | n_H0 | Pastoral% MW p | Effect r | Total BC MW p | Effect r |\n")
cat("|:---|:---:|:---:|:---:|:---:|:---:|:---:|\n")
for (key in names(results_list)) {
  r <- results_list[[key]]
  if (!is.null(r)) {
    sig_p <- ifelse(r$pastoral_mw_p < 0.001, "***",
             ifelse(r$pastoral_mw_p < 0.01, "**",
             ifelse(r$pastoral_mw_p < 0.05, "*", "ns")))
    sig_bc <- ifelse(r$bc_mw_p < 0.001, "***",
              ifelse(r$bc_mw_p < 0.01, "**",
              ifelse(r$bc_mw_p < 0.05, "*", "ns")))
    cat(sprintf("| %s | %d | %d | %.4f %s | %.3f | %.4f %s | %.3f |\n",
                r$label, r$n_h1, r$n_h0,
                r$pastoral_mw_p, sig_p, r$pastoral_r,
                r$bc_mw_p, sig_bc, r$bc_r))
  }
}

cat("\n### Interpretation\n\n")

# Determine best alternative
best_key <- NULL
best_p <- 1
for (key in c("alt1", "alt2", "alt3")) {
  r <- results_list[[key]]
  if (!is.null(r) && r$pastoral_mw_p < best_p) {
    best_p <- r$pastoral_mw_p
    best_key <- key
  }
}

if (!is.null(best_key)) {
  r_best <- results_list[[best_key]]
  cat(sprintf("**Strongest non-circular evidence**: %s (p=%.4f, r=%.3f for pastoral%%)\n\n",
              r_best$label, r_best$pastoral_mw_p, r_best$pastoral_r))
}

cat("#### Circularity Assessment\n\n")
cat("- **Alt1 (indicator presence)**: Low circularity. Binary presence of Plantago/Cerealia\n")
cat("  is conceptually different from quantitative 2SD exceedance of pastoral+arable group.\n")
cat("  However, both rely on anthropogenic pollen taxa, so partial overlap remains.\n\n")
cat("- **Alt2 (temporal clustering)**: Zero circularity. Classification is based purely on\n")
cat("  timing relative to known archaeological dates. If H1_alt2 sites show higher\n")
cat("  pastoral BC%, this cannot be an artifact of the classification method.\n\n")
cat("- **Alt3 (geographic position)**: Zero circularity. Classification uses only\n")
cat("  coordinates and timing, completely independent of pollen composition.\n")
cat("  Limited by coordinate availability.\n\n")

cat("#### Which alternative provides strongest evidence?\n\n")
cat("Alt2 and Alt3 are the most methodologically clean because they are fully\n")
cat("independent of pollen composition. If either shows significant pastoral%\n")
cat("differences, the finding is definitively non-circular.\n\n")
cat("Alt1 reduces circularity substantially (binary presence vs. quantitative\n")
cat("exceedance) but does not eliminate it entirely.\n\n")

# LMM results
if (has_lme4) {
  cat("### Mixed-Effects Models (LMM)\n\n")
  cat("| Classification | pastoral% LMM p | total_bc LMM p |\n")
  cat("|:---|:---:|:---:|\n")
  for (key in names(results_list)) {
    r <- results_list[[key]]
    if (!is.null(r)) {
      p_p <- ifelse(is.na(r$lmm_p_pastoral), "N/A", sprintf("%.4f", r$lmm_p_pastoral))
      p_b <- ifelse(is.na(r$lmm_p_bc), "N/A", sprintf("%.4f", r$lmm_p_bc))
      cat(sprintf("| %s | %s | %s |\n", r$label, p_p, p_b))
    }
  }
  cat("\n")
}

cat("\n---\n")
cat(sprintf("*Generated: %s*\n", Sys.time()))
sink()

cat("\n\nResults written to: shared/phase2_alt_classifications.md\n")
cat("Done.\n")
