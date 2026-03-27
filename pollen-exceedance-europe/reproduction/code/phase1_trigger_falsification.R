#!/usr/bin/env Rscript
# Phase 1: Trigger Falsification Test
# Compare pastoral+arable BC% and total BC between H1 and H3 sites
# Mann-Whitney U test with rank-biserial effect size

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

cat("=" ,rep("=", 70), "\n", sep="")
cat("PHASE 1: TRIGGER FALSIFICATION TEST\n")
cat("H1 (anthropogenic) vs H3 (natural) — Pastoral BC% and Total BC\n")
cat("=" ,rep("=", 70), "\n\n", sep="")

# --- Helper functions ---

run_mw_test <- function(x_h1, x_h3, label) {
  # Mann-Whitney U test with rank-biserial effect size
  n1 <- length(x_h1)
  n2 <- length(x_h3)

  if (n1 < 2 || n2 < 2) {
    cat(sprintf("  %s: Insufficient data (n_H1=%d, n_H3=%d)\n", label, n1, n2))
    return(list(p = NA, r = NA, n1 = n1, n2 = n2,
                med_h1 = median(x_h1), med_h3 = median(x_h3)))
  }

  wt <- wilcox.test(x_h1, x_h3, exact = FALSE, correct = TRUE)
  U <- wt$statistic  # This is the W statistic = U for x_h1
  r_effect <- 1 - 2 * U / (n1 * n2)  # rank-biserial correlation

  cat(sprintf("  %s:\n", label))
  cat(sprintf("    H1: n=%d, median=%.3f, IQR=[%.3f, %.3f]\n",
              n1, median(x_h1), quantile(x_h1, 0.25), quantile(x_h1, 0.75)))
  cat(sprintf("    H3: n=%d, median=%.3f, IQR=[%.3f, %.3f]\n",
              n2, median(x_h3), quantile(x_h3, 0.25), quantile(x_h3, 0.75)))
  cat(sprintf("    Mann-Whitney W=%.1f, p=%.4f\n", U, wt$p.value))
  cat(sprintf("    Rank-biserial r=%.3f", r_effect))
  if (abs(r_effect) < 0.1) cat(" (negligible)")
  else if (abs(r_effect) < 0.3) cat(" (small)")
  else if (abs(r_effect) < 0.5) cat(" (medium)")
  else cat(" (large)")
  cat("\n\n")

  return(list(p = wt$p.value, r = r_effect, n1 = n1, n2 = n2,
              med_h1 = median(x_h1), med_h3 = median(x_h3),
              iqr_h1 = paste0("[", round(quantile(x_h1, 0.25), 3), ", ", round(quantile(x_h1, 0.75), 3), "]"),
              iqr_h3 = paste0("[", round(quantile(x_h3, 0.25), 3), ", ", round(quantile(x_h3, 0.75), 3), "]")))
}

# ============================================================
# 1. UK DATA
# ============================================================
cat("\n--- UK (Britain & Ireland) ---\n\n")

uk <- read.csv("/home/ayu/archeco/shared/signal_driver_investigation.csv",
               stringsAsFactors = FALSE)

# Get pastoral BC contribution per site
uk_pastoral <- uk %>%
  filter(category2 == "pastoral") %>%
  select(sitename, datasetid, pastoral_bc = cat_contribution, hypothesis)

# Get total_bc per site (same for all rows of a site)
uk_total <- uk %>%
  distinct(sitename, datasetid, total_bc, hypothesis)

# Merge
uk_site <- uk_total %>%
  left_join(uk_pastoral %>% select(sitename, datasetid, pastoral_bc),
            by = c("sitename", "datasetid")) %>%
  mutate(pastoral_bc = ifelse(is.na(pastoral_bc), 0, pastoral_bc),
         pastoral_pct = pastoral_bc / total_bc * 100)

cat(sprintf("  Total UK sites: %d (H1: %d, H3: %d, indeterminate: %d)\n",
            nrow(uk_site),
            sum(uk_site$hypothesis == "H1_anthropogenic"),
            sum(grepl("^H3_", uk_site$hypothesis)),
            sum(uk_site$hypothesis == "indeterminate")))
cat("  H3 subtypes:", paste(unique(uk_site$hypothesis[grepl("^H3_", uk_site$hypothesis)]), collapse=", "), "\n")

uk_h1 <- uk_site %>% filter(hypothesis == "H1_anthropogenic")
uk_h3 <- uk_site %>% filter(grepl("^H3_", hypothesis))

uk_pastoral_res <- run_mw_test(uk_h1$pastoral_pct, uk_h3$pastoral_pct,
                                "Pastoral BC%")
uk_bc_res <- run_mw_test(uk_h1$total_bc, uk_h3$total_bc,
                          "Total BC magnitude")


# ============================================================
# 2. SCANDINAVIA DATA
# ============================================================
cat("\n--- Scandinavia ---\n\n")

scan <- read.csv("/home/ayu/archeco/shared/signal_decomposition_scandinavia.csv",
                 stringsAsFactors = FALSE)

# Remove duplicates (same datasetid appearing twice)
scan <- scan %>% distinct(datasetid, .keep_all = TRUE)

cat(sprintf("  Total Scandinavia sites: %d (H1: %d, H3: %d)\n",
            nrow(scan),
            sum(scan$classification == "H1_anthropogenic"),
            sum(scan$classification == "H3_natural")))

scan_h1 <- scan %>% filter(classification == "H1_anthropogenic")
scan_h3 <- scan %>% filter(classification == "H3_natural")

# bc_anthro_pct is the pastoral+arable combined %
scan_pastoral_res <- run_mw_test(scan_h1$bc_anthro_pct, scan_h3$bc_anthro_pct,
                                  "Pastoral+Arable BC%")
scan_bc_res <- run_mw_test(scan_h1$bc_total, scan_h3$bc_total,
                            "Total BC magnitude")


# ============================================================
# 3. CENTRAL EUROPE DATA
# ============================================================
cat("\n--- Central Europe ---\n\n")

ce <- read.csv("/home/ayu/archeco/shared/signal_decomposition_central_europe.csv",
               stringsAsFactors = FALSE)

cat(sprintf("  Total CE sites: %d (H1: %d, H3: %d)\n",
            nrow(ce),
            sum(ce$classification == "H1_anthropogenic"),
            sum(ce$classification == "H3_natural")))

ce_h1 <- ce %>% filter(classification == "H1_anthropogenic")
ce_h3 <- ce %>% filter(classification == "H3_natural")

ce_pastoral_res <- run_mw_test(ce_h1$pastoral_arable_pct, ce_h3$pastoral_arable_pct,
                                "Pastoral+Arable BC%")
ce_bc_res <- run_mw_test(ce_h1$total_bc, ce_h3$total_bc,
                          "Total BC magnitude")


# ============================================================
# 4. COMBINED (ALL REGIONS)
# ============================================================
cat("\n--- COMBINED (All Regions) ---\n\n")

# Standardize and combine
all_data <- bind_rows(
  uk_site %>%
    mutate(classification = ifelse(grepl("^H1_", hypothesis), "H1_anthropogenic",
                            ifelse(grepl("^H3_", hypothesis), "H3_natural", hypothesis))) %>%
    filter(classification %in% c("H1_anthropogenic", "H3_natural")) %>%
    transmute(region = "UK", sitename, datasetid,
              pastoral_arable_pct = pastoral_pct,
              total_bc,
              classification),
  scan %>%
    transmute(region = "Scandinavia", sitename, datasetid,
              pastoral_arable_pct = bc_anthro_pct,
              total_bc = bc_total,
              classification),
  ce %>%
    transmute(region = "Central_Europe", sitename, datasetid,
              pastoral_arable_pct,
              total_bc,
              classification)
)

cat(sprintf("  Total combined sites: %d (H1: %d, H3: %d)\n",
            nrow(all_data),
            sum(all_data$classification == "H1_anthropogenic"),
            sum(all_data$classification == "H3_natural")))

all_h1 <- all_data %>% filter(classification == "H1_anthropogenic")
all_h3 <- all_data %>% filter(classification == "H3_natural")

all_pastoral_res <- run_mw_test(all_h1$pastoral_arable_pct, all_h3$pastoral_arable_pct,
                                 "Pastoral+Arable BC%")
all_bc_res <- run_mw_test(all_h1$total_bc, all_h3$total_bc,
                           "Total BC magnitude")


# ============================================================
# 5. SUMMARY TABLE
# ============================================================
cat("\n")
cat("=" ,rep("=", 70), "\n", sep="")
cat("SUMMARY TABLE\n")
cat("=" ,rep("=", 70), "\n\n")

summary_rows <- list(
  list("UK", "Pastoral BC%", uk_pastoral_res),
  list("UK", "Total BC", uk_bc_res),
  list("Scandinavia", "Pastoral+Arable BC%", scan_pastoral_res),
  list("Scandinavia", "Total BC", scan_bc_res),
  list("Central Europe", "Pastoral+Arable BC%", ce_pastoral_res),
  list("Central Europe", "Total BC", ce_bc_res),
  list("COMBINED", "Pastoral+Arable BC%", all_pastoral_res),
  list("COMBINED", "Total BC", all_bc_res)
)

cat(sprintf("%-18s %-22s %5s %5s %8s %8s %8s %8s\n",
            "Region", "Metric", "n_H1", "n_H3", "med_H1", "med_H3", "p", "r_effect"))
cat(paste(rep("-", 95), collapse=""), "\n")

for (row in summary_rows) {
  res <- row[[3]]
  sig <- ifelse(is.na(res$p), "  ", ifelse(res$p < 0.001, "***",
                ifelse(res$p < 0.01, "** ",
                ifelse(res$p < 0.05, "*  ", "ns "))))
  cat(sprintf("%-18s %-22s %5d %5d %8.3f %8.3f %8.4f %8.3f %s\n",
              row[[1]], row[[2]], res$n1, res$n2,
              res$med_h1, res$med_h3, res$p, res$r, sig))
}

cat("\nSignificance: *** p<0.001, ** p<0.01, * p<0.05, ns p>=0.05\n")
cat("Effect size r: negligible<0.1, small<0.3, medium<0.5, large>=0.5\n")

# ============================================================
# 6. CONCLUSION
# ============================================================
cat("\n")
cat("=" ,rep("=", 70), "\n", sep="")
cat("CONCLUSION\n")
cat("=" ,rep("=", 70), "\n\n")

# Check combined pastoral result
if (!is.na(all_pastoral_res$p) && all_pastoral_res$p < 0.05 && all_pastoral_res$med_h1 > all_pastoral_res$med_h3) {
  cat("TRIGGER MODEL: SUPPORTED\n")
  cat("H1 sites show significantly higher pastoral+arable BC% than H3 sites.\n")
  cat(sprintf("Combined: H1 median=%.2f%% vs H3 median=%.2f%%, p=%.4f, r=%.3f\n",
              all_pastoral_res$med_h1, all_pastoral_res$med_h3,
              all_pastoral_res$p, all_pastoral_res$r))
  trigger_conclusion <- "SUPPORTED"
} else if (!is.na(all_pastoral_res$p) && all_pastoral_res$p >= 0.05) {
  cat("TRIGGER MODEL: FALSIFIED (or unsupported)\n")
  cat("No significant difference in pastoral+arable BC% between H1 and H3 sites.\n")
  cat(sprintf("Combined: H1 median=%.2f%% vs H3 median=%.2f%%, p=%.4f, r=%.3f\n",
              all_pastoral_res$med_h1, all_pastoral_res$med_h3,
              all_pastoral_res$p, all_pastoral_res$r))
  trigger_conclusion <- "FALSIFIED"
} else {
  cat("TRIGGER MODEL: INCONCLUSIVE\n")
  trigger_conclusion <- "INCONCLUSIVE"
}


# ============================================================
# 7. WRITE MARKDOWN REPORT
# ============================================================

report <- c(
  "# Phase 1: Trigger Falsification Test Results",
  "",
  "## Falsification Criterion",
  "If H1 sites' pastoral+arable BC% is NOT significantly higher than H3 sites',",
  "the trigger model (small pastoral disturbance initiates cascading forest reorganization) is unsupported.",
  "",
  "## Methods",
  "- Mann-Whitney U test (two-sided, continuity-corrected)",
  "- Effect size: rank-biserial correlation r = 1 - 2U/(n1*n2)",
  "- Metrics: (1) Pastoral+arable BC% of total Bray-Curtis dissimilarity, (2) Total BC magnitude",
  "",
  "## Results by Region",
  "",
  "### UK (Britain & Ireland)",
  sprintf("- **Pastoral BC%%**: H1 median=%.2f%% (n=%d) vs H3 median=%.2f%% (n=%d), p=%.4f, r=%.3f",
          uk_pastoral_res$med_h1, uk_pastoral_res$n1, uk_pastoral_res$med_h3, uk_pastoral_res$n2,
          uk_pastoral_res$p, uk_pastoral_res$r),
  sprintf("- **Total BC**: H1 median=%.3f (n=%d) vs H3 median=%.3f (n=%d), p=%.4f, r=%.3f",
          uk_bc_res$med_h1, uk_bc_res$n1, uk_bc_res$med_h3, uk_bc_res$n2,
          uk_bc_res$p, uk_bc_res$r),
  "",
  "### Scandinavia",
  sprintf("- **Pastoral+Arable BC%%**: H1 median=%.2f%% (n=%d) vs H3 median=%.2f%% (n=%d), p=%.4f, r=%.3f",
          scan_pastoral_res$med_h1, scan_pastoral_res$n1, scan_pastoral_res$med_h3, scan_pastoral_res$n2,
          scan_pastoral_res$p, scan_pastoral_res$r),
  sprintf("- **Total BC**: H1 median=%.3f (n=%d) vs H3 median=%.3f (n=%d), p=%.4f, r=%.3f",
          scan_bc_res$med_h1, scan_bc_res$n1, scan_bc_res$med_h3, scan_bc_res$n2,
          scan_bc_res$p, scan_bc_res$r),
  "",
  "### Central Europe",
  sprintf("- **Pastoral+Arable BC%%**: H1 median=%.2f%% (n=%d) vs H3 median=%.2f%% (n=%d), p=%.4f, r=%.3f",
          ce_pastoral_res$med_h1, ce_pastoral_res$n1, ce_pastoral_res$med_h3, ce_pastoral_res$n2,
          ce_pastoral_res$p, ce_pastoral_res$r),
  sprintf("- **Total BC**: H1 median=%.3f (n=%d) vs H3 median=%.3f (n=%d), p=%.4f, r=%.3f",
          ce_bc_res$med_h1, ce_bc_res$n1, ce_bc_res$med_h3, ce_bc_res$n2,
          ce_bc_res$p, ce_bc_res$r),
  "",
  "### Combined (All Regions)",
  sprintf("- **Pastoral+Arable BC%%**: H1 median=%.2f%% (n=%d) vs H3 median=%.2f%% (n=%d), p=%.4f, r=%.3f",
          all_pastoral_res$med_h1, all_pastoral_res$n1, all_pastoral_res$med_h3, all_pastoral_res$n2,
          all_pastoral_res$p, all_pastoral_res$r),
  sprintf("- **Total BC**: H1 median=%.3f (n=%d) vs H3 median=%.3f (n=%d), p=%.4f, r=%.3f",
          all_bc_res$med_h1, all_bc_res$n1, all_bc_res$med_h3, all_bc_res$n2,
          all_bc_res$p, all_bc_res$r),
  "",
  "## Summary Table",
  "",
  "| Region | Metric | n_H1 | n_H3 | Median H1 | Median H3 | p-value | Effect r | Sig |",
  "|--------|--------|------|------|-----------|-----------|---------|----------|-----|"
)

for (row in summary_rows) {
  res <- row[[3]]
  sig <- ifelse(is.na(res$p), "-", ifelse(res$p < 0.001, "***",
                ifelse(res$p < 0.01, "**",
                ifelse(res$p < 0.05, "*", "ns"))))
  report <- c(report,
    sprintf("| %s | %s | %d | %d | %.3f | %.3f | %.4f | %.3f | %s |",
            row[[1]], row[[2]], res$n1, res$n2,
            res$med_h1, res$med_h3, res$p, res$r, sig))
}

report <- c(report,
  "",
  "## Conclusion",
  "",
  sprintf("**Trigger model: %s**", trigger_conclusion),
  ""
)

if (trigger_conclusion == "SUPPORTED") {
  report <- c(report,
    "H1 (anthropogenic) sites show significantly higher pastoral+arable BC% than H3 (natural) sites.",
    "This is consistent with the trigger model: pastoral activity contributes more to vegetation",
    "change at anthropogenic sites, supporting the interpretation that pastoral disturbance",
    "initiates cascading forest reorganization.",
    "",
    "Key observations:",
    sprintf("- Combined effect: H1 median %.2f%% vs H3 median %.2f%% (p=%.4f, r=%.3f)",
            all_pastoral_res$med_h1, all_pastoral_res$med_h3,
            all_pastoral_res$p, all_pastoral_res$r)
  )
} else {
  report <- c(report,
    "No significant difference in pastoral+arable BC% between H1 and H3 sites.",
    "The trigger model lacks empirical support from this falsification test."
  )
}

writeLines(report, "/home/ayu/archeco/shared/phase1_trigger_test_results.md")
cat("\nReport saved to: shared/phase1_trigger_test_results.md\n")
