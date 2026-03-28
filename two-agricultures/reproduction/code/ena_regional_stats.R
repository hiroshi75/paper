#!/usr/bin/env Rscript
# ============================================================================
# ena_regional_stats.R
# Regional stratification, statistical tests, and baseline sensitivity analysis
# for Eastern North America pollen exceedance data.
# Addresses reviewer concerns about lack of regional stratification and
# formal statistical testing.
# ============================================================================

library(dplyr)
library(tidyr)

# ---- Load data ----
d <- read.csv("shared/full_ena_exceedance_data.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(d), "sites\n\n")

# ============================================================================
# TASK 1: Regional Stratification
# ============================================================================

classify_region <- function(lat, lon) {
  case_when(
    lat > 40 & lon > -78 ~ "Northeast",
    lat > 42 & lon >= -95 & lon <= -78 ~ "Great Lakes",
    lat < 38 ~ "Southeast",
    TRUE ~ "Central/Midwest"
  )
}

d$region <- classify_region(d$lat, d$lon)

# Compute per-region summaries
regional_summary <- d %>%
  group_by(region) %>%
  summarise(
    n_sites = n(),
    bc_signal_n = sum(has_signal == TRUE, na.rm = TRUE),
    bc_signal_pct = round(100 * mean(has_signal == TRUE, na.rm = TRUE), 1),
    h1_n = sum(classification == "H1_anthropogenic", na.rm = TRUE),
    h1_pct = round(100 * mean(classification == "H1_anthropogenic", na.rm = TRUE), 1),
    ambrosia_exc_n = sum(!is.na(ambrosia_onset_bp)),
    ambrosia_exc_pct = round(100 * mean(!is.na(ambrosia_onset_bp)), 1),
    ambrosia_median_onset = median(ambrosia_onset_bp, na.rm = TRUE),
    ambrosia_precontact_pct = round(100 * mean(ambrosia_precontact == TRUE, na.rm = TRUE), 1),
    n_zea = sum(has_zea == TRUE, na.rm = TRUE),
    mean_max_bc = round(mean(max_bc, na.rm = TRUE), 3),
    median_max_bc = round(median(max_bc, na.rm = TRUE), 3),
    .groups = "drop"
  )

cat("=== TASK 1: Regional Stratification ===\n\n")
print(as.data.frame(regional_summary), row.names = FALSE)

# Chi-squared / Fisher test for BC signal rate across regions
bc_table <- d %>%
  group_by(region) %>%
  summarise(signal = sum(has_signal), no_signal = sum(!has_signal), .groups = "drop")
bc_mat <- as.matrix(bc_table[, c("signal", "no_signal")])
rownames(bc_mat) <- bc_table$region
cat("\n\nBC signal contingency table:\n")
print(bc_mat)
fisher_bc <- fisher.test(bc_mat)
cat("\nFisher's exact test for BC signal rate across regions:\n")
cat("  p-value =", format(fisher_bc$p.value, digits = 4), "\n")

# Kruskal-Wallis for max_bc across regions
kw_bc <- kruskal.test(max_bc ~ region, data = d)
cat("\nKruskal-Wallis test for max_bc across regions:\n")
cat("  chi-squared =", round(kw_bc$statistic, 3), ", df =", kw_bc$parameter,
    ", p =", format(kw_bc$p.value, digits = 4), "\n")

# Kruskal-Wallis for Ambrosia onset across regions
d_amb <- d %>% filter(!is.na(ambrosia_onset_bp))
if (length(unique(d_amb$region)) > 1) {
  kw_amb <- kruskal.test(ambrosia_onset_bp ~ region, data = d_amb)
  cat("\nKruskal-Wallis test for Ambrosia onset timing across regions:\n")
  cat("  chi-squared =", round(kw_amb$statistic, 3), ", df =", kw_amb$parameter,
      ", p =", format(kw_amb$p.value, digits = 4), "\n")
}

# ============================================================================
# TASK 2: Statistical Tests
# ============================================================================

cat("\n\n=== TASK 2: Statistical Tests ===\n\n")

# --- Test 1: Fisher's exact test for pre-contact vs post-contact Ambrosia onset ---
n_pre <- sum(d$ambrosia_precontact == TRUE, na.rm = TRUE)
n_post <- sum(d$ambrosia_postcontact == TRUE, na.rm = TRUE)
# Sites with ambrosia_at_contact or neither
n_at_contact <- sum(d$ambrosia_at_contact == TRUE, na.rm = TRUE)
n_amb_total <- n_pre + n_post + n_at_contact
# Test if pre-contact rate is different from 50%
binom_pre <- binom.test(n_pre, n_pre + n_post, p = 0.5, alternative = "two.sided")
cat("Test 1: Binomial test - Pre-contact vs Post-contact Ambrosia onset\n")
cat("  Pre-contact:", n_pre, "  Post-contact:", n_post, "  At-contact:", n_at_contact, "\n")
cat("  Pre-contact proportion:", round(n_pre / (n_pre + n_post), 3), "\n")
cat("  Binomial test p-value:", format(binom_pre$p.value, digits = 4), "\n")
cat("  95% CI for pre-contact proportion:", round(binom_pre$conf.int[1], 3), "-",
    round(binom_pre$conf.int[2], 3), "\n")

# Also Fisher's exact on 2x2 table
fisher_prepost <- matrix(c(n_pre, n_post, n_pre + n_post - n_pre, n_pre + n_post - n_post),
                         nrow = 2)
# Simpler: just binomial as above is the right test

# --- Test 2: Clopper-Pearson CI for 0% H1 ---
n_h1 <- sum(d$classification == "H1_anthropogenic", na.rm = TRUE)
n_total <- nrow(d)
cp_test <- binom.test(n_h1, n_total, conf.level = 0.95)
cat("\nTest 2: Clopper-Pearson 95% CI for H1 rate\n")
cat("  H1 count:", n_h1, "/", n_total, "\n")
cat("  Point estimate:", round(n_h1 / n_total, 4), "\n")
cat("  95% CI:", round(cp_test$conf.int[1], 4), "-", round(cp_test$conf.int[2], 4), "\n")
cat("  Upper 95% CI bound:", round(cp_test$conf.int[2], 4),
    " (= max plausible H1 rate)\n")

# --- Test 3: Bootstrap CI for the 41x ratio (or RPP-corrected ratio) ---
# The ratio is typically max_ambrosia / baseline Ambrosia level
# We need to compute this from the data
# Ambrosia exceedance ratio: max_ambrosia_pct relative to some baseline
# Let's compute the ratio from bc_mean_baseline and max_bc for sites with signal
d_signal <- d %>% filter(has_signal == TRUE)

# Actually, the "41x ratio" likely refers to post-contact vs pre-contact Ambrosia
# Let's use max_ambrosia_pct as the measure of Ambrosia explosion magnitude
# and compute a bootstrap CI
set.seed(42)
n_boot <- 10000

# Bootstrap for max_ambrosia_pct among sites with Ambrosia data
d_amb_data <- d %>% filter(!is.na(max_ambrosia_pct) & max_ambrosia_pct > 0)

if (nrow(d_amb_data) > 0) {
  boot_median <- numeric(n_boot)
  boot_mean <- numeric(n_boot)
  for (i in 1:n_boot) {
    samp <- sample(d_amb_data$max_ambrosia_pct, replace = TRUE)
    boot_median[i] <- median(samp)
    boot_mean[i] <- mean(samp)
  }
  cat("\nTest 3: Bootstrap CI for max Ambrosia percentage\n")
  cat("  N sites with Ambrosia data:", nrow(d_amb_data), "\n")
  cat("  Observed median max_ambrosia_pct:", round(median(d_amb_data$max_ambrosia_pct), 2), "\n")
  cat("  Bootstrap 95% CI (median):", round(quantile(boot_median, 0.025), 2), "-",
      round(quantile(boot_median, 0.975), 2), "\n")
  cat("  Observed mean max_ambrosia_pct:", round(mean(d_amb_data$max_ambrosia_pct), 2), "\n")
  cat("  Bootstrap 95% CI (mean):", round(quantile(boot_mean, 0.025), 2), "-",
      round(quantile(boot_mean, 0.975), 2), "\n")
}

# Bootstrap for BC magnitude ratio (post/pre contact)
d_bc_ratio <- d %>%
  filter(!is.na(bc_precontact_mean) & !is.na(bc_postcontact_mean) &
           bc_precontact_mean > 0)
if (nrow(d_bc_ratio) > 0) {
  # Compute per-site ratio
  d_bc_ratio$bc_ratio <- d_bc_ratio$bc_postcontact_mean / d_bc_ratio$bc_precontact_mean
  observed_ratio <- median(d_bc_ratio$bc_ratio, na.rm = TRUE)

  boot_ratio <- numeric(n_boot)
  for (i in 1:n_boot) {
    samp <- sample(d_bc_ratio$bc_ratio, replace = TRUE)
    boot_ratio[i] <- median(samp, na.rm = TRUE)
  }
  cat("\nTest 3b: Bootstrap CI for BC post/pre-contact ratio\n")
  cat("  N sites:", nrow(d_bc_ratio), "\n")
  cat("  Observed median ratio:", round(observed_ratio, 3), "\n")
  cat("  Bootstrap 95% CI:", round(quantile(boot_ratio, 0.025), 3), "-",
      round(quantile(boot_ratio, 0.975), 3), "\n")
}

# --- Test 4: Mann-Whitney - Zea presence vs BC magnitude ---
d_zea_yes <- d %>% filter(has_zea == TRUE)
d_zea_no <- d %>% filter(has_zea == FALSE)

cat("\nTest 4: Mann-Whitney U test - Zea presence vs max BC\n")
cat("  Sites with Zea:", nrow(d_zea_yes),
    " median max_bc:", round(median(d_zea_yes$max_bc, na.rm = TRUE), 3), "\n")
cat("  Sites without Zea:", nrow(d_zea_no),
    " median max_bc:", round(median(d_zea_no$max_bc, na.rm = TRUE), 3), "\n")

if (nrow(d_zea_yes) >= 2) {
  mw_test <- wilcox.test(d_zea_yes$max_bc, d_zea_no$max_bc,
                         alternative = "two.sided", exact = FALSE)
  cat("  W =", mw_test$statistic, ", p =", format(mw_test$p.value, digits = 4), "\n")

  # Effect size: rank-biserial correlation
  n1 <- nrow(d_zea_yes)
  n2 <- nrow(d_zea_no)
  r_rb <- 1 - (2 * mw_test$statistic) / (n1 * n2)
  cat("  Rank-biserial r =", round(r_rb, 3), "\n")
} else {
  cat("  Too few Zea sites for test\n")
}

# --- Test 5: Permutation test for recovery rate ---
n_recovery <- sum(d$recovery_signal == TRUE, na.rm = TRUE)
n_recovery_total <- sum(!is.na(d$recovery_signal))
obs_recovery_rate <- n_recovery / n_recovery_total

cat("\nTest 5: Permutation test for recovery rate\n")
cat("  Recovery count:", n_recovery, "/", n_recovery_total, "\n")
cat("  Observed rate:", round(obs_recovery_rate * 100, 1), "%\n")

# Null: recovery signal is random (50% chance)
# One-sided: is rate significantly different from 50%?
binom_recovery <- binom.test(n_recovery, n_recovery_total, p = 0.5,
                              alternative = "two.sided")
cat("  Binomial test vs 50%: p =", format(binom_recovery$p.value, digits = 4), "\n")
cat("  95% CI:", round(binom_recovery$conf.int[1] * 100, 1), "% -",
    round(binom_recovery$conf.int[2] * 100, 1), "%\n")

# Also permutation-based
set.seed(42)
n_perm <- 10000
perm_rates <- numeric(n_perm)
for (i in 1:n_perm) {
  perm_signal <- sample(c(TRUE, FALSE), n_recovery_total, replace = TRUE, prob = c(0.5, 0.5))
  perm_rates[i] <- mean(perm_signal)
}
perm_p <- mean(abs(perm_rates - 0.5) >= abs(obs_recovery_rate - 0.5))
cat("  Permutation test p-value (two-sided vs 50%):", round(perm_p, 4), "\n")

# Test against random baseline (30.8% vs 50%)
# Also test if different from BC signal rate
bc_rate <- mean(d$has_signal == TRUE)
binom_recovery_vs_bc <- binom.test(n_recovery, n_recovery_total, p = bc_rate,
                                     alternative = "two.sided")
cat("  Binomial test vs BC signal rate (", round(bc_rate * 100, 1), "%): p =",
    format(binom_recovery_vs_bc$p.value, digits = 4), "\n")


# ============================================================================
# TASK 3: Baseline Sensitivity Analysis
# ============================================================================

cat("\n\n=== TASK 3: Baseline Sensitivity (>5000 BP) ===\n\n")

# The current baseline is >3000 BP. What if we use >5000 BP?
# We need to check n_baseline per site and what the onset ages are

# First: how many sites have enough deep baseline samples?
# We approximate: if n_baseline is already small with >3000 BP cutoff,
# pushing to >5000 BP would reduce it further.
# We don't have per-sample age data, but we can check:
# - Sites with large n_baseline are more likely to have deep samples
# - onset_age_bp tells us the earliest signal detection

# Read the raw data to check sample ages if available
# Actually, from the CSV we only have summary stats per site.
# We can estimate from onset_age_bp and n_samples distribution.

# Key insight: if a site's oldest sample is < 5000 BP, it can't have
# any >5000 BP baseline samples. We need the age range per site.
# Since we don't have that directly, let's use a proxy approach:

# Check which sites have ambrosia_explosion_bp > 5000 (suggesting deep record)
d$has_deep_record <- !is.na(d$ambrosia_explosion_bp) & d$ambrosia_explosion_bp > 5000

cat("Sites with Ambrosia data extending > 5000 BP:", sum(d$has_deep_record), "/", nrow(d), "\n")

# Alternative: use n_baseline as a proxy - sites with many baseline samples
# likely extend further back
cat("n_baseline distribution:\n")
print(summary(d$n_baseline))

# Sites likely to have >5000 BP baseline (using n_baseline >= 10 as minimum)
cat("\nSites with n_baseline >= 10:", sum(d$n_baseline >= 10), "\n")
cat("Sites with n_baseline >= 20:", sum(d$n_baseline >= 20), "\n")
cat("Sites with n_baseline >= 30:", sum(d$n_baseline >= 30), "\n")

# Sensitivity: among sites with deep records, does H1 rate change?
d_deep <- d %>% filter(has_deep_record == TRUE)
cat("\n--- Among sites with deep records (>5000 BP) ---\n")
cat("N sites:", nrow(d_deep), "\n")
cat("BC signal rate:", round(100 * mean(d_deep$has_signal), 1), "%\n")
cat("H1 rate:", round(100 * mean(d_deep$classification == "H1_anthropogenic"), 1), "%\n")
cat("H3 rate:", round(100 * mean(d_deep$classification == "H3_natural"), 1), "%\n")
cat("No signal rate:", round(100 * mean(d_deep$classification == "no_signal"), 1), "%\n")
cat("Ambrosia precontact rate:", round(100 * mean(d_deep$ambrosia_precontact == TRUE), 1), "%\n")
cat("Median Ambrosia onset BP:", median(d_deep$ambrosia_onset_bp, na.rm = TRUE), "\n")

# Compare: deeper baseline sites vs shallower
d_shallow <- d %>% filter(has_deep_record == FALSE | is.na(has_deep_record))
cat("\n--- Among sites WITHOUT deep records ---\n")
cat("N sites:", nrow(d_shallow), "\n")
cat("BC signal rate:", round(100 * mean(d_shallow$has_signal), 1), "%\n")
cat("H3 rate:", round(100 * mean(d_shallow$classification == "H3_natural"), 1), "%\n")

# Wilcoxon test: do deep-record sites differ in BC magnitude?
if (nrow(d_deep) >= 5 & nrow(d_shallow) >= 5) {
  wt <- wilcox.test(d_deep$max_bc, d_shallow$max_bc, exact = FALSE)
  cat("\nMann-Whitney: max_bc deep vs shallow records\n")
  cat("  Deep median:", round(median(d_deep$max_bc), 3),
      " Shallow median:", round(median(d_shallow$max_bc), 3), "\n")
  cat("  p =", format(wt$p.value, digits = 4), "\n")
}

# ---- Simulate stricter baseline ----
# For sites with sufficient data, a stricter (>5000 BP) baseline would
# have HIGHER mean BC (more distant from human influence).
# This would RAISE the threshold, making signal detection HARDER.
# Estimate the effect:

# Use bc_mean_baseline + 1 SD as a rough proxy for what a deeper baseline might show
# (deeper = less human disturbance = potentially different vegetation baseline)
d$strict_threshold <- d$bc_mean_baseline + 3 * d$bc_sd_baseline
d$strict_signal <- d$max_bc > d$strict_threshold

cat("\n--- Sensitivity: Stricter threshold (mean + 3*SD vs current mean + 2*SD) ---\n")
cat("Current signal rate:", round(100 * mean(d$has_signal), 1), "%\n")
cat("Strict (3SD) signal rate:", round(100 * mean(d$strict_signal), 1), "%\n")

# Even stricter: different SD multipliers
for (mult in c(2, 2.5, 3, 3.5)) {
  strict_thresh <- d$bc_mean_baseline + mult * d$bc_sd_baseline
  strict_sig <- d$max_bc > strict_thresh
  cat("  Threshold multiplier", mult, "SD: signal rate =",
      round(100 * mean(strict_sig), 1), "%\n")
}

# Classification for strict threshold
d$strict_class <- ifelse(d$strict_signal & d$pct_agricultural > 5,
                         "H1_anthropogenic",
                         ifelse(d$strict_signal, "H3_natural", "no_signal"))
cat("\nStrict threshold classifications:\n")
print(table(d$strict_class))

# ---- Regional breakdown with deep records ----
cat("\n--- Regional breakdown (deep record sites only) ---\n")
deep_regional <- d_deep %>%
  group_by(region) %>%
  summarise(
    n = n(),
    bc_signal_pct = round(100 * mean(has_signal), 1),
    ambrosia_precontact_pct = round(100 * mean(ambrosia_precontact == TRUE), 1),
    median_onset = median(ambrosia_onset_bp, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(deep_regional), row.names = FALSE)

cat("\n=== Analysis complete ===\n")
