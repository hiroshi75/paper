#!/usr/bin/env Rscript
# Phase 2: Lag Investigation — Distribution of pastoral exceedance lags
# relative to Neolithic arrival, across UK, Scandinavia, and Central Europe.
#
# Key question: Is the ~3000yr lag between Neolithic arrival and detectable
# pastoral pollen signals a methodological artefact or a real archaeological pattern?

library(dplyr)

# Helper: read CSV (base R, no readr dependency)
read_csv_base <- function(path) {
  read.csv(path, stringsAsFactors = FALSE)
}
# Alias for compatibility
read_csv <- function(path, ...) read_csv_base(path)
write_csv <- function(df, path) write.csv(df, path, row.names = FALSE)

cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("PHASE 2: PASTORAL LAG INVESTIGATION\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n\n")

# --- Neolithic arrival dates (cal BP) ---
neo_uk <- 6000
neo_scand <- 5800
neo_ce <- 7500

# =====================================================================
# ANALYSIS 1: Distribution of functional onset lags by region
# =====================================================================

cat("ANALYSIS 1: DISTRIBUTION OF PASTORAL EXCEEDANCE LAGS\n")
cat("-" |> rep(50) |> paste(collapse=""), "\n\n")

# --- UK ---
uk <- read_csv("/home/ayu/archeco/shared/composite_indices_analysis.csv")
uk_lat <- read_csv("/home/ayu/archeco/shared/signal_decomposition_uk_corrected.csv",
                   show_col_types = FALSE) %>%
  select(datasetid, lat, lon) %>%
  distinct()

# pastoral_first_change_age is the column for UK
uk_pastoral <- uk %>%
  filter(!is.na(pastoral_first_change_age)) %>%
  left_join(uk_lat, by = "datasetid") %>%
  mutate(lag = pastoral_first_change_age - neo_uk,
         region = "UK")

cat(sprintf("UK: %d sites with pastoral exceedance (out of %d total)\n",
            nrow(uk_pastoral), nrow(uk)))

# --- Scandinavia ---
scand <- read_csv("/home/ayu/archeco/shared/signal_decomposition_scandinavia.csv",
                  show_col_types = FALSE)
scand_lat <- read_csv("/home/ayu/archeco/shared/signal_phase_scandinavia.csv",
                      show_col_types = FALSE) %>%
  select(datasetid, lat, lon) %>%
  distinct()

scand_pastoral <- scand %>%
  filter(!is.na(pastoral_exceedance_age), pastoral_detected == TRUE) %>%
  left_join(scand_lat, by = "datasetid") %>%
  mutate(lag = pastoral_exceedance_age - neo_scand,
         region = "Scandinavia")

cat(sprintf("Scandinavia: %d sites with pastoral exceedance (out of %d total)\n",
            nrow(scand_pastoral), nrow(scand)))

# --- Central Europe ---
ce <- read_csv("/home/ayu/archeco/shared/signal_decomposition_central_europe.csv")

ce_pastoral <- ce %>%
  filter(!is.na(pastoral_exceed_age), has_pastoral == TRUE) %>%
  mutate(lag = pastoral_exceed_age - neo_ce,
         region = "Central Europe")

cat(sprintf("Central Europe: %d sites with pastoral exceedance (out of %d total)\n",
            nrow(ce_pastoral), nrow(ce)))

# --- Summary function ---
summarise_lag <- function(df, region_name, neo_date) {
  lags <- df$lag
  cat(sprintf("\n--- %s (Neolithic arrival = %d BP) ---\n", region_name, neo_date))
  cat(sprintf("  N sites with pastoral detection: %d\n", length(lags)))
  cat(sprintf("  Median lag: %.0f yr\n", median(lags)))
  cat(sprintf("  IQR: [%.0f, %.0f] yr\n", quantile(lags, 0.25), quantile(lags, 0.75)))
  cat(sprintf("  Range: [%.0f, %.0f] yr\n", min(lags), max(lags)))
  cat(sprintf("  Mean lag: %.0f yr\n", mean(lags)))
  cat(sprintf("  SD: %.0f yr\n", sd(lags)))
  cat("\n")

  # Note: ages are in cal BP, so "before Neolithic" means lag > 0
  # pastoral_exceedance_age > neolithic_arrival means pastoral signal appears
  # at an older date = before Neolithic arrival
  pct_before <- mean(lags > 0) * 100
  pct_within1000 <- mean(lags <= 0 & lags >= -1000) * 100
  pct_after3000 <- mean(lags < -3000) * 100
  pct_after2000 <- mean(lags < -2000) * 100

  cat(sprintf("  %% sites pastoral BEFORE Neolithic: %.1f%%\n", pct_before))
  cat(sprintf("  %% sites within 1000yr AFTER Neolithic: %.1f%%\n", pct_within1000))
  cat(sprintf("  %% sites 2000-3000yr after Neolithic: %.1f%%\n",
              mean(lags < -2000 & lags >= -3000) * 100))
  cat(sprintf("  %% sites >3000yr AFTER Neolithic: %.1f%%\n", pct_after3000))

  return(data.frame(
    region = region_name,
    n = length(lags),
    median_lag = median(lags),
    iqr_lo = quantile(lags, 0.25),
    iqr_hi = quantile(lags, 0.75),
    min_lag = min(lags),
    max_lag = max(lags),
    pct_before_neo = pct_before,
    pct_within_1000 = pct_within1000,
    pct_after_3000 = pct_after3000
  ))
}

s1 <- summarise_lag(uk_pastoral, "UK", neo_uk)
s2 <- summarise_lag(scand_pastoral, "Scandinavia", neo_scand)
s3 <- summarise_lag(ce_pastoral, "Central Europe", neo_ce)

summary_table <- bind_rows(s1, s2, s3)

cat("\n\nSUMMARY TABLE:\n")
print(summary_table)

# =====================================================================
# ANALYSIS 2: Archaeological context
# =====================================================================

cat("\n\n")
cat("ANALYSIS 2: ARCHAEOLOGICAL CONTEXT\n")
cat("-" |> rep(50) |> paste(collapse=""), "\n\n")

cat("Key archaeological benchmarks for landscape-scale pastoralism:\n\n")

cat("CENTRAL EUROPE:\n")
cat("  - LBK (Linearbandkeramik): ~7500-6500 BP — small-scale garden cultivation\n")
cat("  - Forest clearance was localized, rapid regeneration after abandonment\n")
cat("  - Bronze Age (~4200-3000 BP): expansion of pastoral economy, open landscapes\n")
cat("  - Expected lag for regional pollen signal: ~3000-4000 yr\n")
cat("  - Key refs: Behre 1988, Fyfe et al. 2015, Roberts et al. 2018\n\n")

cat("UK:\n")
cat("  - Early Neolithic: ~6000 BP — small clearances, elm decline ~5700 BP\n")
cat("  - Middle-Late Bronze Age: ~3500-3000 BP — first sustained landscape opening\n")
cat("  - Iron Age onwards: ~2800-2000 BP — extensive pastoral landscapes\n")
cat("  - Expected lag for regional pollen signal: ~2000-3000 yr\n")
cat("  - Key refs: Fyfe et al. 2013, Woodbridge et al. 2014, Roberts et al. 2018\n\n")

cat("SCANDINAVIA:\n")
cat("  - Funnel Beaker (TRB): ~5800-4800 BP — initial farming in south\n")
cat("  - Late Bronze Age/Pre-Roman Iron Age: ~3000-2000 BP — pastoral expansion\n")
cat("  - Much of Scandinavia remained forested until Iron Age\n")
cat("  - Expected lag for regional pollen signal: ~2500-3500 yr\n")
cat("  - Key refs: Berglund et al. 1991, Gaillard et al. 2010\n\n")

# =====================================================================
# ANALYSIS 3: Systematic variation in lag
# =====================================================================

cat("\n")
cat("ANALYSIS 3: SYSTEMATIC VARIATION IN LAG\n")
cat("-" |> rep(50) |> paste(collapse=""), "\n\n")

# Combine all data with lat/lon and total BC
all_data <- bind_rows(
  uk_pastoral %>% select(datasetid, sitename, lat, lon, lag, region),
  scand_pastoral %>% select(datasetid, sitename, lat, lon, lag, region),
  ce_pastoral %>% select(datasetid, sitename, lat, lon, lag, region)
)

# Add total BC where available
# UK doesn't have total_bc in the composite analysis, use the corrected decomposition
uk_bc <- read_csv("/home/ayu/archeco/shared/signal_decomposition_uk_corrected.csv",
                  show_col_types = FALSE) %>%
  select(datasetid, total_bc) %>%
  distinct()

scand_bc <- scand %>% select(datasetid, bc_total) %>% rename(total_bc = bc_total)
ce_bc <- ce %>% select(datasetid, total_bc)

all_bc <- bind_rows(uk_bc, scand_bc, ce_bc)

all_data <- all_data %>%
  left_join(all_bc, by = "datasetid")

cat("3a. Lag vs Latitude\n")
if (sum(!is.na(all_data$lat)) > 5) {
  cor_lat <- cor.test(all_data$lat, all_data$lag, use = "complete.obs")
  cat(sprintf("  Pearson r = %.3f (p = %.4f, n = %d)\n",
              cor_lat$estimate, cor_lat$p.value,
              sum(!is.na(all_data$lat) & !is.na(all_data$lag))))

  # Also by region
  for (reg in unique(all_data$region)) {
    sub <- all_data %>% filter(region == reg, !is.na(lat))
    if (nrow(sub) > 3) {
      ct <- cor.test(sub$lat, sub$lag, use = "complete.obs")
      cat(sprintf("  %s: r = %.3f (p = %.4f, n = %d)\n",
                  reg, ct$estimate, ct$p.value, nrow(sub)))
    }
  }
} else {
  cat("  Insufficient lat data for correlation\n")
}

cat("\n3b. Lag vs Total Bray-Curtis magnitude\n")
if (sum(!is.na(all_data$total_bc)) > 5) {
  cor_bc <- cor.test(all_data$total_bc, all_data$lag, use = "complete.obs")
  cat(sprintf("  Pearson r = %.3f (p = %.4f, n = %d)\n",
              cor_bc$estimate, cor_bc$p.value,
              sum(!is.na(all_data$total_bc) & !is.na(all_data$lag))))

  for (reg in unique(all_data$region)) {
    sub <- all_data %>% filter(region == reg, !is.na(total_bc))
    if (nrow(sub) > 3) {
      ct <- cor.test(sub$total_bc, sub$lag, use = "complete.obs")
      cat(sprintf("  %s: r = %.3f (p = %.4f, n = %d)\n",
                  reg, ct$estimate, ct$p.value, nrow(sub)))
    }
  }
} else {
  cat("  Insufficient BC data for correlation\n")
}

cat("\n3c. Lag interpretation:\n")
cat("  Positive lag = pastoral signal detected BEFORE regional Neolithic 'arrival'\n")
cat("    (could mean: pre-Neolithic disturbance, or earlier local Neolithic)\n")
cat("  Negative lag = pastoral signal detected AFTER Neolithic arrival\n")
cat("    (expected pattern; question is how long after)\n")

# =====================================================================
# SYNTHESIS
# =====================================================================

cat("\n\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("SYNTHESIS\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n\n")

# Compute overall lag
all_lags <- all_data$lag
cat(sprintf("Overall: %d sites across 3 regions\n", length(all_lags)))
cat(sprintf("  Overall median lag: %.0f yr (IQR: %.0f to %.0f)\n",
            median(all_lags), quantile(all_lags, 0.25), quantile(all_lags, 0.75)))

pct_before_all <- mean(all_lags > 0) * 100
pct_within1k_all <- mean(all_lags <= 0 & all_lags >= -1000) * 100
pct_1k_3k_all <- mean(all_lags < -1000 & all_lags >= -3000) * 100
pct_3k_all <- mean(all_lags < -3000) * 100

cat(sprintf("  Before Neolithic: %.1f%%\n", pct_before_all))
cat(sprintf("  Within 1000yr after: %.1f%%\n", pct_within1k_all))
cat(sprintf("  1000-3000yr after: %.1f%%\n", pct_1k_3k_all))
cat(sprintf("  >3000yr after: %.1f%%\n", pct_3k_all))

cat("\n")
cat("INTERPRETATION:\n")
cat("The lag between Neolithic arrival and pastoral pollen exceedance should be\n")
cat("evaluated against what archaeology actually tells us about early farming:\n\n")
cat("1. Early Neolithic farming was small-scale, localized, often shifting.\n")
cat("   Pollen catchment areas (~20-50 km radius for large lakes) would not\n")
cat("   detect small clearances affecting <5% of the landscape.\n\n")
cat("2. The transition to landscape-scale pastoralism (detectable in regional\n")
cat("   pollen records) typically occurs in the Bronze-Iron Age across Europe.\n\n")
cat("3. Therefore, a ~3000yr lag is NOT necessarily a methodological failure.\n")
cat("   It may accurately reflect the time between:\n")
cat("   - 'arrival of Neolithic economy' (small-scale, archaeologically visible)\n")
cat("   - 'pastoral transformation' (landscape-scale, palynologically visible)\n\n")
cat("4. The reviewer's concern is valid in framing terms: the paper should NOT\n")
cat("   claim to detect 'onset of human impact' but rather 'onset of landscape-\n")
cat("   scale pastoral/arable transformation.'\n\n")
cat("5. This distinction is actually a STRENGTH: it identifies a genuine\n")
cat("   archaeological transition (intensification, not arrival) that has\n")
cat("   ecological significance.\n")

# Save combined data
write_csv(all_data, "/home/ayu/archeco/shared/phase2_lag_data_combined.csv")
cat("\nCombined data saved to shared/phase2_lag_data_combined.csv\n")
