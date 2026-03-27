#!/usr/bin/env Rscript
# Phase 2: Poaceae Sensitivity Analysis
# Reviewer concern: Poaceae (grasses) in "pastoral/arable" category may be wild grasses,
# contaminating the pastoral signal.
#
# Approach:
#   1. UK taxon-level data: Remove all Poaceae taxa from pastoral/arable, recalculate
#   2. Scandinavia/CE: Check if taxon-level separation is possible
#   3. Indicator taxa test: Verify non-Poaceae indicators still support pastoral-first

library(dplyr)

cat("=== Phase 2: Poaceae Sensitivity Analysis ===\n\n")

# --- Helper ---
read_csv_safe <- function(path) {
  read.csv(path, stringsAsFactors = FALSE) %>% as_tibble()
}

# ============================================================
# PART 1: UK taxon-level Poaceae removal
# ============================================================
cat("--- PART 1: UK Taxon-Level Poaceae Removal ---\n\n")

taxon <- read_csv_safe("/home/ayu/archeco/shared/signal_taxon_contributions.csv")
ce_decomp <- read_csv_safe("/home/ayu/archeco/shared/signal_decomposition_central_europe.csv")

cat("Total taxon records:", nrow(taxon), "\n")
cat("Unique sites:", length(unique(taxon$datasetid)), "\n")
cat("Unique categories:", paste(sort(unique(taxon$category)), collapse = ", "), "\n\n")

# Identify all Poaceae-related taxa
poaceae_taxa <- taxon %>%
  filter(grepl("Poaceae", taxon, ignore.case = TRUE)) %>%
  pull(taxon) %>%
  unique()
cat("Poaceae-related taxa found:", paste(poaceae_taxa, collapse = ", "), "\n\n")

# Which categories do Poaceae taxa fall into?
poaceae_cats <- taxon %>%
  filter(grepl("Poaceae", taxon, ignore.case = TRUE)) %>%
  count(category, name = "n_records") %>%
  arrange(desc(n_records))
cat("Poaceae records by category:\n")
print(as.data.frame(poaceae_cats))
cat("\n")

# Identify Poaceae records in pastoral or arable categories
poaceae_pastoral <- taxon %>%
  filter(grepl("Poaceae", taxon, ignore.case = TRUE),
         category %in% c("pastoral", "arable"))

cat("Poaceae records in pastoral/arable categories:", nrow(poaceae_pastoral), "\n")
if (nrow(poaceae_pastoral) > 0) {
  cat("Sites with Poaceae in pastoral/arable:\n")
  print(as.data.frame(poaceae_pastoral %>%
    select(sitename, datasetid, taxon, bc_contribution, category, group) %>%
    arrange(desc(bc_contribution))))
  cat("\n")
}

# --- Calculate original pastoral_arable BC% per site ---
# pastoral_arable = sum of bc_contribution for category in c("pastoral", "arable")
# total_bc = sum of all bc_contribution for that site
site_orig <- taxon %>%
  group_by(sitename, datasetid, group) %>%
  summarise(
    total_bc = sum(bc_contribution, na.rm = TRUE),
    pastoral_arable_bc = sum(bc_contribution[category %in% c("pastoral", "arable")], na.rm = TRUE),
    poaceae_pastoral_bc = sum(bc_contribution[category %in% c("pastoral", "arable") &
                                                grepl("Poaceae", taxon, ignore.case = TRUE)], na.rm = TRUE),
    non_poaceae_pastoral_bc = pastoral_arable_bc - poaceae_pastoral_bc,
    .groups = "drop"
  ) %>%
  mutate(
    orig_pastoral_pct = (pastoral_arable_bc / total_bc) * 100,
    excl_pastoral_pct = (non_poaceae_pastoral_bc / total_bc) * 100,
    poaceae_fraction = ifelse(pastoral_arable_bc > 0,
                              poaceae_pastoral_bc / pastoral_arable_bc * 100, 0)
  )

cat("--- Per-site Poaceae contribution to pastoral signal ---\n")
sites_with_poaceae <- site_orig %>%
  filter(poaceae_pastoral_bc > 0) %>%
  arrange(desc(poaceae_fraction))
cat("Sites where Poaceae contributes to pastoral/arable:\n")
print(as.data.frame(sites_with_poaceae %>%
  select(sitename, datasetid, group, orig_pastoral_pct, excl_pastoral_pct, poaceae_fraction)))
cat("\n")

# --- Compare H1 vs H3: original and Poaceae-excluded ---
cat("--- H1 vs H3 Comparison ---\n\n")

# Map group to classification
site_orig <- site_orig %>%
  mutate(classification = case_when(
    grepl("A_anthropogenic", group) ~ "H1_anthropogenic",
    grepl("B_unexplained", group) ~ "H3_natural",  # treat B as natural/H3
    TRUE ~ group
  ))

# Alternatively, use the CE decomposition file classification if available
# The CE file has classification column - let's use that for CE sites
# For UK (taxon file), use group mapping

h1_sites <- site_orig %>% filter(grepl("A_anthropogenic", group))
h3_sites <- site_orig %>% filter(!grepl("A_anthropogenic", group))

cat("H1 (anthropogenic) sites: n =", nrow(h1_sites), "\n")
cat("H3 (natural/other) sites: n =", nrow(h3_sites), "\n\n")

cat("ORIGINAL pastoral_arable_pct:\n")
cat("  H1 mean:", round(mean(h1_sites$orig_pastoral_pct, na.rm = TRUE), 2), "%\n")
cat("  H1 median:", round(median(h1_sites$orig_pastoral_pct, na.rm = TRUE), 2), "%\n")
cat("  H3 mean:", round(mean(h3_sites$orig_pastoral_pct, na.rm = TRUE), 2), "%\n")
cat("  H3 median:", round(median(h3_sites$orig_pastoral_pct, na.rm = TRUE), 2), "%\n")

if (nrow(h1_sites) >= 2 & nrow(h3_sites) >= 2) {
  wt_orig <- wilcox.test(h1_sites$orig_pastoral_pct, h3_sites$orig_pastoral_pct,
                         alternative = "greater")
  cat("  Mann-Whitney (H1 > H3): W =", wt_orig$statistic, ", p =",
      format(wt_orig$p.value, digits = 4), "\n")
}
cat("\n")

cat("POACEAE-EXCLUDED pastoral_arable_pct:\n")
cat("  H1 mean:", round(mean(h1_sites$excl_pastoral_pct, na.rm = TRUE), 2), "%\n")
cat("  H1 median:", round(median(h1_sites$excl_pastoral_pct, na.rm = TRUE), 2), "%\n")
cat("  H3 mean:", round(mean(h3_sites$excl_pastoral_pct, na.rm = TRUE), 2), "%\n")
cat("  H3 median:", round(median(h3_sites$excl_pastoral_pct, na.rm = TRUE), 2), "%\n")

if (nrow(h1_sites) >= 2 & nrow(h3_sites) >= 2) {
  wt_excl <- wilcox.test(h1_sites$excl_pastoral_pct, h3_sites$excl_pastoral_pct,
                         alternative = "greater")
  cat("  Mann-Whitney (H1 > H3): W =", wt_excl$statistic, ", p =",
      format(wt_excl$p.value, digits = 4), "\n")
}
cat("\n")

# Summary stats on Poaceae fraction
cat("--- Poaceae as fraction of pastoral signal (all sites) ---\n")
cat("  Mean:", round(mean(site_orig$poaceae_fraction, na.rm = TRUE), 2), "%\n")
cat("  Median:", round(median(site_orig$poaceae_fraction, na.rm = TRUE), 2), "%\n")
cat("  Sites with any Poaceae in pastoral:", sum(site_orig$poaceae_pastoral_bc > 0), "/",
    nrow(site_orig), "\n")
cat("  Max Poaceae fraction:", round(max(site_orig$poaceae_fraction, na.rm = TRUE), 2), "%\n\n")

# ============================================================
# PART 2: Scandinavia & Central Europe — can we separate Poaceae?
# ============================================================
cat("--- PART 2: Scandinavia & Central Europe ---\n\n")

scand <- read_csv_safe("/home/ayu/archeco/shared/signal_decomposition_scandinavia.csv")

cat("Scandinavia decomposition: columns =", paste(names(scand), collapse = ", "), "\n")
cat("  n =", nrow(scand), "\n")
cat("  Column for anthropogenic: bc_anthro_pct (not pastoral_arable_pct)\n")
cat("  No taxon-level breakdown available → cannot separate Poaceae\n\n")

cat("Central Europe decomposition: columns =", paste(names(ce_decomp), collapse = ", "), "\n")
cat("  n =", nrow(ce_decomp), "\n")
cat("  Column for pastoral: pastoral_arable_pct\n")
cat("  No taxon-level breakdown available → cannot separate Poaceae\n\n")

# Report what we know about pastoral composition from literature
cat("LIMITATION: Scandinavia and CE data lack taxon-level decomposition.\n")
cat("  Poaceae contribution to pastoral_arable indices cannot be quantified.\n")
cat("  However, the pastoral composite index in these regions typically includes:\n")
cat("    Plantago lanceolata, Rumex, Artemisia, Poaceae, Cerealia-type\n")
cat("  Indirect evidence from indicator taxa tests (Part 3) can assess\n")
cat("  whether non-Poaceae pastoral taxa independently support findings.\n\n")

# ============================================================
# PART 3: Indicator taxa timing — non-Poaceae pastoral indicators
# ============================================================
cat("--- PART 3: Non-Poaceae Indicator Taxa Timing ---\n\n")

# UK indicator test
uk_ind <- read_csv_safe("/home/ayu/archeco/shared/signal_phase_indicator_test.csv")
cat("UK indicator test: n =", nrow(uk_ind), "records from",
    length(unique(uk_ind$datasetid)), "sites\n")
cat("Indicators:", paste(unique(uk_ind$indicator), collapse = ", "), "\n\n")

# Non-Poaceae pastoral indicators: Plantago_lanceolata, Plantago_all, Rumex
# Arable: Cerealia
# Other disturbance: Artemisia, Urtica
non_poaceae_pastoral <- c("Plantago_lanceolata", "Plantago_all", "Rumex")
arable_ind <- c("Cerealia")
other_disturbance <- c("Artemisia", "Urtica")

# For each indicator, check: does it appear BEFORE signal onset (positive lag)?
# indicator_before = TRUE means indicator present before compositional change
cat("UK: Indicator appearance relative to signal onset\n")
cat("(indicator_before = TRUE means taxon present before the vegetation shift)\n\n")

for (ind in unique(uk_ind$indicator)) {
  sub <- uk_ind %>% filter(indicator == ind, !is.na(first_appearance_age))
  n_total <- nrow(uk_ind %>% filter(indicator == ind))
  n_present <- nrow(sub)
  n_before <- sum(sub$indicator_before, na.rm = TRUE)
  n_after <- sum(sub$indicator_after, na.rm = TRUE)
  n_sync <- sum(sub$synchronous_500yr, na.rm = TRUE)

  category <- ifelse(ind %in% non_poaceae_pastoral, "PASTORAL (non-Poaceae)",
               ifelse(ind %in% arable_ind, "ARABLE",
               "OTHER_DISTURBANCE"))

  cat(sprintf("  %s [%s]: present at %d/%d sites, before: %d, after: %d, synchronous: %d\n",
              ind, category, n_present, n_total, n_before, n_after, n_sync))
}
cat("\n")

# Scandinavia/CE indicator test
scce_ind <- read_csv_safe("/home/ayu/archeco/shared/indicator_test_scand_ce.csv")
cat("Scandinavia/CE indicator test: n =", nrow(scce_ind), "records from",
    length(unique(scce_ind$datasetid)), "sites\n")
cat("Indicators:", paste(unique(scce_ind$indicator), collapse = ", "), "\n")
cat("Regions:", paste(unique(scce_ind$region), collapse = ", "), "\n\n")

for (reg in unique(scce_ind$region)) {
  cat(sprintf("  Region: %s\n", reg))
  sub_reg <- scce_ind %>% filter(region == reg)

  for (ind in unique(sub_reg$indicator)) {
    sub <- sub_reg %>% filter(indicator == ind, !is.na(first_appearance_age))
    n_total <- nrow(sub_reg %>% filter(indicator == ind))
    n_present <- nrow(sub)

    if (n_present > 0) {
      # Positive lag = indicator before signal onset (in cal BP, higher age = earlier)
      n_before <- sum(sub$lag > 0, na.rm = TRUE)
      n_after <- sum(sub$lag < 0, na.rm = TRUE)
      n_sync <- sum(abs(sub$lag) <= 500, na.rm = TRUE)
      median_lag <- median(sub$lag, na.rm = TRUE)
    } else {
      n_before <- 0; n_after <- 0; n_sync <- 0; median_lag <- NA
    }

    category <- ifelse(ind == "Plantago", "PASTORAL (non-Poaceae)",
                 ifelse(ind == "Rumex", "PASTORAL (non-Poaceae)",
                 ifelse(ind == "Cerealia", "ARABLE",
                 ifelse(ind == "Artemisia", "OTHER_DISTURBANCE", ind))))

    cat(sprintf("    %s [%s]: present at %d/%d sites, before: %d, after: %d, sync(±500yr): %d, median lag: %s yr\n",
                ind, category, n_present, n_total, n_before, n_after, n_sync,
                ifelse(is.na(median_lag), "NA", round(median_lag))))
  }
  cat("\n")
}

# ============================================================
# PART 4: Pastoral-first temporal ordering without Poaceae
# ============================================================
cat("--- PART 4: Does pastoral-first ordering hold without Poaceae? ---\n\n")

# From UK indicator data, check if non-Poaceae pastoral indicators
# (Plantago, Rumex) appear before arable indicators (Cerealia)

# For each site, find earliest non-Poaceae pastoral vs earliest arable
uk_timing <- uk_ind %>%
  filter(!is.na(first_appearance_age)) %>%
  mutate(type = case_when(
    indicator %in% non_poaceae_pastoral ~ "pastoral_nonPoaceae",
    indicator %in% arable_ind ~ "arable",
    TRUE ~ "other"
  )) %>%
  filter(type %in% c("pastoral_nonPoaceae", "arable")) %>%
  group_by(datasetid, sitename, type) %>%
  summarise(earliest = max(first_appearance_age, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = type, values_from = earliest)

cat("Sites with both non-Poaceae pastoral AND arable indicators:\n")
both <- uk_timing %>%
  filter(!is.na(pastoral_nonPoaceae) & !is.na(arable))

if (nrow(both) > 0) {
  both <- both %>%
    mutate(pastoral_first = pastoral_nonPoaceae > arable,
           lag = pastoral_nonPoaceae - arable)

  cat("  n =", nrow(both), "\n")
  cat("  Pastoral (non-Poaceae) first:", sum(both$pastoral_first), "/", nrow(both), "\n")
  cat("  Arable first:", sum(!both$pastoral_first), "/", nrow(both), "\n")
  cat("  Median lag (pastoral - arable):", round(median(both$lag)), "years\n")
  cat("  Mean lag:", round(mean(both$lag)), "years\n\n")

  cat("  Site-by-site detail:\n")
  print(as.data.frame(both %>%
    select(sitename, datasetid, pastoral_nonPoaceae, arable, lag, pastoral_first) %>%
    arrange(desc(lag))))
} else {
  cat("  No sites with both types of indicators.\n")
}
cat("\n")

# Same for Scandinavia/CE
cat("--- Scandinavia/CE: Pastoral (non-Poaceae) vs Arable timing ---\n\n")
scce_timing <- scce_ind %>%
  filter(!is.na(first_appearance_age)) %>%
  mutate(type = case_when(
    indicator %in% c("Plantago", "Rumex") ~ "pastoral_nonPoaceae",
    indicator == "Cerealia" ~ "arable",
    TRUE ~ "other"
  )) %>%
  filter(type %in% c("pastoral_nonPoaceae", "arable")) %>%
  group_by(region, datasetid, sitename, type) %>%
  summarise(earliest = max(first_appearance_age, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = type, values_from = earliest)

for (reg in unique(scce_timing$region)) {
  sub <- scce_timing %>%
    filter(region == reg, !is.na(pastoral_nonPoaceae) & !is.na(arable))

  if (nrow(sub) > 0) {
    sub <- sub %>%
      mutate(pastoral_first = pastoral_nonPoaceae > arable,
             lag = pastoral_nonPoaceae - arable)

    cat(sprintf("  %s: n = %d sites with both indicators\n", reg, nrow(sub)))
    cat(sprintf("    Pastoral (non-Poaceae) first: %d/%d (%.0f%%)\n",
                sum(sub$pastoral_first), nrow(sub),
                100 * sum(sub$pastoral_first) / nrow(sub)))
    cat(sprintf("    Median lag: %d years\n", round(median(sub$lag))))
    cat(sprintf("    Mean lag: %d years\n\n", round(mean(sub$lag))))
  } else {
    cat(sprintf("  %s: no sites with both pastoral (non-Poaceae) and arable indicators\n\n", reg))
  }
}

# ============================================================
# SUMMARY
# ============================================================
cat("\n=== SUMMARY ===\n\n")

cat("1. Poaceae contamination in UK pastoral signal:\n")
cat(sprintf("   - %d/%d sites have Poaceae in pastoral/arable category\n",
            sum(site_orig$poaceae_pastoral_bc > 0), nrow(site_orig)))
cat(sprintf("   - Mean Poaceae fraction of pastoral signal: %.1f%%\n",
            mean(site_orig$poaceae_fraction, na.rm = TRUE)))
cat(sprintf("   - Median Poaceae fraction: %.1f%%\n",
            median(site_orig$poaceae_fraction, na.rm = TRUE)))
cat("\n")

cat("2. H1 vs H3 pastoral signal:\n")
cat(sprintf("   - Original: H1 mean = %.2f%% vs H3 mean = %.2f%%\n",
            mean(h1_sites$orig_pastoral_pct, na.rm = TRUE),
            mean(h3_sites$orig_pastoral_pct, na.rm = TRUE)))
cat(sprintf("   - Poaceae-excluded: H1 mean = %.2f%% vs H3 mean = %.2f%%\n",
            mean(h1_sites$excl_pastoral_pct, na.rm = TRUE),
            mean(h3_sites$excl_pastoral_pct, na.rm = TRUE)))
cat("\n")

cat("3. Scandinavia/CE: Taxon-level Poaceae separation NOT possible from existing data.\n")
cat("   Limitation documented.\n\n")

cat("4. Non-Poaceae pastoral indicators (Plantago, Rumex) timing:\n")
cat("   See detailed results above.\n\n")

cat("Script complete.\n")
