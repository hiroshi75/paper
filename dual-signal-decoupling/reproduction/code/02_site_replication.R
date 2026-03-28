#!/usr/bin/env Rscript
# =============================================================================
# Analysis 2: Site-Level Replication Rate
#
# Paper 3 (Dual Signal) — Additional Analysis
#
# For each individual pollen site (not composite):
#   1. Determine if Cerealia-type appears (>0.01%) BEFORE the site's AP
#      starts declining (defined as 2 consecutive bins of AP decrease)
#   2. Classify each site as:
#      - "indicator first": Cerealia rise precedes AP decline
#      - "AP decline first": AP decline precedes Cerealia rise
#      - "simultaneous": both occur in the same 200-yr bin
#      - "no clear pattern": neither event detected or insufficient data
#
# Dependencies: neotoma2, dplyr, tidyr
# Output: shared/analysis/output/site_replication_results.csv
#         shared/analysis/output/site_replication_summary.csv
# =============================================================================

library(neotoma2)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BIN_WIDTH <- 200
MIN_SAMPLES <- 3        # minimum samples per site to include
CEREALIA_THRESHOLD <- 0.01  # percentage threshold for "Cerealia present"
OUTPUT_DIR <- "output"

# Region definitions (same as Analysis 1)
regions <- list(
  britain = list(
    name = "Britain",
    bbox = c(-6, 50, 2, 57),
    ag_start = 6000
  ),
  scandinavia = list(
    name = "Scandinavia",
    bbox = c(5, 55, 30, 70),
    ag_start = 5500
  ),
  alps = list(
    name = "Alps",
    bbox = c(6, 45, 16, 48),
    ag_start = 7500
  )
)

# Neotoma ecological groups
AP_GROUPS <- c("TRSH", "MANG", "PALM")
NAP_GROUPS <- c("UPHE", "SUCC")

# Agricultural indicator taxa
CEREALIA_TAXA <- c("Cerealia", "Cerealia-type", "Cerealia undiff.",
                   "Triticum", "Hordeum", "Secale", "Avena",
                   "Triticum-type", "Hordeum-type", "Secale-type")

bin_edges <- seq(0, 10000, by = BIN_WIDTH)
bin_mids  <- bin_edges[-length(bin_edges)] + BIN_WIDTH / 2

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

assign_bin <- function(age) {
  idx <- findInterval(age, bin_edges, rightmost.closed = TRUE)
  idx <- pmin(pmax(idx, 1), length(bin_mids))
  bin_mids[idx]
}


#' Fetch and process pollen data for a region.
#' Returns a data frame with per-sample records including
#' site_id, site_name, age, bin, ap_ratio, cerealia_pct
fetch_region_data <- function(bbox) {
  sites <- get_sites(loc = bbox)
  if (length(sites) == 0) return(data.frame())

  datasets <- get_datasets(sites, datasettype = "pollen")
  if (length(datasets) == 0) return(data.frame())

  dl <- get_downloads(datasets)

  all_results <- list()

  for (i in seq_along(dl)) {
    tryCatch({
      samples_df <- samples(dl[[i]])
      if (nrow(samples_df) == 0) next

      sample_results <- samples_df %>%
        mutate(
          is_ap  = ecologicalgroup %in% AP_GROUPS,
          is_nap = ecologicalgroup %in% NAP_GROUPS,
          is_cerealia = variablename %in% CEREALIA_TAXA
        ) %>%
        group_by(sampleid, age) %>%
        summarise(
          ap_sum      = sum(value[is_ap], na.rm = TRUE),
          nap_sum     = sum(value[is_nap], na.rm = TRUE),
          total_sum   = ap_sum + nap_sum,
          cerealia_ct = sum(value[is_cerealia], na.rm = TRUE),
          total_pollen = sum(value, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        filter(total_sum > 10, !is.na(age)) %>%
        mutate(
          ap_ratio     = ap_sum / total_sum,
          cerealia_pct = (cerealia_ct / total_pollen) * 100,
          bin          = assign_bin(age)
        )

      if (nrow(sample_results) < MIN_SAMPLES) next

      sid <- unique(dl[[i]]@sites@sites$siteid)[1]
      sname <- unique(dl[[i]]@sites@sites$sitename)[1]
      sample_results$site_id <- sid
      sample_results$site_name <- sname
      all_results[[length(all_results) + 1]] <- sample_results

    }, error = function(e) {
      message("Skipping site ", i, ": ", e$message)
    })
  }

  if (length(all_results) == 0) return(data.frame())
  bind_rows(all_results)
}


#' For a single site, compute binned AP and Cerealia timeseries
#' and classify the temporal ordering of events.
#'
#' Returns a list with:
#'   classification: one of "indicator first", "AP decline first",
#'                   "simultaneous", "no clear pattern"
#'   cerealia_onset_bin: bin (cal BP) of first Cerealia appearance
#'   ap_decline_bin:     bin (cal BP) of first sustained AP decline
classify_site <- function(site_df) {
  # Compute bin-level averages for this site
  site_bins <- site_df %>%
    group_by(bin) %>%
    summarise(
      ap_ratio     = mean(ap_ratio, na.rm = TRUE),
      cerealia_pct = mean(cerealia_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(bin))   # oldest first (high cal BP)

  if (nrow(site_bins) < 4) {
    return(list(
      classification = "no clear pattern",
      cerealia_onset = NA_real_,
      ap_decline     = NA_real_,
      reason         = "insufficient bins"
    ))
  }

  # --- Find Cerealia onset ---
  # First bin (going forward in time, i.e., decreasing cal BP)
  # where Cerealia > threshold
  cerealia_onset <- NA_real_
  for (j in seq_len(nrow(site_bins))) {
    if (site_bins$cerealia_pct[j] > CEREALIA_THRESHOLD) {
      cerealia_onset <- site_bins$bin[j]
      break
    }
  }

  # --- Find AP decline onset ---
  # First bin where AP decreases for 2 consecutive bins
  # (i.e., 3 consecutive bins with monotonic decline going forward in time)
  # We walk forward in time (rows are oldest-first)
  ap_decline <- NA_real_
  for (j in 1:(nrow(site_bins) - 2)) {
    # Bins are ordered oldest to newest (decreasing cal BP)
    ap1 <- site_bins$ap_ratio[j]
    ap2 <- site_bins$ap_ratio[j + 1]
    ap3 <- site_bins$ap_ratio[j + 2]
    if (ap2 < ap1 && ap3 < ap2) {
      # AP declined over 2 consecutive steps starting at bin j
      ap_decline <- site_bins$bin[j + 1]  # decline starts at 2nd bin
      break
    }
  }

  # --- Classify ---
  if (is.na(cerealia_onset) && is.na(ap_decline)) {
    classification <- "no clear pattern"
    reason <- "neither event detected"
  } else if (is.na(cerealia_onset)) {
    # AP declined but no Cerealia detected
    classification <- "AP decline first"
    reason <- "AP decline without Cerealia"
  } else if (is.na(ap_decline)) {
    # Cerealia appeared but no sustained AP decline
    classification <- "indicator first"
    reason <- "Cerealia without AP decline"
  } else if (cerealia_onset == ap_decline) {
    classification <- "simultaneous"
    reason <- "same bin"
  } else if (cerealia_onset > ap_decline) {
    # cerealia_onset is higher cal BP = older = earlier in time
    classification <- "indicator first"
    reason <- paste0("Cerealia at ", cerealia_onset, " BP, AP decline at ",
                     ap_decline, " BP")
  } else {
    # AP decline is older
    classification <- "AP decline first"
    reason <- paste0("AP decline at ", ap_decline, " BP, Cerealia at ",
                     cerealia_onset, " BP")
  }

  return(list(
    classification = classification,
    cerealia_onset = cerealia_onset,
    ap_decline     = ap_decline,
    reason         = reason
  ))
}

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

site_results_all <- data.frame()
summary_all <- data.frame()

for (reg_key in names(regions)) {
  reg <- regions[[reg_key]]
  cat("\n=== Processing:", reg$name, "===\n")

  # Fetch data
  cat("Fetching pollen data...\n")
  region_df <- fetch_region_data(reg$bbox)

  if (nrow(region_df) == 0) {
    cat("  No data retrieved. Skipping.\n")
    next
  }

  site_ids <- unique(region_df$site_id)
  cat("  Sites with data:", length(site_ids), "\n")

  # Classify each site
  site_classifications <- list()

  for (sid in site_ids) {
    site_df <- region_df %>% filter(site_id == sid)
    sname <- unique(site_df$site_name)[1]
    result <- classify_site(site_df)

    site_classifications[[length(site_classifications) + 1]] <- data.frame(
      region          = reg$name,
      site_id         = sid,
      site_name       = sname,
      n_bins          = n_distinct(site_df$bin),
      classification  = result$classification,
      cerealia_onset  = ifelse(is.na(result$cerealia_onset), NA,
                               result$cerealia_onset),
      ap_decline_bin  = ifelse(is.na(result$ap_decline), NA,
                               result$ap_decline),
      reason          = result$reason,
      stringsAsFactors = FALSE
    )
  }

  site_results <- bind_rows(site_classifications)
  site_results_all <- bind_rows(site_results_all, site_results)

  # Summarise
  n_total <- nrow(site_results)
  counts <- table(site_results$classification)

  summary_row <- data.frame(
    region = reg$name,
    total_sites = n_total,
    n_indicator_first  = as.integer(counts["indicator first"] %||% 0),
    n_AP_decline_first = as.integer(counts["AP decline first"] %||% 0),
    n_simultaneous     = as.integer(counts["simultaneous"] %||% 0),
    n_no_clear_pattern = as.integer(counts["no clear pattern"] %||% 0),
    stringsAsFactors = FALSE
  )

  # Replace NAs from missing factor levels with 0
  summary_row[is.na(summary_row)] <- 0

  summary_row <- summary_row %>%
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

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

cat("\n\n========================================\n")
cat("SITE-LEVEL REPLICATION SUMMARY\n")
cat("========================================\n\n")
print(summary_all, row.names = FALSE)

# Save detailed site-level results
site_file <- file.path(OUTPUT_DIR, "site_replication_results.csv")
write.csv(site_results_all, site_file, row.names = FALSE)
cat("\nSite-level results saved to:", site_file, "\n")

# Save summary table
summary_file <- file.path(OUTPUT_DIR, "site_replication_summary.csv")
write.csv(summary_all, summary_file, row.names = FALSE)
cat("Summary saved to:", summary_file, "\n")

cat("\nAnalysis complete.\n")
