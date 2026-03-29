#!/usr/bin/env Rscript
# =============================================================================
# Pilot Exceedance Analysis: China Pollen Data — Third Continental Test
# =============================================================================
# Tests the identifiability framework (Paper 8) on Chinese pollen data.
# Prediction: European indicators fail, Chinese indicators detect rice agriculture,
# tree exceedance ~93% (universal).
#
# Three indicator sets:
#   Set A (European): Plantago, Rumex, Cerealia-type, Secale
#   Set B (Chinese):  Oryza, Cannabis, Fagopyrum (rice-based agriculture)
#   Set C (Tree):     Quercus, Pinus, Betula, Fagus, Castanopsis,
#                     Cyclobalanopsis, Liquidambar
#
# Author: ArchEco research agent
# Date: 2026-03-28
# =============================================================================

suppressPackageStartupMessages({
  library(neotoma2)
  library(vegan)
})

# ============================================================
# Configuration
# ============================================================
CACHE_FILE <- "/home/ayu/archeco/shared/cache/china_neotoma2_downloads.rds"
OUT_DIR    <- "/home/ayu/archeco/shared"
SCRIPT_DIR <- "/home/ayu/archeco/shared/scripts"
BASELINE_CUTOFF <- 5000  # BP: >5000 = pre-rice-intensification baseline
EXCEEDANCE_SD   <- 2     # threshold: mean + N*SD
MIN_HOLOCENE_SAMPLES <- 5
MIN_BASELINE_SAMPLES <- 1  # at least 1 sample >5000 BP

# Indicator sets (case-insensitive regex matching)
INDICATOR_SETS <- list(
  european = list(
    name = "Set A: European (Behre 1981)",
    patterns = c("^Plantago", "^Rumex", "^Cerealia", "^Secale",
                 "Plantago lanceolata", "Plantago major")
  ),
  chinese = list(
    name = "Set B: Chinese agriculture",
    patterns = c("^Oryza", "Oryza-type", "^Cannabis", "^Fagopyrum",
                 "^Buckwheat", "^Humulus", "Cannabis/Humulus")
  ),
  tree = list(
    name = "Set C: Tree taxa",
    patterns = c("^Quercus", "^Pinus", "^Betula", "^Fagus",
                 "^Castanopsis", "^Cyclobalanopsis", "^Liquidambar",
                 "^Cryptomeria", "^Cunninghamia", "^Podocarpus",
                 "^Tsuga", "^Picea", "^Abies", "^Ulmus", "^Zelkova")
  )
)

dir.create(dirname(CACHE_FILE), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Step 1: Download Chinese pollen data from Neotoma
# ============================================================
cat("=== Step 1: Getting pollen data from Neotoma (China) ===\n")

if (file.exists(CACHE_FILE)) {
  cat("  Loading from cache...\n")
  all_samples <- readRDS(CACHE_FILE)
  cat(sprintf("  Cached: %d rows from %d sites\n",
              nrow(all_samples), length(unique(all_samples$siteid))))
} else {
  cat("  Querying Neotoma for China pollen sites...\n")
  cat("  Bounding box: lon 75-135, lat 20-55\n")

  # China bounding box
  sites_obj <- get_sites(loc = c(75, 20, 135, 55),
                         datasettype = "pollen",
                         limit = 200)
  cat(sprintf("  Found %d sites\n", length(sites_obj)))

  if (length(sites_obj) == 0) {
    stop("No sites found. Check Neotoma connectivity.")
  }

  # Download all
  cat(sprintf("  Downloading %d sites...\n", length(sites_obj)))
  ds_obj <- get_downloads(sites_obj)

  # Extract samples
  cat("  Extracting sample data...\n")
  all_samples <- samples(ds_obj)

  # Keep only pollen datasets
  all_samples <- all_samples[all_samples$datasettype == "pollen", ]

  cat(sprintf("  Total pollen samples: %d rows, %d sites\n",
              nrow(all_samples), length(unique(all_samples$siteid))))

  # Cache
  saveRDS(all_samples, CACHE_FILE)
  cat("  Cached to disk.\n")
}

# ============================================================
# Step 2: Quality filter
# ============================================================
cat("\n=== Step 2: Filtering sites ===\n")

site_summary <- do.call(rbind, lapply(split(all_samples, all_samples$siteid), function(df) {
  ages <- unique(df$age)
  ages <- ages[!is.na(ages)]
  data.frame(
    siteid   = df$siteid[1],
    sitename = df$sitename[1],
    lat      = df$lat[1],
    lon      = df$long[1],
    n_ages   = length(ages),
    age_min  = if (length(ages) > 0) min(ages) else NA,
    age_max  = if (length(ages) > 0) max(ages) else NA,
    has_baseline = any(ages > BASELINE_CUTOFF),
    n_baseline   = sum(ages > BASELINE_CUTOFF),
    n_holocene   = sum(ages >= 0 & ages <= 11700),
    stringsAsFactors = FALSE
  )
}))

# Filter: at least 1 baseline sample and >=5 Holocene samples
good <- site_summary$has_baseline &
        site_summary$n_baseline >= MIN_BASELINE_SAMPLES &
        site_summary$n_holocene >= MIN_HOLOCENE_SAMPLES &
        site_summary$n_ages >= 5

good_sites <- site_summary[good, ]

cat(sprintf("  Total sites: %d\n", nrow(site_summary)))
cat(sprintf("  Sites with >=1 sample >5000 BP: %d\n", sum(site_summary$has_baseline)))
cat(sprintf("  Sites with >=5 Holocene samples: %d\n", sum(site_summary$n_holocene >= MIN_HOLOCENE_SAMPLES)))
cat(sprintf("  Sites passing both filters: %d\n", nrow(good_sites)))

if (nrow(good_sites) == 0) {
  cat("\n  WARNING: No sites pass filters. Relaxing to n_holocene >= 3...\n")
  good <- site_summary$has_baseline &
          site_summary$n_baseline >= MIN_BASELINE_SAMPLES &
          site_summary$n_holocene >= 3 &
          site_summary$n_ages >= 4
  good_sites <- site_summary[good, ]
  cat(sprintf("  Sites passing relaxed filters: %d\n", nrow(good_sites)))
}

cat("\n  Passing sites:\n")
for (i in 1:min(nrow(good_sites), 30)) {
  s <- good_sites[i, ]
  cat(sprintf("    %s (id=%d): lat=%.2f, lon=%.2f, n_ages=%d, baseline=%d, holocene=%d\n",
              s$sitename, s$siteid, s$lat, s$lon, s$n_ages, s$n_baseline, s$n_holocene))
}
if (nrow(good_sites) > 30) cat(sprintf("    ... and %d more\n", nrow(good_sites) - 30))

# ============================================================
# Step 3: Taxon inventory — check what indicators are present
# ============================================================
cat("\n=== Step 3: Taxon inventory ===\n")

good_ids <- good_sites$siteid
good_samples <- all_samples[all_samples$siteid %in% good_ids, ]
all_taxa <- sort(unique(good_samples$variablename))

cat(sprintf("  Total unique taxa across %d sites: %d\n", length(good_ids), length(all_taxa)))

# Check each indicator set
for (set_name in names(INDICATOR_SETS)) {
  iset <- INDICATOR_SETS[[set_name]]
  cat(sprintf("\n  %s:\n", iset$name))
  found_taxa <- character()
  for (pat in iset$patterns) {
    matches <- all_taxa[grepl(pat, all_taxa, ignore.case = TRUE)]
    if (length(matches) > 0) {
      found_taxa <- c(found_taxa, matches)
      for (m in matches) {
        n_sites_with <- length(unique(good_samples$siteid[
          grepl(paste0("^", gsub("([.+?^${}()|\\[\\]])", "\\\\\\1", m), "$"),
                good_samples$variablename, ignore.case = TRUE) &
          good_samples$value > 0]))
        cat(sprintf("    [FOUND] %-30s  (%d sites)\n", m, n_sites_with))
      }
    }
  }
  if (length(found_taxa) == 0) {
    cat("    [NONE FOUND]\n")
  }
}

# ============================================================
# Step 4: Exceedance analysis per site per indicator set
# ============================================================
cat("\n=== Step 4: Exceedance analysis ===\n")

# Function: check if any taxon in an indicator set exceeds baseline threshold
exceedance_for_set <- function(sid, set_patterns, mat_prop, ages, taxa) {
  # Find columns matching patterns
  matched_cols <- integer()
  for (pat in set_patterns) {
    idx <- grep(pat, taxa, ignore.case = TRUE)
    matched_cols <- c(matched_cols, idx)
  }
  matched_cols <- unique(matched_cols)

  if (length(matched_cols) == 0) {
    return(list(testable = FALSE, exceedance = FALSE, onset_bp = NA,
                n_taxa_matched = 0, taxa_matched = ""))
  }

  # Get the subset
  sub_mat <- mat_prop[, matched_cols, drop = FALSE]
  taxa_names <- taxa[matched_cols]

  # Baseline and test periods
  baseline_idx <- which(ages > BASELINE_CUTOFF)
  test_idx     <- which(ages <= BASELINE_CUTOFF & ages >= 0)

  if (length(baseline_idx) < 1 || length(test_idx) < 1) {
    return(list(testable = FALSE, exceedance = FALSE, onset_bp = NA,
                n_taxa_matched = length(matched_cols),
                taxa_matched = paste(taxa_names, collapse = "; ")))
  }

  # Per-taxon: baseline mean + 2SD threshold
  any_exceedance <- FALSE
  earliest_onset <- -Inf

  for (j in seq_len(ncol(sub_mat))) {
    baseline_vals <- sub_mat[baseline_idx, j]
    test_vals     <- sub_mat[test_idx, j]
    test_ages_j   <- ages[test_idx]

    b_mean <- mean(baseline_vals, na.rm = TRUE)
    b_sd   <- sd(baseline_vals, na.rm = TRUE)
    if (is.na(b_sd)) b_sd <- 0

    threshold <- b_mean + EXCEEDANCE_SD * b_sd

    # If baseline is all zeros, threshold = 0; any positive value exceeds
    if (b_mean == 0 && b_sd == 0) {
      exc <- test_vals > 0
    } else {
      exc <- test_vals > threshold
    }

    if (any(exc, na.rm = TRUE)) {
      any_exceedance <- TRUE
      onset <- max(test_ages_j[exc], na.rm = TRUE)
      if (onset > earliest_onset) earliest_onset <- onset
    }
  }

  list(
    testable = TRUE,
    exceedance = any_exceedance,
    onset_bp = if (any_exceedance) earliest_onset else NA,
    n_taxa_matched = length(matched_cols),
    taxa_matched = paste(taxa_names, collapse = "; ")
  )
}

# Main analysis loop
results_list <- list()
n_fail <- 0

for (i in seq_along(good_ids)) {
  sid <- good_ids[i]
  tryCatch({
    df <- all_samples[all_samples$siteid == sid, ]
    df <- df[!is.na(df$age) & !is.na(df$value), ]

    if (nrow(df) < 10) { n_fail <- n_fail + 1; next }

    meta <- good_sites[good_sites$siteid == sid, ]

    # Build sample-taxon matrix
    ages <- sort(unique(df$age))
    ages <- ages[!is.na(ages)]
    taxa <- unique(df$variablename)
    taxa <- taxa[!is.na(taxa)]

    if (length(ages) < 4 || length(taxa) < 3) { n_fail <- n_fail + 1; next }

    mat <- matrix(0, nrow = length(ages), ncol = length(taxa),
                  dimnames = list(seq_along(ages), taxa))

    ai <- match(df$age, ages)
    ti <- match(df$variablename, taxa)
    valid <- !is.na(ai) & !is.na(ti)
    for (r in which(valid)) {
      mat[ai[r], ti[r]] <- mat[ai[r], ti[r]] + df$value[r]
    }

    # Remove non-pollen columns
    bad_cols <- grep("^(Sample|Concentration|Influx|Lycopodium|Spike|Microsphere|Charcoal|Carbon)",
                     taxa, ignore.case = TRUE)
    if (length(bad_cols) > 0) {
      mat <- mat[, -bad_cols, drop = FALSE]
      taxa <- taxa[-bad_cols]
    }

    # Remove very rare taxa
    taxa_presence <- colSums(mat > 0)
    keep <- taxa_presence >= 2
    mat <- mat[, keep, drop = FALSE]
    taxa <- taxa[keep]

    if (ncol(mat) < 3) { n_fail <- n_fail + 1; next }

    # Convert to proportions
    row_sums <- rowSums(mat)
    mat_prop <- mat / ifelse(row_sums == 0, 1, row_sums)

    # Run exceedance for each indicator set
    res_eur <- exceedance_for_set(sid, INDICATOR_SETS$european$patterns, mat_prop, ages, taxa)
    res_chi <- exceedance_for_set(sid, INDICATOR_SETS$chinese$patterns, mat_prop, ages, taxa)
    res_tree <- exceedance_for_set(sid, INDICATOR_SETS$tree$patterns, mat_prop, ages, taxa)

    # BC dissimilarity exceedance (whole assemblage, as in ENA pilot)
    baseline_idx <- which(ages > BASELINE_CUTOFF)
    test_idx     <- which(ages <= BASELINE_CUTOFF & ages >= 0)

    bc_signal <- FALSE
    bc_onset  <- NA
    max_bc    <- NA

    if (length(baseline_idx) >= 1 && length(test_idx) >= 1) {
      baseline_centroid <- colMeans(mat_prop[baseline_idx, , drop = FALSE])
      combined <- rbind(mat_prop, centroid = baseline_centroid)
      bc_all <- as.matrix(vegdist(combined, method = "bray"))
      centroid_row <- nrow(combined)
      bc_from_centroid <- bc_all[1:(centroid_row - 1), centroid_row]

      bc_baseline <- bc_from_centroid[baseline_idx]
      bc_mean <- mean(bc_baseline, na.rm = TRUE)
      bc_sd   <- sd(bc_baseline, na.rm = TRUE)
      if (!is.na(bc_sd) && bc_sd > 0) {
        threshold <- bc_mean + EXCEEDANCE_SD * bc_sd
        bc_test <- bc_from_centroid[test_idx]
        test_ages <- ages[test_idx]
        exceedances <- bc_test > threshold
        bc_signal <- any(exceedances, na.rm = TRUE)
        bc_onset <- if (bc_signal) max(test_ages[exceedances], na.rm = TRUE) else NA
        max_bc <- max(bc_test, na.rm = TRUE)
      }
    }

    results_list[[length(results_list) + 1]] <- data.frame(
      siteid = as.integer(sid),
      sitename = as.character(meta$sitename[1]),
      lat = as.numeric(meta$lat[1]),
      lon = as.numeric(meta$lon[1]),
      n_ages = length(ages),
      n_baseline = length(baseline_idx),
      n_test = length(test_idx),
      n_taxa = ncol(mat),
      # European indicators
      eur_testable = res_eur$testable,
      eur_exceedance = res_eur$exceedance,
      eur_onset_bp = res_eur$onset_bp,
      eur_n_taxa = res_eur$n_taxa_matched,
      eur_taxa = res_eur$taxa_matched,
      # Chinese indicators
      chi_testable = res_chi$testable,
      chi_exceedance = res_chi$exceedance,
      chi_onset_bp = res_chi$onset_bp,
      chi_n_taxa = res_chi$n_taxa_matched,
      chi_taxa = res_chi$taxa_matched,
      # Tree indicators
      tree_testable = res_tree$testable,
      tree_exceedance = res_tree$exceedance,
      tree_onset_bp = res_tree$onset_bp,
      tree_n_taxa = res_tree$n_taxa_matched,
      tree_taxa = res_tree$taxa_matched,
      # Whole-assemblage BC
      bc_signal = bc_signal,
      bc_onset_bp = bc_onset,
      max_bc = if (!is.na(max_bc)) round(max_bc, 4) else NA,
      row.names = NULL,
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    cat(sprintf("  ERROR site %d: %s\n", sid, e$message))
    n_fail <<- n_fail + 1
  })
}

cat(sprintf("  Succeeded: %d, Failed: %d\n", length(results_list), n_fail))

if (length(results_list) == 0) {
  cat("\n  FATAL: No sites produced results.\n")
  quit(status = 1)
}

results <- do.call(rbind, results_list)
rownames(results) <- NULL

# ============================================================
# Step 5: Summary statistics
# ============================================================
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("PILOT CHINA EXCEEDANCE ANALYSIS RESULTS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

n_total <- nrow(results)

# European indicators
n_eur_testable <- sum(results$eur_testable)
n_eur_exc      <- sum(results$eur_exceedance)
eur_rate       <- if (n_eur_testable > 0) n_eur_exc / n_eur_testable * 100 else 0
eur_rate_total <- n_eur_exc / n_total * 100

# Chinese indicators
n_chi_testable <- sum(results$chi_testable)
n_chi_exc      <- sum(results$chi_exceedance)
chi_rate       <- if (n_chi_testable > 0) n_chi_exc / n_chi_testable * 100 else 0
chi_rate_total <- n_chi_exc / n_total * 100

# Tree indicators
n_tree_testable <- sum(results$tree_testable)
n_tree_exc      <- sum(results$tree_exceedance)
tree_rate       <- if (n_tree_testable > 0) n_tree_exc / n_tree_testable * 100 else 0
tree_rate_total <- n_tree_exc / n_total * 100

# BC whole-assemblage
n_bc_signal <- sum(results$bc_signal)
bc_rate     <- n_bc_signal / n_total * 100

cat(sprintf("1. INDICATOR SET EXCEEDANCE RATES\n\n"))
cat(sprintf("   %-35s  Testable  Exceedance  Rate(testable)  Rate(total)\n", "Indicator Set"))
cat(sprintf("   %-35s  --------  ----------  --------------  -----------\n", ""))
cat(sprintf("   %-35s  %4d      %4d        %5.1f%%          %5.1f%%\n",
            "Set A: European (Behre 1981)", n_eur_testable, n_eur_exc, eur_rate, eur_rate_total))
cat(sprintf("   %-35s  %4d      %4d        %5.1f%%          %5.1f%%\n",
            "Set B: Chinese agriculture", n_chi_testable, n_chi_exc, chi_rate, chi_rate_total))
cat(sprintf("   %-35s  %4d      %4d        %5.1f%%          %5.1f%%\n",
            "Set C: Tree taxa", n_tree_testable, n_tree_exc, tree_rate, tree_rate_total))
cat(sprintf("   %-35s  %4d      %4d        %5.1f%%          %5.1f%%\n",
            "Whole-assemblage BC", n_total, n_bc_signal, bc_rate, bc_rate))

cat(sprintf("\n2. ONSET TIMING (among exceedance sites)\n\n"))

for (set_label in c("European", "Chinese", "Tree")) {
  col <- switch(set_label,
                European = "eur_onset_bp",
                Chinese  = "chi_onset_bp",
                Tree     = "tree_onset_bp")
  onsets <- results[[col]]
  onsets <- onsets[!is.na(onsets)]
  if (length(onsets) > 0) {
    cat(sprintf("   %s: median=%.0f BP, mean=%.0f BP, range=%d-%d BP (n=%d)\n",
                set_label, median(onsets), mean(onsets),
                min(onsets), max(onsets), length(onsets)))
  } else {
    cat(sprintf("   %s: no exceedances\n", set_label))
  }
}

cat(sprintf("\n3. INDICATOR DEPENDENCY GAP\n"))
rates <- c(eur_rate_total, chi_rate_total, tree_rate_total)
gap <- max(rates) - min(rates)
cat(sprintf("   Max detection: %.1f%% (%s)\n",
            max(rates), c("European", "Chinese", "Tree")[which.max(rates)]))
cat(sprintf("   Min detection: %.1f%% (%s)\n",
            min(rates), c("European", "Chinese", "Tree")[which.min(rates)]))
cat(sprintf("   Gap: %.1f percentage points\n", gap))

cat(sprintf("\n4. COMPARISON WITH PREDICTIONS\n\n"))
cat(sprintf("   | Indicator     | Europe | ENA     | China predicted | China observed |\n"))
cat(sprintf("   |---------------|--------|---------|-----------------|----------------|\n"))
cat(sprintf("   | European      | 87%%    | 16.1%%  | ~0%%             | %.1f%%          |\n", eur_rate_total))
cat(sprintf("   | Region-spec.  | --     | 26.8%%  | >0%%?            | %.1f%%          |\n", chi_rate_total))
cat(sprintf("   | Tree          | 93.3%% | 26.8%%  | ~93%%?           | %.1f%%          |\n", tree_rate_total))

# ============================================================
# Step 6: Save outputs
# ============================================================
cat("\n=== Saving outputs ===\n")

# CSV
csv_path <- file.path(OUT_DIR, "pilot_china_exceedance_data.csv")
write.csv(results, csv_path, row.names = FALSE)
cat(sprintf("  Saved CSV: %s\n", csv_path))

# Markdown report
report_lines <- c(
  "# Pilot China Exceedance Analysis Results",
  "",
  sprintf("**Date**: %s", Sys.Date()),
  sprintf("**Sites analyzed**: %d (of %d downloaded, %d in bounding box)",
          n_total, nrow(site_summary), nrow(site_summary)),
  "**Framework**: Paper 8 indicator-dependent exceedance",
  sprintf("**Baseline**: >%d BP (pre-rice-intensification)", BASELINE_CUTOFF),
  sprintf("**Threshold**: mean + %dSD of baseline per taxon", EXCEEDANCE_SD),
  "",
  "## Summary: Indicator Set Exceedance Rates",
  "",
  "| Indicator Set | N Testable | N Exceedance | Rate (testable) | Rate (total) |",
  "|---|---|---|---|---|",
  sprintf("| Set A: European (Behre 1981) | %d | %d | %.1f%% | %.1f%% |",
          n_eur_testable, n_eur_exc, eur_rate, eur_rate_total),
  sprintf("| Set B: Chinese agriculture | %d | %d | %.1f%% | %.1f%% |",
          n_chi_testable, n_chi_exc, chi_rate, chi_rate_total),
  sprintf("| Set C: Tree taxa | %d | %d | %.1f%% | %.1f%% |",
          n_tree_testable, n_tree_exc, tree_rate, tree_rate_total),
  sprintf("| Whole-assemblage BC | %d | %d | %.1f%% | %.1f%% |",
          n_total, n_bc_signal, bc_rate, bc_rate),
  "",
  "## Cross-Continental Comparison",
  "",
  "| Indicator | Europe | ENA | China (predicted) | China (observed) |",
  "|---|---|---|---|---|",
  sprintf("| European | 87%% | 16.1%% | ~0%% | **%.1f%%** |", eur_rate_total),
  sprintf("| Region-specific | -- | 26.8%% (EAC) | >0%% | **%.1f%%** |", chi_rate_total),
  sprintf("| Tree | 93.3%% | 26.8%% | ~93%% | **%.1f%%** |", tree_rate_total),
  "",
  sprintf("## Indicator Dependency Gap: %.1f pp", gap),
  "",
  sprintf("- Max detection: %.1f%% (%s)", max(rates), c("European", "Chinese", "Tree")[which.max(rates)]),
  sprintf("- Min detection: %.1f%% (%s)", min(rates), c("European", "Chinese", "Tree")[which.min(rates)]),
  "",
  "## Onset Timing",
  ""
)

for (set_label in c("European", "Chinese", "Tree")) {
  col <- switch(set_label, European = "eur_onset_bp", Chinese = "chi_onset_bp", Tree = "tree_onset_bp")
  onsets <- results[[col]]
  onsets <- onsets[!is.na(onsets)]
  if (length(onsets) > 0) {
    report_lines <- c(report_lines,
      sprintf("### %s indicators (n=%d exceedances)", set_label, length(onsets)),
      sprintf("- Median onset: %.0f BP", median(onsets)),
      sprintf("- Mean onset: %.0f BP", mean(onsets)),
      sprintf("- Range: %d - %d BP", min(onsets), max(onsets)),
      "")
  } else {
    report_lines <- c(report_lines,
      sprintf("### %s indicators: no exceedances", set_label), "")
  }
}

# Taxon inventory section
report_lines <- c(report_lines,
  "## Taxon Inventory",
  "",
  "### European indicator taxa found",
  ""
)
eur_found <- unique(unlist(lapply(results$eur_taxa[results$eur_taxa != ""], function(x) strsplit(x, "; ")[[1]])))
if (length(eur_found) > 0) {
  for (t in eur_found) report_lines <- c(report_lines, sprintf("- %s", t))
} else {
  report_lines <- c(report_lines, "- None found")
}

report_lines <- c(report_lines, "", "### Chinese indicator taxa found", "")
chi_found <- unique(unlist(lapply(results$chi_taxa[results$chi_taxa != ""], function(x) strsplit(x, "; ")[[1]])))
if (length(chi_found) > 0) {
  for (t in chi_found) report_lines <- c(report_lines, sprintf("- %s", t))
} else {
  report_lines <- c(report_lines, "- None found")
}

report_lines <- c(report_lines, "", "### Tree indicator taxa found", "")
tree_found <- unique(unlist(lapply(results$tree_taxa[results$tree_taxa != ""], function(x) strsplit(x, "; ")[[1]])))
if (length(tree_found) > 0) {
  for (t in sort(tree_found)) report_lines <- c(report_lines, sprintf("- %s", t))
} else {
  report_lines <- c(report_lines, "- None found")
}

# Site-level table
report_lines <- c(report_lines,
  "",
  "## Site-Level Results",
  "",
  "| Site | Lat | Lon | Eur exc | Eur onset | Chi exc | Chi onset | Tree exc | Tree onset | BC signal |",
  "|------|-----|-----|---------|-----------|---------|-----------|----------|------------|-----------|"
)

for (i in 1:nrow(results)) {
  r <- results[i, ]
  report_lines <- c(report_lines,
    sprintf("| %s | %.2f | %.2f | %s | %s | %s | %s | %s | %s | %s |",
            substr(r$sitename, 1, 25),
            r$lat, r$lon,
            ifelse(r$eur_exceedance, "YES", ifelse(r$eur_testable, "no", "-")),
            ifelse(is.na(r$eur_onset_bp), "-", as.character(round(r$eur_onset_bp))),
            ifelse(r$chi_exceedance, "YES", ifelse(r$chi_testable, "no", "-")),
            ifelse(is.na(r$chi_onset_bp), "-", as.character(round(r$chi_onset_bp))),
            ifelse(r$tree_exceedance, "YES", ifelse(r$tree_testable, "no", "-")),
            ifelse(is.na(r$tree_onset_bp), "-", as.character(round(r$tree_onset_bp))),
            ifelse(r$bc_signal, "YES", "no")))
}

report_lines <- c(report_lines,
  "",
  "## Interpretation Notes",
  "",
  "- **Baseline**: >5000 BP chosen to predate rice intensification in most Chinese regions",
  "- **European indicators**: Plantago, Rumex, Cerealia-type, Secale -- developed for European pastoral/cereal systems",
  "- **Chinese indicators**: Oryza/Oryza-type, Cannabis, Fagopyrum -- markers of East Asian rice/hemp/buckwheat agriculture",
  "- **Tree taxa**: Universal change detectors (non-diagnostic); expected to show high exceedance everywhere",
  "- **Indicator dependency gap**: difference between max and min detection rates across indicator sets",
  "- This is the third continental test (after Europe and ENA) of the identifiability framework"
)

report_path <- file.path(OUT_DIR, "pilot_china_exceedance.md")
writeLines(report_lines, report_path)
cat(sprintf("  Saved report: %s\n", report_path))

cat("\n=== DONE ===\n")
