#!/usr/bin/env Rscript
# =============================================================================
# Pilot Exceedance Analysis: Eastern North America Pollen Data
# =============================================================================
# Applies Paper 6's Bray-Curtis dissimilarity + exceedance framework to ENA
# Key signals: maize agriculture (~2000-3000 BP), 1492 European contact
#
# Author: ArchEco research agent
# Date: 2026-03-28
# =============================================================================

suppressPackageStartupMessages({
  library(neotoma2)
  library(vegan)
  library(parallel)
})

# ============================================================
# Configuration
# ============================================================
CACHE_FILE <- "/home/ayu/archeco/shared/cache/ena_neotoma2_downloads.rds"
OUT_DIR    <- "/home/ayu/archeco/shared"
N_CORES    <- 4L
BASELINE_CUTOFF <- 3000  # BP: everything older = pre-maize baseline
EXCEEDANCE_SD   <- 2     # threshold: mean + N*SD
CONTACT_YEAR    <- 458   # 1492 CE in cal BP
N_SITES_TARGET  <- 150   # How many sites to analyze

# Functional groups for ENA context
FUNC_GROUPS <- list(
  agricultural  = c("Zea", "Chenopodium", "Amaranthus", "Cucurbita",
                     "Helianthus", "Iva"),
  disturbance   = c("Ambrosia", "Artemisia", "Plantago", "Rumex",
                     "Xanthium", "Poaceae undiff"),
  deciduous     = c("Quercus", "Fagus", "Acer", "Ulmus", "Carya", "Betula",
                     "Fraxinus", "Tilia", "Juglans", "Castanea",
                     "Liquidambar", "Ostrya", "Carpinus"),
  coniferous    = c("Pinus", "Picea", "Tsuga", "Abies", "Juniperus", "Thuja")
)

dir.create(dirname(CACHE_FILE), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Step 1: Download data from Neotoma (with caching)
# ============================================================
cat("=== Step 1: Getting pollen data from Neotoma ===\n")

if (file.exists(CACHE_FILE)) {
  cat("  Loading from cache...\n")
  all_samples <- readRDS(CACHE_FILE)
  cat(sprintf("  Cached: %d rows from %d sites\n",
              nrow(all_samples), length(unique(all_samples$siteid))))
} else {
  cat("  Querying Neotoma for ENA pollen sites...\n")
  sites_obj <- get_sites(loc = c(-95, 25, -65, 50),
                         datasettype = "pollen",
                         limit = 250)
  cat(sprintf("  Found %d sites\n", length(sites_obj)))

  # Download in batches
  n_fetch <- min(N_SITES_TARGET, length(sites_obj))
  cat(sprintf("  Downloading %d sites...\n", n_fetch))

  ds_obj <- get_downloads(sites_obj[1:n_fetch])

  # Extract all samples into one big data frame
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
# Step 2: Pre-filter sites
# ============================================================
cat("\n=== Step 2: Filtering sites ===\n")

# Identify sites with sufficient data:
# - age range spanning baseline (>3000 BP) and test (<3000 BP)
# - at least 10 samples total
site_summary <- do.call(rbind, lapply(split(all_samples, all_samples$siteid), function(df) {
  ages <- unique(df$age)
  ages <- ages[!is.na(ages)]
  data.frame(
    siteid   = df$siteid[1],
    sitename = df$sitename[1],
    lat      = df$lat[1],
    lon      = df$long[1],
    n_ages   = length(ages),
    age_min  = min(ages),
    age_max  = max(ages),
    has_baseline = any(ages > BASELINE_CUTOFF),
    has_test     = any(ages <= BASELINE_CUTOFF & ages >= 0),
    n_baseline   = sum(ages > BASELINE_CUTOFF),
    n_test       = sum(ages <= BASELINE_CUTOFF & ages >= 0),
    has_ambrosia = any(grepl("Ambrosia", df$variablename, ignore.case = TRUE)),
    has_zea      = any(grepl("^Zea", df$variablename, ignore.case = TRUE)),
    stringsAsFactors = FALSE
  )
}))

# Filter
good <- site_summary$has_baseline & site_summary$has_test &
        site_summary$n_baseline >= 3 & site_summary$n_test >= 3 &
        site_summary$n_ages >= 10

good_sites <- site_summary[good, ]
cat(sprintf("  Total sites: %d\n", nrow(site_summary)))
cat(sprintf("  Sites passing quality filter: %d\n", nrow(good_sites)))
cat(sprintf("  With Ambrosia: %d\n", sum(good_sites$has_ambrosia)))
cat(sprintf("  With Zea: %d\n", sum(good_sites$has_zea)))

# ============================================================
# Step 3: Exceedance analysis per site (PARALLELIZED)
# ============================================================
cat("\n=== Step 3: Exceedance analysis ===\n")

classify_taxon <- function(taxon) {
  for (grp in names(FUNC_GROUPS)) {
    for (pattern in FUNC_GROUPS[[grp]]) {
      if (grepl(paste0("^", pattern), taxon, ignore.case = TRUE)) {
        return(grp)
      }
    }
  }
  return("other")
}

analyze_one_site <- function(sid, verbose = FALSE) {
  tryCatch({
    df <- all_samples[all_samples$siteid == sid, ]
    df <- df[!is.na(df$age) & !is.na(df$value) & df$value > 0, ]

    if (nrow(df) < 20) { if(verbose) cat(sprintf("    site %d: <20 rows after filter (%d)\n", sid, nrow(df))); return(NULL) }

    meta <- site_summary[site_summary$siteid == sid, ]

    # Build sample-taxon matrix: rows = unique ages, cols = taxa
    df <- df[!is.na(df$age) & is.finite(df$age), ]
    ages <- sort(unique(df$age))
    ages <- ages[!is.na(ages)]
    taxa <- unique(df$variablename)
    taxa <- taxa[!is.na(taxa)]

    if (length(ages) < 5 || length(taxa) < 5) { if(verbose) cat(sprintf("    site %d: too few ages(%d) or taxa(%d)\n", sid, length(ages), length(taxa))); return(NULL) }

    # Pivot to wide: use simple integer row names to avoid issues
    mat <- matrix(0, nrow = length(ages), ncol = length(taxa),
                  dimnames = list(seq_along(ages), taxa))

    # Vectorized fill using match()
    ai <- match(df$age, ages)
    ti <- match(df$variablename, taxa)
    valid <- !is.na(ai) & !is.na(ti)
    for (r in which(valid)) {
      mat[ai[r], ti[r]] <- mat[ai[r], ti[r]] + df$value[r]
    }

    # Remove non-pollen columns (e.g. "Sample quantity", "Concentration")
    # Keep only taxa with ecologicalgroup containing "TRSH" or "UPHE" or "VACR" etc.
    # Simpler: just remove obvious non-count columns
    bad_cols <- grep("^(Sample|Concentration|Influx|Lycopodium|Spike|Microsphere)",
                     taxa, ignore.case = TRUE)
    if (length(bad_cols) > 0) {
      mat <- mat[, -bad_cols, drop = FALSE]
      taxa <- taxa[-bad_cols]
    }

    # Remove taxa present in fewer than 2 samples
    taxa_presence <- colSums(mat > 0)
    keep <- taxa_presence >= 2
    mat <- mat[, keep, drop = FALSE]
    taxa <- taxa[keep]

    if (ncol(mat) < 5) { if(verbose) cat(sprintf("    site %d: <5 taxa after filter (%d)\n", sid, ncol(mat))); return(NULL) }

    # Convert to proportions
    row_sums <- rowSums(mat)
    zero_rows <- row_sums == 0
    mat_prop <- mat / ifelse(row_sums == 0, 1, row_sums)

    # Split into baseline (>3000 BP) and test periods
    baseline_idx <- which(ages > BASELINE_CUTOFF)
    test_idx     <- which(ages <= BASELINE_CUTOFF & ages >= 0)

    if (length(baseline_idx) < 3 || length(test_idx) < 3) { if(verbose) cat(sprintf("    site %d: insufficient baseline(%d) or test(%d)\n", sid, length(baseline_idx), length(test_idx))); return(NULL) }

    # Baseline centroid (mean proportions)
    baseline_centroid <- colMeans(mat_prop[baseline_idx, , drop = FALSE])

    # BC dissimilarity of each sample from centroid
    combined <- rbind(mat_prop, centroid = baseline_centroid)
    bc_all <- as.matrix(vegdist(combined, method = "bray"))

    centroid_row <- nrow(combined)
    bc_from_centroid <- bc_all[1:(centroid_row - 1), centroid_row]

    # Baseline BC stats
    bc_baseline <- bc_from_centroid[baseline_idx]
    bc_mean <- mean(bc_baseline, na.rm = TRUE)
    bc_sd   <- sd(bc_baseline, na.rm = TRUE)
    if (is.na(bc_sd) || bc_sd == 0) { if(verbose) cat(sprintf("    site %d: bc_sd=0 or NA\n", sid)); return(NULL) }
    threshold <- bc_mean + EXCEEDANCE_SD * bc_sd

    # Test period exceedances
    bc_test <- bc_from_centroid[test_idx]
    test_ages <- ages[test_idx]
    exceedances <- bc_test > threshold

    has_signal <- any(exceedances, na.rm = TRUE)
    onset_age <- if (has_signal) max(test_ages[exceedances], na.rm = TRUE) else NA

    # Max BC in test period
    max_bc <- max(bc_test, na.rm = TRUE)
    max_bc_age <- test_ages[which.max(bc_test)]

    # Functional decomposition: compare recent (<500 BP) to baseline
    recent_idx <- which(ages <= 500 & ages >= 0)
    if (length(recent_idx) < 1) recent_idx <- test_idx[length(test_idx)]

    recent_mean <- colMeans(mat_prop[recent_idx, , drop = FALSE])
    delta <- recent_mean - baseline_centroid

    taxon_groups <- sapply(taxa, classify_taxon)
    group_change <- tapply(abs(delta), taxon_groups, sum, na.rm = TRUE)
    total_change <- sum(abs(delta))
    group_pct <- if (total_change > 0) group_change / total_change * 100 else setNames(rep(0, 5), c("agricultural","coniferous","deciduous","disturbance","other"))

    # Ambrosia analysis
    amb_cols <- grep("^Ambrosia", taxa, ignore.case = TRUE)
    amb_ts <- if (length(amb_cols) > 0) rowSums(mat_prop[, amb_cols, drop = FALSE]) else rep(0, length(ages))

    amb_baseline <- amb_ts[baseline_idx]
    amb_mean <- mean(amb_baseline, na.rm = TRUE)
    amb_sd   <- sd(amb_baseline, na.rm = TRUE)
    amb_thresh <- amb_mean + EXCEEDANCE_SD * ifelse(is.na(amb_sd) || amb_sd == 0, 0.01, amb_sd)

    amb_test <- amb_ts[test_idx]
    amb_exc  <- amb_test > amb_thresh
    ambrosia_onset <- if (any(amb_exc, na.rm = TRUE)) max(test_ages[amb_exc], na.rm = TRUE) else NA

    # Ambrosia near 1492 contact
    near_contact <- test_ages >= 258 & test_ages <= 658
    ambrosia_at_contact <- if (any(near_contact) && length(amb_cols) > 0) {
      any(amb_test[near_contact] > amb_thresh, na.rm = TRUE)
    } else FALSE

    # Zea mays detection
    zea_cols <- grep("^Zea", taxa, ignore.case = TRUE)
    has_zea <- length(zea_cols) > 0 && any(mat[, zea_cols] > 0)
    zea_first_age <- if (has_zea) {
      zea_rows <- rowSums(mat[, zea_cols, drop = FALSE]) > 0
      if (any(zea_rows)) max(ages[zea_rows]) else NA
    } else NA

    # Post-1492 recovery
    post_contact <- which(ages < CONTACT_YEAR & ages >= 0)
    pre_contact  <- which(ages >= CONTACT_YEAR & ages <= 1000)
    recovery_signal <- if (length(post_contact) >= 2 && length(pre_contact) >= 2) {
      bc_post <- mean(bc_from_centroid[post_contact], na.rm = TRUE)
      bc_pre  <- mean(bc_from_centroid[pre_contact], na.rm = TRUE)
      bc_post < bc_pre
    } else NA

    # Classification
    human_pct <- sum(group_pct[c("agricultural", "disturbance")], na.rm = TRUE)
    classification <- if (!has_signal) "no_signal"
                      else if (human_pct > 50) "H1_human"
                      else "H3_natural"

    # Safely extract group percentages (some groups may not exist)
    get_pct <- function(g) { v <- group_pct[g]; if(is.null(v) || is.na(v)) 0 else round(v, 1) }

    data.frame(
      siteid = as.integer(sid),
      sitename = as.character(meta$sitename[1]),
      lat = as.numeric(meta$lat[1]),
      lon = as.numeric(meta$lon[1]),
      n_samples = length(ages),
      n_baseline = length(baseline_idx),
      n_test = length(test_idx),
      n_taxa = ncol(mat),
      bc_mean_baseline = round(bc_mean, 4),
      bc_sd_baseline = round(bc_sd, 4),
      threshold = round(threshold, 4),
      has_signal = has_signal,
      onset_age_bp = onset_age,
      max_bc = round(max_bc, 4),
      max_bc_age = max_bc_age,
      classification = classification,
      pct_agricultural = get_pct("agricultural"),
      pct_disturbance = get_pct("disturbance"),
      pct_deciduous = get_pct("deciduous"),
      pct_coniferous = get_pct("coniferous"),
      pct_other = get_pct("other"),
      ambrosia_onset_bp = ambrosia_onset,
      ambrosia_at_contact = ambrosia_at_contact,
      has_zea = has_zea,
      zea_first_age_bp = zea_first_age,
      recovery_signal = recovery_signal,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    if(verbose) cat(sprintf("    site %d ERROR: %s\n", sid, e$message))
    NULL
  })
}

# Run analysis in parallel
target_sids <- good_sites$siteid
cat(sprintf("  Analyzing %d sites with %d cores...\n", length(target_sids), N_CORES))
t0 <- Sys.time()

# mclapply needs the data in the global env (it forks)
# First pass: sequential with error capture to identify failures
results_list <- list()
n_fail <- 0
for (i in seq_along(target_sids)) {
  sid <- target_sids[i]
  r <- tryCatch(
    analyze_one_site(sid, verbose = TRUE),
    error = function(e) {
      cat(sprintf("  ERROR site %d (%s): %s\n", sid,
                  good_sites$sitename[good_sites$siteid == sid], e$message))
      NULL
    }
  )
  if (!is.null(r)) {
    results_list[[length(results_list) + 1]] <- r
  } else {
    n_fail <- n_fail + 1
  }
}
cat(sprintf("  Succeeded: %d, Failed: %d\n", length(results_list), n_fail))

results <- do.call(rbind, results_list)
rownames(results) <- NULL

t1 <- Sys.time()
cat(sprintf("  Analysis completed in %.1f seconds\n", as.numeric(t1 - t0, units = "secs")))
cat(sprintf("  Sites with valid results: %d / %d\n", nrow(results), length(target_sids)))

# ============================================================
# Step 4: Summary statistics
# ============================================================
cat("\n" , paste(rep("=", 60), collapse=""), "\n")
cat("PILOT ENA EXCEEDANCE ANALYSIS RESULTS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

n_total <- nrow(results)
n_signal <- sum(results$has_signal, na.rm = TRUE)

cat(sprintf("1. SIGNAL DETECTION\n"))
cat(sprintf("   Total sites analyzed: %d\n", n_total))
cat(sprintf("   Sites with BC exceedance signal: %d (%.1f%%)\n",
            n_signal, n_signal / n_total * 100))

signal_sites <- results[results$has_signal == TRUE, ]

cat(sprintf("\n2. ONSET TIMING\n"))
if (nrow(signal_sites) > 0) {
  cat(sprintf("   Median onset age: %.0f BP\n", median(signal_sites$onset_age_bp, na.rm = TRUE)))
  cat(sprintf("   Mean onset age: %.0f BP\n", mean(signal_sites$onset_age_bp, na.rm = TRUE)))
  cat(sprintf("   Range: %d - %d BP\n",
              min(signal_sites$onset_age_bp, na.rm = TRUE),
              max(signal_sites$onset_age_bp, na.rm = TRUE)))

  pre_maize <- sum(signal_sites$onset_age_bp > 2000, na.rm = TRUE)
  maize_era <- sum(signal_sites$onset_age_bp <= 2000 & signal_sites$onset_age_bp > 1000, na.rm = TRUE)
  intensif  <- sum(signal_sites$onset_age_bp <= 1000 & signal_sites$onset_age_bp > CONTACT_YEAR, na.rm = TRUE)
  contact   <- sum(signal_sites$onset_age_bp <= CONTACT_YEAR, na.rm = TRUE)
  cat(sprintf("   Pre-maize (>2000 BP): %d sites\n", pre_maize))
  cat(sprintf("   Maize era (2000-1000 BP): %d sites\n", maize_era))
  cat(sprintf("   Intensification (1000-458 BP): %d sites\n", intensif))
  cat(sprintf("   Post-contact (<458 BP): %d sites\n", contact))
} else {
  cat("   No signals detected.\n")
}

cat(sprintf("\n3. FUNCTIONAL DECOMPOSITION (mean %% across signal sites)\n"))
if (nrow(signal_sites) > 0) {
  cat(sprintf("   Agricultural: %.1f%%\n", mean(signal_sites$pct_agricultural, na.rm = TRUE)))
  cat(sprintf("   Disturbance:  %.1f%%\n", mean(signal_sites$pct_disturbance, na.rm = TRUE)))
  cat(sprintf("   Deciduous:    %.1f%%\n", mean(signal_sites$pct_deciduous, na.rm = TRUE)))
  cat(sprintf("   Coniferous:   %.1f%%\n", mean(signal_sites$pct_coniferous, na.rm = TRUE)))
  cat(sprintf("   Other:        %.1f%%\n", mean(signal_sites$pct_other, na.rm = TRUE)))
}

cat(sprintf("\n4. EXCEEDANCE CLASSIFICATION\n"))
cat(sprintf("   H1 (human-driven): %d (%.1f%%)\n",
            sum(results$classification == "H1_human", na.rm = TRUE),
            sum(results$classification == "H1_human", na.rm = TRUE) / n_total * 100))
cat(sprintf("   H3 (natural):      %d (%.1f%%)\n",
            sum(results$classification == "H3_natural", na.rm = TRUE),
            sum(results$classification == "H3_natural", na.rm = TRUE) / n_total * 100))
cat(sprintf("   No signal:         %d (%.1f%%)\n",
            sum(results$classification == "no_signal", na.rm = TRUE),
            sum(results$classification == "no_signal", na.rm = TRUE) / n_total * 100))

cat(sprintf("\n5. AMBROSIA EXCEEDANCE TIMING\n"))
amb_sites <- results[!is.na(results$ambrosia_onset_bp), ]
if (nrow(amb_sites) > 0) {
  cat(sprintf("   Sites with Ambrosia exceedance: %d (%.1f%%)\n",
              nrow(amb_sites), nrow(amb_sites) / n_total * 100))
  cat(sprintf("   Median onset: %.0f BP\n", median(amb_sites$ambrosia_onset_bp, na.rm = TRUE)))
  amb_precolumbian <- sum(amb_sites$ambrosia_onset_bp > CONTACT_YEAR, na.rm = TRUE)
  amb_postcontact  <- sum(amb_sites$ambrosia_onset_bp <= CONTACT_YEAR, na.rm = TRUE)
  cat(sprintf("   Pre-Columbian onset (>458 BP): %d\n", amb_precolumbian))
  cat(sprintf("   Post-contact onset (<=458 BP): %d\n", amb_postcontact))
  cat(sprintf("   Ambrosia spike near 1492: %d sites (%.1f%%)\n",
              sum(results$ambrosia_at_contact, na.rm = TRUE),
              sum(results$ambrosia_at_contact, na.rm = TRUE) / n_total * 100))
} else {
  cat("   No Ambrosia exceedances detected.\n")
}

cat(sprintf("\n6. ZEA MAYS PRESENCE\n"))
zea_results <- results[results$has_zea == TRUE, ]
cat(sprintf("   Sites with Zea pollen: %d (%.1f%%)\n",
            nrow(zea_results), nrow(zea_results) / n_total * 100))
if (nrow(zea_results) > 0) {
  cat(sprintf("   Median first appearance: %.0f BP\n",
              median(zea_results$zea_first_age_bp, na.rm = TRUE)))
  cat(sprintf("   Range: %d - %d BP\n",
              min(zea_results$zea_first_age_bp, na.rm = TRUE),
              max(zea_results$zea_first_age_bp, na.rm = TRUE)))
}

cat(sprintf("\n7. POST-1492 RECOVERY SIGNAL\n"))
rec_data <- results[!is.na(results$recovery_signal), ]
n_recovery <- 0
if (nrow(rec_data) > 0) {
  n_recovery <- sum(rec_data$recovery_signal, na.rm = TRUE)
  cat(sprintf("   Sites testable for recovery: %d\n", nrow(rec_data)))
  cat(sprintf("   Sites showing recovery (BC decrease after 1492): %d (%.1f%%)\n",
              n_recovery, n_recovery / nrow(rec_data) * 100))
} else {
  cat("   Insufficient temporal resolution to test\n")
}

# ============================================================
# Step 5: Save outputs
# ============================================================
cat("\n=== Saving outputs ===\n")

# CSV
csv_path <- file.path(OUT_DIR, "pilot_ena_exceedance_data.csv")
write.csv(results, csv_path, row.names = FALSE)
cat(sprintf("  Saved CSV: %s\n", csv_path))

# Markdown report
report_lines <- c(
  "# Pilot ENA Exceedance Analysis Results",
  "",
  sprintf("**Date**: %s", Sys.Date()),
  sprintf("**Sites analyzed**: %d (of %d downloaded)", n_total, nrow(site_summary)),
  sprintf("**Framework**: Bray-Curtis dissimilarity exceedance (Paper 6 method)"),
  sprintf("**Baseline**: >%d BP (pre-maize)", BASELINE_CUTOFF),
  sprintf("**Threshold**: mean + %dSD of baseline BC", EXCEEDANCE_SD),
  "",
  "## Key Findings",
  "",
  sprintf("### 1. Signal Detection Rate: %.1f%% (%d/%d sites)",
          n_signal / n_total * 100, n_signal, n_total),
  ""
)

if (nrow(signal_sites) > 0) {
  report_lines <- c(report_lines,
    sprintf("### 2. Onset Timing: Median = %.0f BP",
            median(signal_sites$onset_age_bp, na.rm = TRUE)),
    "",
    "| Period | N sites | % of signal sites |",
    "|--------|---------|-------------------|",
    sprintf("| Pre-maize (>2000 BP) | %d | %.1f%% |", pre_maize, pre_maize/nrow(signal_sites)*100),
    sprintf("| Maize era (2000-1000 BP) | %d | %.1f%% |", maize_era, maize_era/nrow(signal_sites)*100),
    sprintf("| Intensification (1000-458 BP) | %d | %.1f%% |", intensif, intensif/nrow(signal_sites)*100),
    sprintf("| Post-contact (<458 BP) | %d | %.1f%% |", contact, contact/nrow(signal_sites)*100),
    "",
    "### 3. Functional Decomposition (mean % of total compositional change)",
    "",
    "| Category | Mean % |",
    "|----------|--------|",
    sprintf("| Agricultural (Zea, Chenopodium, etc.) | %.1f%% |", mean(signal_sites$pct_agricultural, na.rm=TRUE)),
    sprintf("| Disturbance (Ambrosia, Artemisia, etc.) | %.1f%% |", mean(signal_sites$pct_disturbance, na.rm=TRUE)),
    sprintf("| Deciduous trees | %.1f%% |", mean(signal_sites$pct_deciduous, na.rm=TRUE)),
    sprintf("| Coniferous trees | %.1f%% |", mean(signal_sites$pct_coniferous, na.rm=TRUE)),
    sprintf("| Other | %.1f%% |", mean(signal_sites$pct_other, na.rm=TRUE)),
    ""
  )
}

report_lines <- c(report_lines,
  "### 4. Exceedance Classification",
  "",
  "| Class | N | % |",
  "|-------|---|---|",
  sprintf("| H1 (human-driven) | %d | %.1f%% |",
          sum(results$classification == "H1_human"), sum(results$classification == "H1_human")/n_total*100),
  sprintf("| H3 (natural) | %d | %.1f%% |",
          sum(results$classification == "H3_natural"), sum(results$classification == "H3_natural")/n_total*100),
  sprintf("| No signal | %d | %.1f%% |",
          sum(results$classification == "no_signal"), sum(results$classification == "no_signal")/n_total*100),
  ""
)

if (nrow(amb_sites) > 0) {
  report_lines <- c(report_lines,
    "### 5. Ambrosia Exceedance",
    "",
    sprintf("- Sites with Ambrosia exceedance: %d (%.1f%%)", nrow(amb_sites), nrow(amb_sites)/n_total*100),
    sprintf("- Median Ambrosia onset: %.0f BP", median(amb_sites$ambrosia_onset_bp, na.rm=TRUE)),
    sprintf("- Pre-Columbian onset (>458 BP): %d | Post-contact (<=458 BP): %d",
            amb_precolumbian, amb_postcontact),
    sprintf("- Ambrosia spike near 1492: %d sites", sum(results$ambrosia_at_contact, na.rm=TRUE)),
    ""
  )
}

report_lines <- c(report_lines,
  sprintf("### 6. Zea mays: %d sites with Zea pollen (%.1f%%)",
          nrow(zea_results), nrow(zea_results)/n_total*100),
  ""
)

if (nrow(zea_results) > 0) {
  report_lines <- c(report_lines,
    sprintf("- Median first Zea appearance: %.0f BP", median(zea_results$zea_first_age_bp, na.rm=TRUE)),
    sprintf("- Range: %d - %d BP",
            min(zea_results$zea_first_age_bp, na.rm=TRUE), max(zea_results$zea_first_age_bp, na.rm=TRUE)),
    ""
  )
}

if (nrow(rec_data) > 0) {
  report_lines <- c(report_lines,
    "### 7. Post-1492 Recovery",
    "",
    sprintf("- Testable sites: %d", nrow(rec_data)),
    sprintf("- Recovery signal (BC decrease after contact): %d (%.1f%%)",
            n_recovery, n_recovery/nrow(rec_data)*100),
    ""
  )
}

report_lines <- c(report_lines,
  "",
  "## Interpretation Notes",
  "",
  "- **Baseline**: pre-3000 BP chosen to predate maize agriculture in ENA",
  "- **Onset age**: oldest sample exceeding threshold (first departure from baseline)",
  "- **H1 vs H3**: H1 = agricultural + disturbance taxa > 50% of change; H3 = tree/natural taxa dominate",
  "- **Ambrosia**: key disturbance indicator; sharp rise could indicate either pre-Columbian agriculture or post-European land clearance",
  "- **Recovery**: BC returning toward baseline after 1492 suggests vegetation recovery following ~90% Indigenous population collapse",
  "",
  "## Site-Level Results",
  "",
  "| Site | Lat | Lon | Signal | Onset BP | Class | Ambrosia onset | Zea | Recovery |",
  "|------|-----|-----|--------|----------|-------|----------------|-----|----------|"
)

for (i in 1:nrow(results)) {
  r <- results[i, ]
  report_lines <- c(report_lines,
    sprintf("| %s | %.2f | %.2f | %s | %s | %s | %s | %s | %s |",
            substr(r$sitename, 1, 25),
            r$lat, r$lon,
            ifelse(r$has_signal, "YES", "no"),
            ifelse(is.na(r$onset_age_bp), "-", as.character(round(r$onset_age_bp))),
            r$classification,
            ifelse(is.na(r$ambrosia_onset_bp), "-", as.character(round(r$ambrosia_onset_bp))),
            ifelse(r$has_zea, "YES", "-"),
            ifelse(is.na(r$recovery_signal), "?",
                   ifelse(r$recovery_signal, "YES", "no"))))
}

report_path <- file.path(OUT_DIR, "pilot_ena_exceedance_results.md")
writeLines(report_lines, report_path)
cat(sprintf("  Saved report: %s\n", report_path))

cat("\n=== DONE ===\n")
