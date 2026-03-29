#!/usr/bin/env Rscript
# =============================================================================
# Pilot Exceedance Analysis: Sub-Saharan Africa Pollen Data — Fifth Continental Test
# =============================================================================
# Tests the identifiability framework (Paper 8) on Sub-Saharan African pollen data.
# Prediction: European indicators low (~15-20%), African crop indicators (oil palm,
# sorghum, millet) show positive detection, tree exceedance high (universal).
#
# Key question: Africa's indigenous crops (yam, sorghum, millet, oil palm) have
# very different pollen characteristics from European cereals. Oil palm (Elaeis
# guineensis) is insect-pollinated but produces abundant pollen. Sorghum/millet
# are wind-pollinated grasses — can they be distinguished from wild Poaceae?
#
# Three indicator sets:
#   Set A (European): Plantago, Rumex, Cerealia-type, Secale
#   Set B (African): Elaeis guineensis, Sorghum-type, Pennisetum-type, Poaceae
#   Set C (Tree): Major African tree genera
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
CACHE_FILE <- "/home/ayu/archeco/shared/cache/africa_neotoma2_downloads.rds"
OUT_DIR    <- "/home/ayu/archeco/shared"
SCRIPT_DIR <- "/home/ayu/archeco/shared/scripts"
BASELINE_CUTOFF <- 5000  # BP: >5000 = pre-agriculture baseline
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
  african = list(
    name = "Set B: African agriculture",
    patterns = c("^Elaeis", "Elaeis guineensis", "^Sorghum", "Sorghum.type",
                 "^Pennisetum", "Pennisetum.type", "^Poaceae",
                 "^Cenchrus", "^Eleusine", "^Digitaria",
                 "^Dioscorea", "^Vigna", "^Musa")
  ),
  tree = list(
    name = "Set C: Tree taxa (African)",
    patterns = c("^Quercus", "^Olea", "^Podocarpus", "^Juniperus",
                 "^Combretaceae", "^Terminalia", "^Brachystegia",
                 "^Isoberlinia", "^Julbernardia",
                 "^Moraceae", "^Urticaceae", "^Myrtaceae",
                 "^Celtis", "^Alchornea", "^Macaranga",
                 "^Hagenia", "^Afrocrania",
                 "^Syzygium", "^Rapanea", "^Myrsine",
                 "^Ericaceae", "^Acacia", "^Vachellia",
                 "^Uapaca", "^Sapotaceae", "^Euphorbiaceae",
                 "^Rubiaceae", "^Pinus", "^Eucalyptus",
                 "^Caesalpiniaceae", "^Fabaceae")
  )
)

dir.create(dirname(CACHE_FILE), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Step 1: Download African pollen data from Neotoma
# ============================================================
cat("=== Step 1: Getting pollen data from Neotoma (Sub-Saharan Africa) ===\n")

if (file.exists(CACHE_FILE)) {
  cat("  Loading from cache...\n")
  all_samples <- readRDS(CACHE_FILE)
  cat(sprintf("  Cached: %d rows from %d sites\n",
              nrow(all_samples), length(unique(all_samples$siteid))))
} else {
  cat("  Querying Neotoma for African pollen sites...\n")
  cat("  Bounding box: lon -20 to 55, lat -35 to 15\n")

  # Sub-Saharan Africa bounding box
  # lon -20 (West Africa coast) to 55 (East Africa/Horn), lat -35 (S Africa) to 15 (Sahel)
  sites_obj <- NULL
  for (attempt in 1:10) {
    tryCatch({
      sites_obj <- get_sites(loc = c(-20, -35, 55, 15),
                             datasettype = "pollen",
                             limit = 300)
      break
    }, error = function(e) {
      cat(sprintf("  get_sites attempt %d failed: %s\n", attempt, conditionMessage(e)))
      if (attempt < 10) {
        wait <- min(30 * attempt, 180)
        cat(sprintf("  Waiting %ds before retry...\n", wait))
        Sys.sleep(wait)
      }
    })
  }
  if (is.null(sites_obj)) stop("Failed to get sites after 10 attempts.")
  cat(sprintf("  Found %d sites\n", length(sites_obj)))

  if (length(sites_obj) == 0) {
    stop("No sites found. Check Neotoma connectivity.")
  }

  # Download in batches to avoid 502 errors
  n_sites <- length(sites_obj)
  batch_size <- 25
  n_batches <- ceiling(n_sites / batch_size)
  cat(sprintf("  Downloading %d sites in %d batches of %d...\n",
              n_sites, n_batches, batch_size))

  all_samples_list <- list()

  for (b in seq_len(n_batches)) {
    start_i <- (b - 1) * batch_size + 1
    end_i   <- min(b * batch_size, n_sites)
    cat(sprintf("  Batch %d/%d (sites %d-%d)...\n", b, n_batches, start_i, end_i))

    batch_sites <- sites_obj[start_i:end_i]

    # Retry up to 3 times
    success <- FALSE
    for (attempt in 1:3) {
      tryCatch({
        ds_batch <- get_downloads(batch_sites)
        samp_batch <- samples(ds_batch)
        samp_batch <- samp_batch[samp_batch$datasettype == "pollen", ]
        all_samples_list[[b]] <- samp_batch
        cat(sprintf("    Got %d rows\n", nrow(samp_batch)))
        success <- TRUE
        break
      }, error = function(e) {
        cat(sprintf("    Attempt %d failed: %s\n", attempt, conditionMessage(e)))
        if (attempt < 3) {
          cat("    Waiting 10s before retry...\n")
          Sys.sleep(10)
        }
      })
    }
    if (!success) {
      cat(sprintf("    SKIPPED batch %d after 3 failures\n", b))
    }

    # Brief pause between batches
    if (b < n_batches) Sys.sleep(2)
  }

  all_samples <- do.call(rbind, all_samples_list[!sapply(all_samples_list, is.null)])

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

# Geographic subregions for Africa
good_sites$region <- ifelse(good_sites$lon < 10 & good_sites$lat > -5, "West_Africa",
                     ifelse(good_sites$lon >= 10 & good_sites$lon < 35 & good_sites$lat > -5, "Central_Africa",
                     ifelse(good_sites$lon >= 35 & good_sites$lat > -5, "East_Africa",
                            "Southern_Africa")))

cat("\n  Regional distribution:\n")
for (reg in unique(good_sites$region)) {
  cat(sprintf("    %s: %d sites\n", reg, sum(good_sites$region == reg)))
}

cat("\n  Passing sites:\n")
for (i in 1:min(nrow(good_sites), 40)) {
  s <- good_sites[i, ]
  cat(sprintf("    %s (id=%d): lat=%.2f, lon=%.2f, n_ages=%d, baseline=%d, holocene=%d [%s]\n",
              s$sitename, s$siteid, s$lat, s$lon, s$n_ages, s$n_baseline, s$n_holocene, s$region))
}
if (nrow(good_sites) > 40) cat(sprintf("    ... and %d more\n", nrow(good_sites) - 40))

# ============================================================
# Step 3: Taxon inventory — check what indicators are present
# ============================================================
cat("\n=== Step 3: Taxon inventory ===\n")

good_ids <- good_sites$siteid
good_samples <- all_samples[all_samples$siteid %in% good_ids, ]
all_taxa <- sort(unique(good_samples$variablename))

cat(sprintf("  Total unique taxa across %d sites: %d\n", length(good_ids), length(all_taxa)))

# Special search: Elaeis guineensis (oil palm) — key African crop indicator
cat("\n  === SPECIAL: Elaeis guineensis (oil palm) search ===\n")
elaeis_matches <- all_taxa[grepl("elaeis|oil.palm|palm", all_taxa, ignore.case = TRUE)]
if (length(elaeis_matches) > 0) {
  for (m in elaeis_matches) {
    n_sites <- length(unique(good_samples$siteid[
      grepl(paste0("^", gsub("([.+?^${}()|\\[\\]])", "\\\\\\1", m), "$"),
            good_samples$variablename, ignore.case = TRUE) &
      good_samples$value > 0]))
    cat(sprintf("    [FOUND] %-35s  (%d sites with >0 counts)\n", m, n_sites))
  }
} else {
  cat("    [NONE FOUND] — No Elaeis/palm taxa in dataset\n")
}

# Special search: Sorghum / Pennisetum / millet
cat("\n  === SPECIAL: Sorghum/Pennisetum/millet search ===\n")
cereal_matches <- all_taxa[grepl("sorghum|pennisetum|millet|cenchrus|eleusine|digitaria", all_taxa, ignore.case = TRUE)]
if (length(cereal_matches) > 0) {
  for (m in cereal_matches) {
    n_sites <- length(unique(good_samples$siteid[
      grepl(paste0("^", gsub("([.+?^${}()|\\[\\]])", "\\\\\\1", m), "$"),
            good_samples$variablename, ignore.case = TRUE) &
      good_samples$value > 0]))
    cat(sprintf("    [FOUND] %-35s  (%d sites with >0 counts)\n", m, n_sites))
  }
} else {
  cat("    [NONE FOUND] — No Sorghum/Pennisetum/millet taxa in dataset\n")
}

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
        cat(sprintf("    [FOUND] %-35s  (%d sites)\n", m, n_sites_with))
      }
    }
  }
  if (length(found_taxa) == 0) {
    cat("    [NONE FOUND]\n")
  }
}

# Also list any agriculture-related taxa we might have missed
cat("\n  === Broad agriculture-related taxa scan ===\n")
agri_patterns <- c("elaeis", "palm", "sorghum", "pennisetum", "millet",
                    "dioscorea", "yam", "vigna", "cowpea", "musa", "banana",
                    "plantain", "ensete", "teff", "eragrostis",
                    "gossypium", "cotton", "ricinus", "castor",
                    "coffea", "coffee", "cola", "kola",
                    "manihot", "cassava", "zea", "maize",
                    "arachis", "groundnut", "sesamum", "sesame")
for (pat in agri_patterns) {
  matches <- all_taxa[grepl(pat, all_taxa, ignore.case = TRUE)]
  if (length(matches) > 0) {
    for (m in matches) {
      n_sites_with <- length(unique(good_samples$siteid[
        grepl(paste0("^", gsub("([.+?^${}()|\\[\\]])", "\\\\\\1", m), "$"),
              good_samples$variablename, ignore.case = TRUE) &
        good_samples$value > 0]))
      cat(sprintf("    [AGRI] %-35s  (%d sites)\n", m, n_sites_with))
    }
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
    bad_cols <- grep("^(Sample|Concentration|Influx|Lycopodium|Spike|Microsphere|Charcoal|Carbon|Indeterminate|Unknown|Damaged)",
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
    res_eur  <- exceedance_for_set(sid, INDICATOR_SETS$european$patterns, mat_prop, ages, taxa)
    res_afr  <- exceedance_for_set(sid, INDICATOR_SETS$african$patterns, mat_prop, ages, taxa)
    res_tree <- exceedance_for_set(sid, INDICATOR_SETS$tree$patterns, mat_prop, ages, taxa)

    # --- Special: Elaeis-only exceedance ---
    elaeis_patterns <- c("^Elaeis", "Elaeis guineensis")
    res_elaeis <- exceedance_for_set(sid, elaeis_patterns, mat_prop, ages, taxa)

    # BC dissimilarity exceedance (whole assemblage)
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
      region = as.character(meta$region[1]),
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
      # African indicators
      afr_testable = res_afr$testable,
      afr_exceedance = res_afr$exceedance,
      afr_onset_bp = res_afr$onset_bp,
      afr_n_taxa = res_afr$n_taxa_matched,
      afr_taxa = res_afr$taxa_matched,
      # Elaeis-only
      elaeis_testable = res_elaeis$testable,
      elaeis_exceedance = res_elaeis$exceedance,
      elaeis_onset_bp = res_elaeis$onset_bp,
      elaeis_n_taxa = res_elaeis$n_taxa_matched,
      elaeis_taxa = res_elaeis$taxa_matched,
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
cat("PILOT SUB-SAHARAN AFRICA EXCEEDANCE ANALYSIS RESULTS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

n_total <- nrow(results)

# European indicators
n_eur_testable <- sum(results$eur_testable)
n_eur_exc      <- sum(results$eur_exceedance)
eur_rate       <- if (n_eur_testable > 0) n_eur_exc / n_eur_testable * 100 else 0
eur_rate_total <- n_eur_exc / n_total * 100

# African indicators
n_afr_testable <- sum(results$afr_testable)
n_afr_exc      <- sum(results$afr_exceedance)
afr_rate       <- if (n_afr_testable > 0) n_afr_exc / n_afr_testable * 100 else 0
afr_rate_total <- n_afr_exc / n_total * 100

# Elaeis-only
n_elaeis_testable <- sum(results$elaeis_testable)
n_elaeis_exc      <- sum(results$elaeis_exceedance)
elaeis_rate       <- if (n_elaeis_testable > 0) n_elaeis_exc / n_elaeis_testable * 100 else 0
elaeis_rate_total <- n_elaeis_exc / n_total * 100

# Tree indicators
n_tree_testable <- sum(results$tree_testable)
n_tree_exc      <- sum(results$tree_exceedance)
tree_rate       <- if (n_tree_testable > 0) n_tree_exc / n_tree_testable * 100 else 0
tree_rate_total <- n_tree_exc / n_total * 100

# BC whole-assemblage
n_bc_signal <- sum(results$bc_signal)
bc_rate     <- n_bc_signal / n_total * 100

cat(sprintf("1. INDICATOR SET EXCEEDANCE RATES (N=%d sites)\n\n", n_total))
cat(sprintf("   %-40s  Testable  Exceedance  Rate(testable)  Rate(total)\n", "Indicator Set"))
cat(sprintf("   %-40s  --------  ----------  --------------  -----------\n", ""))
cat(sprintf("   %-40s  %4d      %4d        %5.1f%%          %5.1f%%\n",
            "Set A: European (Behre 1981)", n_eur_testable, n_eur_exc, eur_rate, eur_rate_total))
cat(sprintf("   %-40s  %4d      %4d        %5.1f%%          %5.1f%%\n",
            "Set B: African agriculture", n_afr_testable, n_afr_exc, afr_rate, afr_rate_total))
cat(sprintf("   %-40s  %4d      %4d        %5.1f%%          %5.1f%%\n",
            "  B-sub: Elaeis guineensis ONLY", n_elaeis_testable, n_elaeis_exc, elaeis_rate, elaeis_rate_total))
cat(sprintf("   %-40s  %4d      %4d        %5.1f%%          %5.1f%%\n",
            "Set C: Tree taxa", n_tree_testable, n_tree_exc, tree_rate, tree_rate_total))
cat(sprintf("   %-40s  %4d      %4d        %5.1f%%          %5.1f%%\n",
            "Whole-assemblage BC", n_total, n_bc_signal, bc_rate, bc_rate))

# Regional breakdown
cat(sprintf("\n2. REGIONAL BREAKDOWN\n\n"))
for (reg in c("West_Africa", "Central_Africa", "East_Africa", "Southern_Africa")) {
  reg_data <- results[results$region == reg, ]
  if (nrow(reg_data) == 0) next
  n_reg <- nrow(reg_data)
  cat(sprintf("   %s (n=%d):\n", gsub("_", " ", reg), n_reg))
  cat(sprintf("     European: %d/%d (%.1f%%)\n",
              sum(reg_data$eur_exceedance), n_reg, sum(reg_data$eur_exceedance)/n_reg*100))
  cat(sprintf("     African: %d/%d (%.1f%%)\n",
              sum(reg_data$afr_exceedance), n_reg, sum(reg_data$afr_exceedance)/n_reg*100))
  cat(sprintf("     Elaeis only: %d/%d (%.1f%%)\n",
              sum(reg_data$elaeis_exceedance), n_reg, sum(reg_data$elaeis_exceedance)/n_reg*100))
  cat(sprintf("     Tree: %d/%d (%.1f%%)\n",
              sum(reg_data$tree_exceedance), n_reg, sum(reg_data$tree_exceedance)/n_reg*100))
  cat(sprintf("     BC: %d/%d (%.1f%%)\n",
              sum(reg_data$bc_signal), n_reg, sum(reg_data$bc_signal)/n_reg*100))
}

cat(sprintf("\n3. ONSET TIMING (among exceedance sites)\n\n"))

for (set_label in c("European", "African", "Elaeis-only", "Tree")) {
  col <- switch(set_label,
                European     = "eur_onset_bp",
                African      = "afr_onset_bp",
                `Elaeis-only` = "elaeis_onset_bp",
                Tree         = "tree_onset_bp")
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

cat(sprintf("\n4. INDICATOR DEPENDENCY GAP\n"))
rates <- c(eur_rate_total, afr_rate_total, tree_rate_total)
gap <- max(rates) - min(rates)
labels <- c("European", "African", "Tree")
cat(sprintf("   Max detection: %.1f%% (%s)\n",
            max(rates), labels[which.max(rates)]))
cat(sprintf("   Min detection: %.1f%% (%s)\n",
            min(rates), labels[which.min(rates)]))
cat(sprintf("   Gap: %.1f percentage points\n", gap))

cat(sprintf("\n5. KEY QUESTION: ELAEIS GUINEENSIS (OIL PALM) DETECTION\n\n"))
cat(sprintf("   Elaeis guineensis is the signature West African cultivar.\n"))
cat(sprintf("   Though insect-pollinated, it produces abundant pollen and is\n"))
cat(sprintf("   morphologically distinctive (triporate, reticulate).\n"))
cat(sprintf("   Expansion of oil palm in pollen records is a classic marker of\n"))
cat(sprintf("   early West African agriculture and forest disturbance.\n\n"))
cat(sprintf("   Elaeis testable at: %d sites (of %d total)\n", n_elaeis_testable, n_total))
cat(sprintf("   Elaeis exceedance: %d sites (%.1f%% of testable, %.1f%% of total)\n",
            n_elaeis_exc, elaeis_rate, elaeis_rate_total))
if (n_elaeis_exc > 0) {
  elaeis_onsets <- results$elaeis_onset_bp[!is.na(results$elaeis_onset_bp)]
  cat(sprintf("   Earliest Elaeis onset: %d BP\n", max(elaeis_onsets)))
  cat(sprintf("   Median Elaeis onset: %.0f BP\n", median(elaeis_onsets)))
}

cat(sprintf("\n6. AFRICAN CROP POLLEN CHALLENGE\n\n"))
cat(sprintf("   Unlike European cereals (large, distinctive pollen), African staple crops pose\n"))
cat(sprintf("   identification challenges:\n"))
cat(sprintf("   - Sorghum/millet: Poaceae pollen, hard to distinguish from wild grasses\n"))
cat(sprintf("   - Yam (Dioscorea): Very low pollen production, rarely preserved\n"))
cat(sprintf("   - Oil palm (Elaeis): Best diagnostic — distinctive morphology, abundant pollen\n"))
cat(sprintf("   This makes Set B heavily reliant on Poaceae (indirect proxy) and Elaeis.\n"))

cat(sprintf("\n7. CROSS-CONTINENTAL COMPARISON (5 continents)\n\n"))
cat(sprintf("   | Indicator     | Europe | ENA      | China  | LatAm    | Africa (observed) |\n"))
cat(sprintf("   |---------------|--------|----------|--------|----------|-------------------|\n"))
cat(sprintf("   | European      | 87%%    | 16.1%%   | 15.6%% | 41%%     | **%.1f%%**         |\n", eur_rate_total))
cat(sprintf("   | Region-spec.  | --     | 26.8%%   | 13.3%% | 41%%     | **%.1f%%**         |\n", afr_rate_total))
cat(sprintf("   | Elaeis only   | --     | --       | --     | --       | **%.1f%%**         |\n", elaeis_rate_total))
cat(sprintf("   | Tree          | 93.3%% | 94%%     | 91.1%% | 88%%     | **%.1f%%**         |\n", tree_rate_total))

# ============================================================
# Step 6: Save outputs
# ============================================================
cat("\n=== Saving outputs ===\n")

# CSV
csv_path <- file.path(OUT_DIR, "pilot_africa_exceedance_data.csv")
write.csv(results, csv_path, row.names = FALSE)
cat(sprintf("  Saved CSV: %s\n", csv_path))

# Markdown report
report_lines <- c(
  "# Pilot Sub-Saharan Africa Exceedance Analysis Results",
  "",
  sprintf("**Date**: %s", Sys.Date()),
  sprintf("**Sites analyzed**: %d (of %d downloaded)", n_total, nrow(site_summary)),
  "**Framework**: Paper 8 indicator-dependent exceedance",
  sprintf("**Baseline**: >%d BP (pre-agriculture baseline)", BASELINE_CUTOFF),
  sprintf("**Threshold**: mean + %dSD of baseline per taxon", EXCEEDANCE_SD),
  "**Region**: Sub-Saharan Africa (lon -20 to 55, lat -35 to 15)",
  "",
  "## Summary: Indicator Set Exceedance Rates",
  "",
  "| Indicator Set | N Testable | N Exceedance | Rate (testable) | Rate (total) |",
  "|---|---|---|---|---|",
  sprintf("| Set A: European (Behre 1981) | %d | %d | %.1f%% | %.1f%% |",
          n_eur_testable, n_eur_exc, eur_rate, eur_rate_total),
  sprintf("| Set B: African agriculture | %d | %d | %.1f%% | %.1f%% |",
          n_afr_testable, n_afr_exc, afr_rate, afr_rate_total),
  sprintf("| B-sub: Elaeis guineensis ONLY | %d | %d | %.1f%% | %.1f%% |",
          n_elaeis_testable, n_elaeis_exc, elaeis_rate, elaeis_rate_total),
  sprintf("| Set C: Tree taxa | %d | %d | %.1f%% | %.1f%% |",
          n_tree_testable, n_tree_exc, tree_rate, tree_rate_total),
  sprintf("| Whole-assemblage BC | %d | %d | %.1f%% | %.1f%% |",
          n_total, n_bc_signal, bc_rate, bc_rate),
  "",
  "## Cross-Continental Comparison (5 continents)",
  "",
  "| Indicator | Europe | ENA | China | LatAm | Africa |",
  "|---|---|---|---|---|---|",
  sprintf("| European | 87%% | 16.1%% | 15.6%% | 41%% | **%.1f%%** |", eur_rate_total),
  sprintf("| Region-specific | -- | 26.8%% (EAC) | 13.3%% (rice) | 41%% (Meso) | **%.1f%%** (Afr) |", afr_rate_total),
  sprintf("| Signature crop | -- | -- | -- | Zea | **Elaeis %.1f%%** |", elaeis_rate_total),
  sprintf("| Tree | 93.3%% | 94%% | 91.1%% | 88%% | **%.1f%%** |", tree_rate_total),
  "",
  sprintf("## Indicator Dependency Gap: %.1f pp", gap),
  "",
  sprintf("- Max detection: %.1f%% (%s)", max(rates), labels[which.max(rates)]),
  sprintf("- Min detection: %.1f%% (%s)", min(rates), labels[which.min(rates)]),
  "",
  "## Regional Breakdown",
  ""
)

for (reg in c("West_Africa", "Central_Africa", "East_Africa", "Southern_Africa")) {
  reg_data <- results[results$region == reg, ]
  if (nrow(reg_data) == 0) next
  n_reg <- nrow(reg_data)
  report_lines <- c(report_lines,
    sprintf("### %s (n=%d)", gsub("_", " ", reg), n_reg),
    sprintf("- European: %d/%d (%.1f%%)", sum(reg_data$eur_exceedance), n_reg, sum(reg_data$eur_exceedance)/n_reg*100),
    sprintf("- African: %d/%d (%.1f%%)", sum(reg_data$afr_exceedance), n_reg, sum(reg_data$afr_exceedance)/n_reg*100),
    sprintf("- Elaeis only: %d/%d (%.1f%%)", sum(reg_data$elaeis_exceedance), n_reg, sum(reg_data$elaeis_exceedance)/n_reg*100),
    sprintf("- Tree: %d/%d (%.1f%%)", sum(reg_data$tree_exceedance), n_reg, sum(reg_data$tree_exceedance)/n_reg*100),
    sprintf("- BC: %d/%d (%.1f%%)", sum(reg_data$bc_signal), n_reg, sum(reg_data$bc_signal)/n_reg*100),
    "")
}

report_lines <- c(report_lines,
  "## Onset Timing",
  ""
)

for (set_label in c("European", "African", "Elaeis-only", "Tree")) {
  col <- switch(set_label, European = "eur_onset_bp", African = "afr_onset_bp",
                `Elaeis-only` = "elaeis_onset_bp", Tree = "tree_onset_bp")
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

# Key question section
report_lines <- c(report_lines,
  "## Key Question: Elaeis guineensis (Oil Palm) Detection",
  "",
  "Elaeis guineensis is the critical test case for the identifiability framework in Africa:",
  "- **Morphologically distinctive**: Triporate, reticulate pollen not confused with other palms",
  "- **Abundant pollen production**: Despite insect-pollination, produces large quantities",
  "- **Classic agriculture marker**: Oil palm expansion in West African pollen records is well-documented",
  "- **Cultivation vs. wild**: Both wild and cultivated forms present — expansion indicates human management",
  "",
  sprintf("**Result**: Elaeis testable at %d sites, exceedance at %d (%.1f%% of testable, %.1f%% of total)",
          n_elaeis_testable, n_elaeis_exc, elaeis_rate, elaeis_rate_total),
  "",
  "## African Crop Pollen Challenge",
  "",
  "African staple crops pose unique pollen identification challenges:",
  "- **Sorghum/millet (Poaceae)**: Cannot be reliably distinguished from wild grasses at family level",
  "- **Yam (Dioscorea)**: Extremely low pollen production, insect-pollinated, rarely in records",
  "- **Oil palm (Elaeis)**: Best diagnostic crop indicator for Africa",
  "- **Implication**: Set B is heavily reliant on Poaceae (indirect proxy) and Elaeis",
  "- This taxonomic resolution problem is analogous to Oryza-type in China (Poaceae confusion)",
  ""
)

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

report_lines <- c(report_lines, "", "### African indicator taxa found", "")
afr_found <- unique(unlist(lapply(results$afr_taxa[results$afr_taxa != ""], function(x) strsplit(x, "; ")[[1]])))
if (length(afr_found) > 0) {
  for (t in afr_found) report_lines <- c(report_lines, sprintf("- %s", t))
} else {
  report_lines <- c(report_lines, "- None found")
}

report_lines <- c(report_lines, "", "### Elaeis taxa found", "")
elaeis_found <- unique(unlist(lapply(results$elaeis_taxa[results$elaeis_taxa != ""], function(x) strsplit(x, "; ")[[1]])))
if (length(elaeis_found) > 0) {
  for (t in elaeis_found) report_lines <- c(report_lines, sprintf("- %s", t))
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
  "| Site | Region | Lat | Lon | Eur | Afr | Elaeis | Tree | BC |",
  "|------|--------|-----|-----|-----|-----|--------|------|----|"
)

for (i in 1:nrow(results)) {
  r <- results[i, ]
  report_lines <- c(report_lines,
    sprintf("| %s | %s | %.1f | %.1f | %s | %s | %s | %s | %s |",
            substr(r$sitename, 1, 20),
            substr(r$region, 1, 8),
            r$lat, r$lon,
            ifelse(r$eur_exceedance, "YES", ifelse(r$eur_testable, "no", "-")),
            ifelse(r$afr_exceedance, "YES", ifelse(r$afr_testable, "no", "-")),
            ifelse(r$elaeis_exceedance, "YES", ifelse(r$elaeis_testable, "no", "-")),
            ifelse(r$tree_exceedance, "YES", ifelse(r$tree_testable, "no", "-")),
            ifelse(r$bc_signal, "YES", "no")))
}

report_lines <- c(report_lines,
  "",
  "## Interpretation Notes",
  "",
  "- **Baseline**: >5000 BP chosen to predate intensified agriculture in most African regions",
  "  (note: early oil palm management in West Africa may date to ~4000-5000 BP)",
  "- **European indicators**: Plantago, Rumex, Cerealia-type, Secale -- developed for European pastoral/cereal systems",
  "- **African indicators**: Elaeis guineensis, Sorghum-type, Pennisetum-type, Poaceae, Dioscorea, Vigna, Musa",
  "- **Elaeis sub-analysis**: Isolated test of the most diagnostic African crop pollen",
  "- **Tree taxa**: African forest/savanna genera expected to show high exceedance (universal change detectors)",
  "- **Indicator dependency gap**: difference between max and min detection rates across indicator sets",
  "- This is the fifth continental test (after Europe, ENA, China, LatAm) of the identifiability framework",
  "- **Neotoma coverage caveat**: African pollen database coverage is historically sparse compared to",
  "  Europe and North America — results may be biased toward East African highlands where most coring occurred"
)

report_path <- file.path(OUT_DIR, "pilot_africa_exceedance.md")
writeLines(report_lines, report_path)
cat(sprintf("  Saved report: %s\n", report_path))

cat("\n=== DONE ===\n")
