#!/usr/bin/env Rscript
# =============================================================================
# Full Exceedance Analysis: Eastern North America Pollen Data
# =============================================================================
# Expands pilot (96 sites) to full dataset with special focus on:
# - All 33 known Zea mays sites
# - Detailed Ambrosia two-phase analysis
# - Pre-Columbian vs post-contact agriculture comparison
# - Post-1492 vegetation recovery test
#
# Author: ArchEco research agent
# Date: 2026-03-28
# =============================================================================

suppressPackageStartupMessages({
  library(neotoma2)
  library(vegan)
  library(parallel)
  library(jsonlite)
})

# ============================================================
# Configuration
# ============================================================
CACHE_FILE <- "/home/ayu/archeco/shared/cache/ena_neotoma2_downloads.rds"
ZEA_FILE   <- "/home/ayu/archeco/shared/zea_mays_sites.json"
OUT_DIR    <- "/home/ayu/archeco/shared"
N_CORES    <- 4L
BASELINE_CUTOFF <- 3000  # BP: everything older = pre-maize baseline
EXCEEDANCE_SD   <- 2     # threshold: mean + N*SD
CONTACT_YEAR    <- 458   # 1492 CE in cal BP
EUROPEAN_AGRI   <- 200   # ~1750 CE in cal BP — European agricultural expansion

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
# Step 1: Download data — merge cache + new Zea sites + extra ENA sites
# ============================================================
cat("=== Step 1: Getting pollen data from Neotoma ===\n")

# Load existing cache
cached_samples <- readRDS(CACHE_FILE)
cached_ids <- unique(cached_samples$siteid)
cat(sprintf("  Cached: %d rows from %d sites\n", nrow(cached_samples), length(cached_ids)))

# Identify Zea mays sites that need downloading
zea_data <- fromJSON(ZEA_FILE)
zea_ids <- as.integer(names(zea_data))
missing_zea <- zea_ids[!zea_ids %in% cached_ids]
cat(sprintf("  Zea mays sites: %d total, %d already cached, %d to download\n",
            length(zea_ids), sum(zea_ids %in% cached_ids), length(missing_zea)))

# Download missing Zea sites in batches
new_samples_list <- list()

if (length(missing_zea) > 0) {
  batch_size <- 10
  n_batches <- ceiling(length(missing_zea) / batch_size)

  for (b in seq_len(n_batches)) {
    start <- (b - 1) * batch_size + 1
    end <- min(b * batch_size, length(missing_zea))
    batch_ids <- missing_zea[start:end]

    cat(sprintf("  Downloading Zea batch %d/%d (sites: %s)...\n",
                b, n_batches, paste(batch_ids, collapse=",")))

    tryCatch({
      sites_obj <- get_sites(batch_ids)
      if (length(sites_obj) > 0) {
        ds_obj <- get_downloads(sites_obj)
        samp <- samples(ds_obj)
        samp <- samp[samp$datasettype == "pollen", ]
        if (nrow(samp) > 0) {
          new_samples_list[[length(new_samples_list) + 1]] <- samp
          cat(sprintf("    Got %d rows from %d sites\n", nrow(samp), length(unique(samp$siteid))))
        }
      }
    }, error = function(e) {
      cat(sprintf("    Batch %d failed: %s\n", b, e$message))
      # Try individual downloads
      for (sid in batch_ids) {
        tryCatch({
          s <- get_sites(sid)
          if (length(s) > 0) {
            d <- get_downloads(s)
            samp <- samples(d)
            samp <- samp[samp$datasettype == "pollen", ]
            if (nrow(samp) > 0) {
              new_samples_list[[length(new_samples_list) + 1]] <- samp
              cat(sprintf("      Individual site %d: %d rows\n", sid, nrow(samp)))
            }
          }
        }, error = function(e2) {
          cat(sprintf("      Site %d failed: %s\n", sid, e2$message))
        })
      }
    })

    Sys.sleep(1)  # Be polite to the API
  }
}

# Also try to get additional ENA sites beyond the original 150
cat("  Querying Neotoma for additional ENA sites...\n")
tryCatch({
  # Query broadly, hoping to get sites beyond the first 150
  extra_sites <- get_sites(loc = c(-95, 25, -65, 50),
                           datasettype = "pollen",
                           limit = 300)
  cat(sprintf("  Query returned %d sites\n", length(extra_sites)))

  # Get site IDs from the object
  extra_site_ids <- as.numeric(getids(extra_sites))
  new_extra <- extra_site_ids[!extra_site_ids %in% cached_ids & !extra_site_ids %in% missing_zea]

  if (length(new_extra) > 0) {
    # Download up to 100 new sites
    n_new <- min(100, length(new_extra))
    cat(sprintf("  Downloading %d new ENA sites...\n", n_new))

    # Download in batches
    for (b in seq_len(ceiling(n_new / 20))) {
      start <- (b - 1) * 20 + 1
      end <- min(b * 20, n_new)
      batch_idx <- start:end

      tryCatch({
        batch_sites <- extra_sites[which(extra_site_ids %in% new_extra[batch_idx])]
        if (length(batch_sites) > 0) {
          ds <- get_downloads(batch_sites)
          samp <- samples(ds)
          samp <- samp[samp$datasettype == "pollen", ]
          if (nrow(samp) > 0) {
            new_samples_list[[length(new_samples_list) + 1]] <- samp
            cat(sprintf("    Extra batch %d: %d rows from %d sites\n",
                        b, nrow(samp), length(unique(samp$siteid))))
          }
        }
      }, error = function(e) {
        cat(sprintf("    Extra batch %d failed: %s\n", b, e$message))
      })

      Sys.sleep(1)
    }
  }
}, error = function(e) {
  cat(sprintf("  Extra site query failed: %s (continuing with cache + Zea)\n", e$message))
})

# Merge all data
if (length(new_samples_list) > 0) {
  # Harmonize columns before rbind
  new_samples <- do.call(rbind, new_samples_list)

  # Ensure column compatibility
  common_cols <- intersect(names(cached_samples), names(new_samples))
  all_samples <- rbind(cached_samples[, common_cols], new_samples[, common_cols])

  # Update cache
  saveRDS(all_samples, CACHE_FILE)
  cat(sprintf("  Updated cache: %d rows from %d sites\n",
              nrow(all_samples), length(unique(all_samples$siteid))))
} else {
  all_samples <- cached_samples
  cat("  No new sites downloaded; using cache only.\n")
}

cat(sprintf("  TOTAL: %d rows, %d unique sites\n",
            nrow(all_samples), length(unique(all_samples$siteid))))

# ============================================================
# Step 2: Pre-filter sites
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
    age_min  = if(length(ages)>0) min(ages) else NA,
    age_max  = if(length(ages)>0) max(ages) else NA,
    has_baseline = any(ages > BASELINE_CUTOFF),
    has_test     = any(ages <= BASELINE_CUTOFF & ages >= 0),
    n_baseline   = sum(ages > BASELINE_CUTOFF),
    n_test       = sum(ages <= BASELINE_CUTOFF & ages >= 0),
    has_ambrosia = any(grepl("Ambrosia", df$variablename, ignore.case = TRUE)),
    has_zea      = any(grepl("^Zea", df$variablename, ignore.case = TRUE)),
    is_zea_target = df$siteid[1] %in% zea_ids,
    stringsAsFactors = FALSE
  )
}))

# Quality filter: baseline + test period, sufficient samples
good <- site_summary$has_baseline & site_summary$has_test &
        site_summary$n_baseline >= 3 & site_summary$n_test >= 3 &
        site_summary$n_ages >= 10

# Also include Zea target sites with relaxed criteria
zea_relaxed <- site_summary$is_zea_target & site_summary$has_test & site_summary$n_ages >= 5

good_sites <- site_summary[good | zea_relaxed, ]
cat(sprintf("  Total sites in data: %d\n", nrow(site_summary)))
cat(sprintf("  Sites passing standard filter: %d\n", sum(good)))
cat(sprintf("  Zea target sites (relaxed filter): %d\n", sum(zea_relaxed & !good)))
cat(sprintf("  Total sites for analysis: %d\n", nrow(good_sites)))
cat(sprintf("  With Ambrosia: %d\n", sum(good_sites$has_ambrosia)))
cat(sprintf("  With Zea pollen: %d\n", sum(good_sites$has_zea)))
cat(sprintf("  Zea target sites in analysis: %d\n", sum(good_sites$is_zea_target)))

# ============================================================
# Step 3: Core exceedance analysis per site
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

    if (nrow(df) < 10) return(NULL)

    meta <- site_summary[site_summary$siteid == sid, ]

    df <- df[!is.na(df$age) & is.finite(df$age), ]
    ages <- sort(unique(df$age))
    ages <- ages[!is.na(ages)]
    taxa <- unique(df$variablename)
    taxa <- taxa[!is.na(taxa)]

    if (length(ages) < 5 || length(taxa) < 5) return(NULL)

    # Build sample-taxon matrix
    mat <- matrix(0, nrow = length(ages), ncol = length(taxa),
                  dimnames = list(seq_along(ages), taxa))

    ai <- match(df$age, ages)
    ti <- match(df$variablename, taxa)
    valid <- !is.na(ai) & !is.na(ti)
    for (r in which(valid)) {
      mat[ai[r], ti[r]] <- mat[ai[r], ti[r]] + df$value[r]
    }

    # Remove non-pollen columns
    bad_cols <- grep("^(Sample|Concentration|Influx|Lycopodium|Spike|Microsphere)",
                     taxa, ignore.case = TRUE)
    if (length(bad_cols) > 0) {
      mat <- mat[, -bad_cols, drop = FALSE]
      taxa <- taxa[-bad_cols]
    }

    # Remove rare taxa
    taxa_presence <- colSums(mat > 0)
    keep <- taxa_presence >= 2
    mat <- mat[, keep, drop = FALSE]
    taxa <- taxa[keep]

    if (ncol(mat) < 5) return(NULL)

    # Proportions
    row_sums <- rowSums(mat)
    mat_prop <- mat / ifelse(row_sums == 0, 1, row_sums)

    # Split baseline/test
    baseline_idx <- which(ages > BASELINE_CUTOFF)
    test_idx     <- which(ages <= BASELINE_CUTOFF & ages >= 0)

    # For Zea target sites, allow smaller baseline
    min_baseline <- if (sid %in% zea_ids) 2 else 3
    min_test <- if (sid %in% zea_ids) 2 else 3

    if (length(baseline_idx) < min_baseline || length(test_idx) < min_test) return(NULL)

    # Baseline centroid
    baseline_centroid <- colMeans(mat_prop[baseline_idx, , drop = FALSE])

    # BC dissimilarity from centroid
    combined <- rbind(mat_prop, centroid = baseline_centroid)
    bc_all <- as.matrix(vegdist(combined, method = "bray"))
    centroid_row <- nrow(combined)
    bc_from_centroid <- bc_all[1:(centroid_row - 1), centroid_row]

    # Baseline stats
    bc_baseline <- bc_from_centroid[baseline_idx]
    bc_mean <- mean(bc_baseline, na.rm = TRUE)
    bc_sd   <- sd(bc_baseline, na.rm = TRUE)
    if (is.na(bc_sd) || bc_sd == 0) return(NULL)
    threshold <- bc_mean + EXCEEDANCE_SD * bc_sd

    # Test period exceedances
    bc_test <- bc_from_centroid[test_idx]
    test_ages <- ages[test_idx]
    exceedances <- bc_test > threshold

    has_signal <- any(exceedances, na.rm = TRUE)
    onset_age <- if (has_signal) max(test_ages[exceedances], na.rm = TRUE) else NA
    max_bc <- max(bc_test, na.rm = TRUE)
    max_bc_age <- test_ages[which.max(bc_test)]

    # ---- Functional decomposition ----
    recent_idx <- which(ages <= 500 & ages >= 0)
    if (length(recent_idx) < 1) recent_idx <- test_idx[length(test_idx)]

    recent_mean <- colMeans(mat_prop[recent_idx, , drop = FALSE])
    delta <- recent_mean - baseline_centroid

    taxon_groups <- sapply(taxa, classify_taxon)
    group_change <- tapply(abs(delta), taxon_groups, sum, na.rm = TRUE)
    total_change <- sum(abs(delta))
    group_pct <- if (total_change > 0) group_change / total_change * 100 else setNames(rep(0, 5), c("agricultural","coniferous","deciduous","disturbance","other"))

    get_pct <- function(g) { v <- group_pct[g]; if(is.null(v) || is.na(v)) 0 else round(v, 1) }

    # ---- AMBROSIA DETAILED ANALYSIS (Part 3) ----
    amb_cols <- grep("^Ambrosia", taxa, ignore.case = TRUE)
    amb_ts <- if (length(amb_cols) > 0) rowSums(mat_prop[, amb_cols, drop = FALSE]) else rep(0, length(ages))

    amb_baseline <- amb_ts[baseline_idx]
    amb_mean <- mean(amb_baseline, na.rm = TRUE)
    amb_sd   <- sd(amb_baseline, na.rm = TRUE)
    if (is.na(amb_sd) || amb_sd == 0) amb_sd <- 0.01

    # Criterion 1: Ambrosia > 5%
    amb_5pct <- amb_ts > 0.05
    # Criterion 2: Ambrosia > baseline + 3SD
    amb_3sd_thresh <- amb_mean + 3 * amb_sd
    amb_3sd <- amb_ts > amb_3sd_thresh
    # Explosion = either criterion
    amb_explosion <- amb_5pct | amb_3sd

    # Standard 2SD exceedance
    amb_2sd_thresh <- amb_mean + 2 * amb_sd
    amb_exc_2sd <- amb_ts > amb_2sd_thresh

    # Find ambrosia onset (2SD threshold, test period only)
    amb_test_exc <- amb_exc_2sd[test_idx]
    ambrosia_onset <- if (any(amb_test_exc, na.rm = TRUE)) max(test_ages[amb_test_exc], na.rm = TRUE) else NA

    # Find explosion point (any age)
    amb_explosion_ages <- ages[amb_explosion]
    amb_explosion_onset <- if (length(amb_explosion_ages) > 0) max(amb_explosion_ages) else NA

    # Classify Ambrosia: pre-contact vs post-contact
    amb_precontact <- !is.na(amb_explosion_onset) && amb_explosion_onset > CONTACT_YEAR
    amb_postcontact <- !is.na(amb_explosion_onset) && amb_explosion_onset <= CONTACT_YEAR

    # Two-phase Ambrosia detection (for post-contact sites)
    # Phase 1: initial rise; Phase 2: secondary European agricultural rise ~200 BP
    amb_phase2 <- FALSE
    amb_phase2_age <- NA
    if (!is.na(ambrosia_onset) && ambrosia_onset <= CONTACT_YEAR) {
      # Look for a second rise after ~200 BP
      very_recent_idx <- which(ages <= EUROPEAN_AGRI & ages >= 0)
      if (length(very_recent_idx) >= 2) {
        # Is there a significant increase in Ambrosia in the most recent period vs. 458-200 BP?
        mid_idx <- which(ages <= CONTACT_YEAR & ages > EUROPEAN_AGRI)
        if (length(mid_idx) >= 1) {
          amb_mid <- mean(amb_ts[mid_idx], na.rm = TRUE)
          amb_recent <- mean(amb_ts[very_recent_idx], na.rm = TRUE)
          if (amb_recent > amb_mid * 1.5 && amb_recent > 0.03) {
            amb_phase2 <- TRUE
            # Find the onset of phase 2
            for (vi in sort(very_recent_idx, decreasing = TRUE)) {
              if (amb_ts[vi] > amb_mid * 1.5) {
                amb_phase2_age <- ages[vi]
              } else break
            }
          }
        }
      }
    }

    # Ambrosia at contact window
    near_contact <- test_ages >= 258 & test_ages <= 658
    ambrosia_at_contact <- if (any(near_contact) && length(amb_cols) > 0) {
      any(amb_ts[test_idx][near_contact] > amb_2sd_thresh, na.rm = TRUE)
    } else FALSE

    # Max Ambrosia % and its age
    max_amb_pct <- max(amb_ts, na.rm = TRUE) * 100
    max_amb_age <- if (max_amb_pct > 0) ages[which.max(amb_ts)] else NA

    # ---- ZEA MAYS DETAILED ANALYSIS ----
    zea_cols <- grep("^Zea", taxa, ignore.case = TRUE)
    has_zea <- length(zea_cols) > 0 && any(mat[, zea_cols] > 0)
    zea_first_age <- NA
    zea_max_pct <- 0
    zea_max_age <- NA
    zea_exceeds_baseline <- FALSE
    zea_exceedance_age <- NA

    if (has_zea) {
      zea_ts <- if (length(zea_cols) == 1) mat_prop[, zea_cols] else rowSums(mat_prop[, zea_cols, drop = FALSE])
      zea_raw <- if (length(zea_cols) == 1) mat[, zea_cols] else rowSums(mat[, zea_cols, drop = FALSE])

      zea_present <- zea_raw > 0
      if (any(zea_present)) {
        zea_first_age <- max(ages[zea_present])
        zea_max_pct <- max(zea_ts, na.rm = TRUE) * 100
        zea_max_age <- ages[which.max(zea_ts)]

        # Does Zea exceed baseline?
        zea_baseline <- zea_ts[baseline_idx]
        zea_bl_mean <- mean(zea_baseline, na.rm = TRUE)
        zea_bl_sd <- sd(zea_baseline, na.rm = TRUE)
        if (is.na(zea_bl_sd) || zea_bl_sd == 0) zea_bl_sd <- 0.001
        zea_thresh <- zea_bl_mean + 2 * zea_bl_sd

        zea_test_exc <- zea_ts[test_idx] > zea_thresh
        if (any(zea_test_exc, na.rm = TRUE)) {
          zea_exceeds_baseline <- TRUE
          zea_exceedance_age <- max(test_ages[zea_test_exc], na.rm = TRUE)
        }
      }
    }

    # ---- CHENOPODIUM / OTHER AGRICULTURAL INDICATORS ----
    chen_cols <- grep("^(Chenopodium|Amaranthus)", taxa, ignore.case = TRUE)
    has_chenopod <- length(chen_cols) > 0 && any(mat[, chen_cols] > 0)
    chen_max_pct <- 0
    if (has_chenopod) {
      chen_ts <- if (length(chen_cols) == 1) mat_prop[, chen_cols] else rowSums(mat_prop[, chen_cols, drop = FALSE])
      chen_max_pct <- max(chen_ts, na.rm = TRUE) * 100
    }

    # ---- POST-1492 RECOVERY TEST (Part 5) ----
    # Window: 0-458 BP (post-contact) vs 458-1000 BP (pre-contact)
    post_contact_idx <- which(ages < CONTACT_YEAR & ages >= 0)
    pre_contact_idx  <- which(ages >= CONTACT_YEAR & ages <= 1000)

    recovery_signal <- NA
    bc_precontact_mean <- NA
    bc_postcontact_mean <- NA
    bc_rerise_mean <- NA
    has_rerise <- NA

    if (length(post_contact_idx) >= 2 && length(pre_contact_idx) >= 2) {
      bc_pre  <- mean(bc_from_centroid[pre_contact_idx], na.rm = TRUE)
      bc_post <- mean(bc_from_centroid[post_contact_idx], na.rm = TRUE)
      recovery_signal <- bc_post < bc_pre
      bc_precontact_mean <- round(bc_pre, 4)
      bc_postcontact_mean <- round(bc_post, 4)

      # Re-rise test: does BC increase again after ~200 BP (European re-clearing)?
      rerise_idx <- which(ages <= EUROPEAN_AGRI & ages >= 0)
      mid_post_idx <- which(ages > EUROPEAN_AGRI & ages < CONTACT_YEAR)
      if (length(rerise_idx) >= 1 && length(mid_post_idx) >= 1) {
        bc_rerise <- mean(bc_from_centroid[rerise_idx], na.rm = TRUE)
        bc_mid <- mean(bc_from_centroid[mid_post_idx], na.rm = TRUE)
        bc_rerise_mean <- round(bc_rerise, 4)
        has_rerise <- bc_rerise > bc_mid
      }
    }

    # ---- Classification ----
    human_pct <- sum(group_pct[c("agricultural", "disturbance")], na.rm = TRUE)
    classification <- if (!has_signal) "no_signal"
                      else if (human_pct > 50) "H1_human"
                      else "H3_natural"

    data.frame(
      siteid = as.integer(sid),
      sitename = as.character(meta$sitename[1]),
      lat = as.numeric(meta$lat[1]),
      lon = as.numeric(meta$lon[1]),
      is_zea_target = sid %in% zea_ids,
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
      # Ambrosia detail
      ambrosia_onset_bp = ambrosia_onset,
      ambrosia_explosion_bp = amb_explosion_onset,
      ambrosia_precontact = amb_precontact,
      ambrosia_postcontact = amb_postcontact,
      ambrosia_at_contact = ambrosia_at_contact,
      ambrosia_phase2 = amb_phase2,
      ambrosia_phase2_age = amb_phase2_age,
      max_ambrosia_pct = round(max_amb_pct, 2),
      max_ambrosia_age = max_amb_age,
      # Zea detail
      has_zea = has_zea,
      zea_first_age_bp = zea_first_age,
      zea_max_pct = round(zea_max_pct, 3),
      zea_max_age_bp = zea_max_age,
      zea_exceeds_baseline = zea_exceeds_baseline,
      zea_exceedance_age_bp = zea_exceedance_age,
      # Other agricultural
      has_chenopod = has_chenopod,
      chenopod_max_pct = round(chen_max_pct, 3),
      # Recovery
      recovery_signal = recovery_signal,
      bc_precontact_mean = bc_precontact_mean,
      bc_postcontact_mean = bc_postcontact_mean,
      bc_rerise_mean = bc_rerise_mean,
      has_rerise = has_rerise,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    if(verbose) cat(sprintf("    site %d ERROR: %s\n", sid, e$message))
    NULL
  })
}

# Run analysis
target_sids <- good_sites$siteid
cat(sprintf("  Analyzing %d sites...\n", length(target_sids)))
t0 <- Sys.time()

results_list <- list()
n_fail <- 0
for (i in seq_along(target_sids)) {
  sid <- target_sids[i]
  if (i %% 25 == 0) cat(sprintf("    %d/%d...\n", i, length(target_sids)))
  r <- tryCatch(
    analyze_one_site(sid, verbose = TRUE),
    error = function(e) {
      cat(sprintf("  ERROR site %d: %s\n", sid, e$message))
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

# ============================================================
# Step 4: Ambrosia detailed output (Part 3)
# ============================================================
cat("\n=== Step 4: Ambrosia detailed analysis ===\n")

amb_detail <- results[, c("siteid", "sitename", "lat", "lon",
                           "ambrosia_onset_bp", "ambrosia_explosion_bp",
                           "ambrosia_precontact", "ambrosia_postcontact",
                           "ambrosia_at_contact", "ambrosia_phase2",
                           "ambrosia_phase2_age", "max_ambrosia_pct",
                           "max_ambrosia_age")]

amb_csv_path <- file.path(OUT_DIR, "ena_ambrosia_analysis.csv")
write.csv(amb_detail, amb_csv_path, row.names = FALSE)
cat(sprintf("  Saved Ambrosia detail: %s\n", amb_csv_path))

# ============================================================
# Step 5: Summary statistics and report
# ============================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("FULL ENA EXCEEDANCE ANALYSIS RESULTS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

n_total <- nrow(results)
n_signal <- sum(results$has_signal, na.rm = TRUE)
signal_sites <- results[results$has_signal == TRUE, ]

cat(sprintf("1. DATASET\n"))
cat(sprintf("   Total sites analyzed: %d\n", n_total))
cat(sprintf("   Zea target sites included: %d\n", sum(results$is_zea_target)))
cat(sprintf("   Sites with BC exceedance signal: %d (%.1f%%)\n",
            n_signal, n_signal / n_total * 100))

cat(sprintf("\n2. ONSET TIMING\n"))
if (nrow(signal_sites) > 0) {
  cat(sprintf("   Median onset: %.0f BP\n", median(signal_sites$onset_age_bp, na.rm = TRUE)))
  cat(sprintf("   Range: %d - %d BP\n",
              min(signal_sites$onset_age_bp, na.rm = TRUE),
              max(signal_sites$onset_age_bp, na.rm = TRUE)))

  pre_maize <- sum(signal_sites$onset_age_bp > 2000, na.rm = TRUE)
  maize_era <- sum(signal_sites$onset_age_bp <= 2000 & signal_sites$onset_age_bp > 1000, na.rm = TRUE)
  intensif  <- sum(signal_sites$onset_age_bp <= 1000 & signal_sites$onset_age_bp > CONTACT_YEAR, na.rm = TRUE)
  contact   <- sum(signal_sites$onset_age_bp <= CONTACT_YEAR, na.rm = TRUE)
  cat(sprintf("   Pre-maize (>2000 BP): %d\n", pre_maize))
  cat(sprintf("   Maize era (2000-1000 BP): %d\n", maize_era))
  cat(sprintf("   Intensification (1000-458 BP): %d\n", intensif))
  cat(sprintf("   Post-contact (<458 BP): %d\n", contact))
}

cat(sprintf("\n3. CLASSIFICATION\n"))
cat(sprintf("   H1 (human): %d (%.1f%%)\n",
            sum(results$classification == "H1_human"), sum(results$classification == "H1_human")/n_total*100))
cat(sprintf("   H3 (natural): %d (%.1f%%)\n",
            sum(results$classification == "H3_natural"), sum(results$classification == "H3_natural")/n_total*100))
cat(sprintf("   No signal: %d (%.1f%%)\n",
            sum(results$classification == "no_signal"), sum(results$classification == "no_signal")/n_total*100))

cat(sprintf("\n4. FUNCTIONAL DECOMPOSITION (signal sites)\n"))
if (nrow(signal_sites) > 0) {
  cat(sprintf("   Agricultural: %.1f%%\n", mean(signal_sites$pct_agricultural, na.rm = TRUE)))
  cat(sprintf("   Disturbance:  %.1f%%\n", mean(signal_sites$pct_disturbance, na.rm = TRUE)))
  cat(sprintf("   Deciduous:    %.1f%%\n", mean(signal_sites$pct_deciduous, na.rm = TRUE)))
  cat(sprintf("   Coniferous:   %.1f%%\n", mean(signal_sites$pct_coniferous, na.rm = TRUE)))
  cat(sprintf("   Other:        %.1f%%\n", mean(signal_sites$pct_other, na.rm = TRUE)))
}

cat(sprintf("\n5. AMBROSIA ANALYSIS\n"))
amb_sites <- results[!is.na(results$ambrosia_onset_bp), ]
if (nrow(amb_sites) > 0) {
  cat(sprintf("   Sites with Ambrosia exceedance: %d (%.1f%%)\n", nrow(amb_sites), nrow(amb_sites)/n_total*100))
  cat(sprintf("   Median onset: %.0f BP\n", median(amb_sites$ambrosia_onset_bp, na.rm=TRUE)))

  amb_pre <- sum(amb_sites$ambrosia_onset_bp > CONTACT_YEAR, na.rm = TRUE)
  amb_post <- sum(amb_sites$ambrosia_onset_bp <= CONTACT_YEAR, na.rm = TRUE)
  cat(sprintf("   Pre-Columbian onset: %d\n", amb_pre))
  cat(sprintf("   Post-contact onset: %d\n", amb_post))
  cat(sprintf("   Ambrosia near 1492: %d\n", sum(results$ambrosia_at_contact, na.rm=TRUE)))
  cat(sprintf("   Phase 2 rise (European agri ~200 BP): %d\n", sum(results$ambrosia_phase2, na.rm=TRUE)))

  # Ambrosia explosion classification
  explosion_sites <- results[!is.na(results$ambrosia_explosion_bp), ]
  if (nrow(explosion_sites) > 0) {
    cat(sprintf("   Explosion (>5%% or >baseline+3SD): %d sites\n", nrow(explosion_sites)))
    cat(sprintf("     Pre-contact explosion: %d\n", sum(explosion_sites$ambrosia_precontact)))
    cat(sprintf("     Post-contact explosion: %d\n", sum(explosion_sites$ambrosia_postcontact)))
  }
}

cat(sprintf("\n6. ZEA MAYS\n"))
zea_results <- results[results$has_zea == TRUE, ]
cat(sprintf("   Sites with Zea pollen: %d (%.1f%%)\n", nrow(zea_results), nrow(zea_results)/n_total*100))
if (nrow(zea_results) > 0) {
  cat(sprintf("   Median first appearance: %.0f BP\n", median(zea_results$zea_first_age_bp, na.rm=TRUE)))
  cat(sprintf("   Range: %d - %d BP\n",
              min(zea_results$zea_first_age_bp, na.rm=TRUE), max(zea_results$zea_first_age_bp, na.rm=TRUE)))
  cat(sprintf("   Zea exceeding baseline: %d sites\n", sum(zea_results$zea_exceeds_baseline)))
  cat(sprintf("   Max Zea %% across all sites: %.2f%%\n", max(zea_results$zea_max_pct, na.rm=TRUE)))
}

# Zea target sites report
zea_target_results <- results[results$is_zea_target, ]
cat(sprintf("   Zea TARGET sites analyzed: %d\n", nrow(zea_target_results)))
if (nrow(zea_target_results) > 0) {
  cat(sprintf("   Of these, Zea detected: %d\n", sum(zea_target_results$has_zea)))
  cat(sprintf("   Of these, BC exceedance: %d\n", sum(zea_target_results$has_signal)))
}

cat(sprintf("\n7. OTHER AGRICULTURAL INDICATORS\n"))
cat(sprintf("   Sites with Chenopodium/Amaranthus: %d (%.1f%%)\n",
            sum(results$has_chenopod), sum(results$has_chenopod)/n_total*100))
chen_sites <- results[results$has_chenopod & results$chenopod_max_pct > 1, ]
cat(sprintf("   Sites with Chenopod >1%%: %d\n", nrow(chen_sites)))

cat(sprintf("\n8. POST-1492 RECOVERY\n"))
rec_data <- results[!is.na(results$recovery_signal), ]
if (nrow(rec_data) > 0) {
  n_recovery <- sum(rec_data$recovery_signal, na.rm = TRUE)
  cat(sprintf("   Testable sites: %d\n", nrow(rec_data)))
  cat(sprintf("   Recovery signal (BC decrease after 1492): %d (%.1f%%)\n",
              n_recovery, n_recovery/nrow(rec_data)*100))

  rerise_data <- rec_data[!is.na(rec_data$has_rerise), ]
  if (nrow(rerise_data) > 0) {
    n_rerise <- sum(rerise_data$has_rerise, na.rm = TRUE)
    cat(sprintf("   Sites with European re-clearing rise: %d (%.1f%%)\n",
                n_rerise, n_rerise/nrow(rerise_data)*100))

    # Recovery then rerise = "two agricultures" pattern
    both <- rerise_data[rerise_data$recovery_signal == TRUE & rerise_data$has_rerise == TRUE, ]
    cat(sprintf("   Recovery + re-rise (two agricultures): %d sites\n", nrow(both)))
  }
}

# ============================================================
# Step 6: Two Agricultures Comparison (Part 4)
# ============================================================
cat(sprintf("\n9. TWO AGRICULTURES COMPARISON\n"))

# Pre-Columbian agriculture signal
precolumbian_agri <- sum(results$has_zea | (results$has_chenopod & results$chenopod_max_pct > 0.5), na.rm = TRUE)
cat(sprintf("   Sites with pre-Columbian crop pollen (Zea/Chenopod): %d (%.1f%%)\n",
            precolumbian_agri, precolumbian_agri/n_total*100))

# Post-contact agriculture signal (Ambrosia explosion post-contact)
postcontact_agri <- sum(results$ambrosia_postcontact, na.rm = TRUE)
cat(sprintf("   Sites with post-contact Ambrosia explosion: %d (%.1f%%)\n",
            postcontact_agri, postcontact_agri/n_total*100))

# Natural vegetation signal
nat_signal <- sum(results$has_signal & results$classification == "H3_natural", na.rm = TRUE)
cat(sprintf("   Sites with natural vegetation change: %d (%.1f%%)\n",
            nat_signal, nat_signal/n_total*100))

cat(sprintf("\n   SIGNAL STRENGTH COMPARISON:\n"))
cat(sprintf("   Mean max Ambrosia %%: %.2f%% (all sites)\n", mean(results$max_ambrosia_pct, na.rm=TRUE)))
if (nrow(zea_results) > 0) {
  cat(sprintf("   Mean max Zea %%: %.3f%% (Zea sites only)\n", mean(zea_results$zea_max_pct, na.rm=TRUE)))
}
cat(sprintf("   Ratio Ambrosia/Zea signal: Ambrosia is %dx stronger\n",
            round(mean(results$max_ambrosia_pct, na.rm=TRUE) / max(0.001, mean(zea_results$zea_max_pct, na.rm=TRUE)))))

# ============================================================
# Step 7: Save all outputs
# ============================================================
cat("\n=== Saving outputs ===\n")

# Main CSV
csv_path <- file.path(OUT_DIR, "full_ena_exceedance_data.csv")
write.csv(results, csv_path, row.names = FALSE)
cat(sprintf("  Saved: %s\n", csv_path))

# ---- Markdown Report ----
report <- c(
  "# Full ENA Exceedance Analysis Results",
  "",
  sprintf("**Date**: %s", Sys.Date()),
  sprintf("**Sites analyzed**: %d", n_total),
  sprintf("**Zea target sites**: %d (%d with Zea detected)", sum(results$is_zea_target), sum(zea_target_results$has_zea)),
  sprintf("**Framework**: Bray-Curtis dissimilarity exceedance (Paper 6 method)"),
  sprintf("**Baseline**: >%d BP (pre-maize)", BASELINE_CUTOFF),
  sprintf("**Threshold**: mean + %dSD of baseline BC", EXCEEDANCE_SD),
  "",
  "## 1. Signal Detection",
  "",
  sprintf("- BC exceedance signal: **%d/%d sites (%.1f%%)**", n_signal, n_total, n_signal/n_total*100),
  ""
)

if (nrow(signal_sites) > 0) {
  report <- c(report,
    "## 2. Onset Timing",
    "",
    sprintf("- Median onset: **%.0f BP**", median(signal_sites$onset_age_bp, na.rm=TRUE)),
    "",
    "| Period | N | % of signal |",
    "|--------|---|-------------|",
    sprintf("| Pre-maize (>2000 BP) | %d | %.1f%% |", pre_maize, pre_maize/nrow(signal_sites)*100),
    sprintf("| Maize era (2000-1000 BP) | %d | %.1f%% |", maize_era, maize_era/nrow(signal_sites)*100),
    sprintf("| Intensification (1000-458 BP) | %d | %.1f%% |", intensif, intensif/nrow(signal_sites)*100),
    sprintf("| Post-contact (<458 BP) | %d | %.1f%% |", contact, contact/nrow(signal_sites)*100),
    ""
  )
}

report <- c(report,
  "## 3. Classification",
  "",
  "| Class | N | % |",
  "|-------|---|---|",
  sprintf("| H1 (human-driven) | %d | %.1f%% |",
          sum(results$classification == "H1_human"), sum(results$classification == "H1_human")/n_total*100),
  sprintf("| H3 (natural) | %d | %.1f%% |",
          sum(results$classification == "H3_natural"), sum(results$classification == "H3_natural")/n_total*100),
  sprintf("| No signal | %d | %.1f%% |",
          sum(results$classification == "no_signal"), sum(results$classification == "no_signal")/n_total*100),
  "",
  "## 4. Functional Decomposition (signal sites)",
  "",
  "| Category | Mean % |",
  "|----------|--------|"
)

if (nrow(signal_sites) > 0) {
  report <- c(report,
    sprintf("| Agricultural | %.1f%% |", mean(signal_sites$pct_agricultural, na.rm=TRUE)),
    sprintf("| Disturbance | %.1f%% |", mean(signal_sites$pct_disturbance, na.rm=TRUE)),
    sprintf("| Deciduous | %.1f%% |", mean(signal_sites$pct_deciduous, na.rm=TRUE)),
    sprintf("| Coniferous | %.1f%% |", mean(signal_sites$pct_coniferous, na.rm=TRUE)),
    sprintf("| Other | %.1f%% |", mean(signal_sites$pct_other, na.rm=TRUE))
  )
}

report <- c(report, "",
  "## 5. Ambrosia Analysis",
  ""
)

if (nrow(amb_sites) > 0) {
  report <- c(report,
    sprintf("- Ambrosia exceedance: **%d sites (%.1f%%)**", nrow(amb_sites), nrow(amb_sites)/n_total*100),
    sprintf("- Median onset: **%.0f BP**", median(amb_sites$ambrosia_onset_bp, na.rm=TRUE)),
    sprintf("- Pre-Columbian onset (>458 BP): %d", amb_pre),
    sprintf("- Post-contact onset (<=458 BP): %d", amb_post),
    sprintf("- Phase 2 European agricultural rise: %d sites", sum(results$ambrosia_phase2, na.rm=TRUE)),
    "",
    "### Ambrosia Classification",
    "",
    "| Category | N | % |",
    "|----------|---|---|",
    sprintf("| Pre-contact Ambrosia rise (>458 BP) | %d | %.1f%% |",
            amb_pre, amb_pre/nrow(amb_sites)*100),
    sprintf("| Post-contact only (<=458 BP) | %d | %.1f%% |",
            amb_post, amb_post/nrow(amb_sites)*100),
    sprintf("| Two-phase (initial + European) | %d | %.1f%% |",
            sum(results$ambrosia_phase2, na.rm=TRUE),
            sum(results$ambrosia_phase2, na.rm=TRUE)/nrow(amb_sites)*100),
    ""
  )
}

report <- c(report,
  "## 6. Zea mays Analysis",
  "",
  sprintf("- Sites with Zea pollen: **%d (%.1f%%)**", nrow(zea_results), nrow(zea_results)/n_total*100),
  ""
)

if (nrow(zea_results) > 0) {
  report <- c(report,
    sprintf("- Median first appearance: %.0f BP", median(zea_results$zea_first_age_bp, na.rm=TRUE)),
    sprintf("- Range: %d - %d BP",
            min(zea_results$zea_first_age_bp, na.rm=TRUE), max(zea_results$zea_first_age_bp, na.rm=TRUE)),
    sprintf("- Zea exceeding baseline: %d sites", sum(zea_results$zea_exceeds_baseline)),
    sprintf("- Max Zea %%: %.2f%%", max(zea_results$zea_max_pct, na.rm=TRUE)),
    "",
    "### Zea mays Site Details",
    "",
    "| Site | Lat | Lon | First Zea (BP) | Max Zea % | Exceeds BL | BC Signal |",
    "|------|-----|-----|----------------|-----------|------------|-----------|"
  )
  for (i in 1:nrow(zea_results)) {
    z <- zea_results[i, ]
    report <- c(report,
      sprintf("| %s | %.2f | %.2f | %s | %.3f%% | %s | %s |",
              substr(z$sitename, 1, 25), z$lat, z$lon,
              ifelse(is.na(z$zea_first_age_bp), "-", as.character(z$zea_first_age_bp)),
              z$zea_max_pct,
              ifelse(z$zea_exceeds_baseline, "YES", "no"),
              ifelse(z$has_signal, z$classification, "none")))
  }
  report <- c(report, "")
}

# Zea target sites that were NOT analyzable
zea_not_in <- zea_ids[!zea_ids %in% results$siteid]
if (length(zea_not_in) > 0) {
  report <- c(report,
    sprintf("### Zea target sites not analyzed (%d sites)", length(zea_not_in)),
    "",
    "These sites either failed download, had insufficient samples, or lacked baseline/test periods:",
    sprintf("Site IDs: %s", paste(zea_not_in, collapse = ", ")),
    ""
  )
}

report <- c(report,
  "## 7. Two Agricultures Comparison",
  "",
  "| Signal Type | Metric | Value |",
  "|-------------|--------|-------|",
  sprintf("| Pre-Columbian crops | Sites with Zea/Chenopod | %d (%.1f%%) |",
          precolumbian_agri, precolumbian_agri/n_total*100),
  sprintf("| Pre-Columbian crops | Mean max Zea %% | %.3f%% |",
          ifelse(nrow(zea_results)>0, mean(zea_results$zea_max_pct, na.rm=TRUE), 0)),
  sprintf("| Post-contact weeds | Sites with Ambrosia explosion | %d (%.1f%%) |",
          postcontact_agri, postcontact_agri/n_total*100),
  sprintf("| Post-contact weeds | Mean max Ambrosia %% | %.2f%% |", mean(results$max_ambrosia_pct, na.rm=TRUE)),
  sprintf("| Natural change | BC exceedance (natural) | %d (%.1f%%) |",
          nat_signal, nat_signal/n_total*100),
  "",
  "**Key finding**: Post-contact European agriculture (via Ambrosia) creates a pollen signal",
  sprintf("~%.0fx stronger than pre-Columbian maize agriculture (via Zea pollen).",
          round(mean(results$max_ambrosia_pct, na.rm=TRUE) / max(0.001, ifelse(nrow(zea_results)>0, mean(zea_results$zea_max_pct, na.rm=TRUE), 0.001)))),
  ""
)

report <- c(report,
  "## 8. Post-1492 Recovery Test",
  ""
)

if (nrow(rec_data) > 0) {
  n_recovery <- sum(rec_data$recovery_signal, na.rm = TRUE)
  report <- c(report,
    sprintf("- Testable sites: %d", nrow(rec_data)),
    sprintf("- Recovery (BC decrease after contact): **%d (%.1f%%)**",
            n_recovery, n_recovery/nrow(rec_data)*100)
  )

  rerise_data <- rec_data[!is.na(rec_data$has_rerise), ]
  if (nrow(rerise_data) > 0) {
    n_rerise <- sum(rerise_data$has_rerise, na.rm = TRUE)
    both <- rerise_data[rerise_data$recovery_signal == TRUE & rerise_data$has_rerise == TRUE, ]
    report <- c(report,
      sprintf("- European re-clearing rise: %d (%.1f%%)", n_rerise, n_rerise/nrow(rerise_data)*100),
      sprintf("- **Recovery then re-rise (two agricultures pattern): %d sites**", nrow(both)),
      ""
    )

    # For sites with pre-contact Ambrosia exceedance
    precontact_amb <- results[results$ambrosia_precontact == TRUE & !is.na(results$recovery_signal), ]
    if (nrow(precontact_amb) > 0) {
      n_rec_precontact <- sum(precontact_amb$recovery_signal, na.rm = TRUE)
      report <- c(report,
        "### Recovery in pre-contact Ambrosia sites",
        "",
        sprintf("- Pre-contact Ambrosia sites with recovery data: %d", nrow(precontact_amb)),
        sprintf("- Of these, showing recovery: %d (%.1f%%)",
                n_rec_precontact, n_rec_precontact/nrow(precontact_amb)*100),
        ""
      )
    }
  }
}

# Site-level table
report <- c(report,
  "## Site-Level Results (all sites)",
  "",
  "| Site | Lat | Lon | Signal | Onset BP | Class | Amb onset | Amb phase2 | Zea | Recovery | Re-rise |",
  "|------|-----|-----|--------|----------|-------|-----------|------------|-----|----------|---------|"
)

for (i in 1:nrow(results)) {
  r <- results[i, ]
  report <- c(report,
    sprintf("| %s | %.2f | %.2f | %s | %s | %s | %s | %s | %s | %s | %s |",
            substr(r$sitename, 1, 22),
            r$lat, r$lon,
            ifelse(r$has_signal, "YES", "no"),
            ifelse(is.na(r$onset_age_bp), "-", as.character(round(r$onset_age_bp))),
            r$classification,
            ifelse(is.na(r$ambrosia_onset_bp), "-", as.character(round(r$ambrosia_onset_bp))),
            ifelse(r$ambrosia_phase2, "YES", "-"),
            ifelse(r$has_zea, "YES", "-"),
            ifelse(is.na(r$recovery_signal), "?", ifelse(r$recovery_signal, "YES", "no")),
            ifelse(is.na(r$has_rerise), "?", ifelse(r$has_rerise, "YES", "no"))))
}

report_path <- file.path(OUT_DIR, "full_ena_exceedance_results.md")
writeLines(report, report_path)
cat(sprintf("  Saved report: %s\n", report_path))

cat("\n=== DONE ===\n")
