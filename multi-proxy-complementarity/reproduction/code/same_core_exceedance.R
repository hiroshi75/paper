#!/usr/bin/env Rscript
# =============================================================================
# Same-Core Exceedance Analysis
# =============================================================================
# For sites with BOTH pollen and charcoal at the same siteid,
# run exceedance analysis on both proxies and build 2x2 table.
# Uses neotoma2 R package for reliable data extraction.
# =============================================================================

library(neotoma2)
library(dplyr)
library(tidyr)

output_lines <- c()
log <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  output_lines <<- c(output_lines, msg)
}

# Load the site summary from Step 1
site_summary <- readRDS("/home/ayu/archeco/shared/cache/same_site_pollen_charcoal.rds")
cat("Total same-site pairs:", nrow(site_summary), "\n")

log("# Same-Core Pollen + Charcoal Exceedance Analysis")
log("")
log(paste0("Total same-site pairs in Neotoma: ", nrow(site_summary)))
log("")

# =============================================================================
# Download ALL same-site pairs via neotoma2 (batch by site ID)
# =============================================================================
log("## Phase 1: Download and assess temporal coverage")
log("")

# Process all sites (not just first 50)
all_results <- list()
all_samples <- list()

for (i in seq_len(nrow(site_summary))) {
  site <- site_summary[i, ]
  sid <- site$siteid

  cache_file <- sprintf("/home/ayu/archeco/shared/cache/neotoma2_site_%d.rds", sid)

  cat(sprintf("  [%d/%d] Site %d: %s", i, nrow(site_summary), sid, site$sitename))

  tryCatch({
    samp <- NULL

    if (file.exists(cache_file)) {
      samp <- readRDS(cache_file)
      cat(" (cached)")
    } else {
      s <- get_sites(sid)
      ds <- get_datasets(s)
      dl <- get_downloads(ds)
      samp <- samples(dl)
      saveRDS(samp, cache_file)
    }

    if (is.null(samp) || nrow(samp) == 0) {
      cat(" - no data\n")
      next
    }

    # Split by type
    char_samp <- samp %>% filter(datasettype == "charcoal")
    poll_samp <- samp %>% filter(grepl("^pollen$", datasettype))

    if (nrow(char_samp) == 0 || nrow(poll_samp) == 0) {
      cat(" - missing one type\n")
      next
    }

    # Age ranges
    char_ages <- range(char_samp$age, na.rm = TRUE)
    poll_ages <- range(poll_samp$age, na.rm = TRUE)

    overlap_min <- max(char_ages[1], poll_ages[1])
    overlap_max <- min(char_ages[2], poll_ages[2])
    has_overlap <- overlap_min < overlap_max
    has_baseline <- has_overlap && overlap_max >= 5000

    all_results[[length(all_results) + 1]] <- data.frame(
      siteid = sid, sitename = site$sitename,
      lat = site$lat, lon = site$lon, region = site$region,
      char_young = char_ages[1], char_old = char_ages[2],
      char_n = nrow(char_samp),
      poll_young = poll_ages[1], poll_old = poll_ages[2],
      poll_n = nrow(poll_samp),
      has_overlap = has_overlap, has_baseline_5k = has_baseline,
      stringsAsFactors = FALSE
    )

    all_samples[[as.character(sid)]] <- list(charcoal = char_samp, pollen = poll_samp)

    cat(sprintf(" - Char:%d-%d(%d), Poll:%d-%d(%d), BL>5k:%s\n",
                round(char_ages[1]), round(char_ages[2]), nrow(char_samp),
                round(poll_ages[1]), round(poll_ages[2]), nrow(poll_samp),
                if(has_baseline) "YES" else "NO"))

  }, error = function(e) {
    cat(sprintf(" - ERROR: %s\n", e$message))
  })

  Sys.sleep(0.5)
}

results_df <- bind_rows(all_results)
saveRDS(results_df, "/home/ayu/archeco/shared/cache/same_site_temporal.rds")

log(paste0("Successfully processed: ", nrow(results_df), " sites"))
log(paste0("With temporal overlap: ", sum(results_df$has_overlap, na.rm = TRUE)))
log(paste0("With baseline >5000 BP: ", sum(results_df$has_baseline_5k, na.rm = TRUE)))

log("")
log("### Sites with baseline >5000 BP (feasible for exceedance)")
log("")
log("| Site | Region | Lat | Lon | Char Range | Poll Range | Char N | Poll N |")
log("|------|--------|-----|-----|-----------|-----------|--------|--------|")
feasible <- results_df %>% filter(has_baseline_5k == TRUE)
for (i in seq_len(nrow(feasible))) {
  r <- feasible[i, ]
  log(sprintf("| %s | %s | %.1f | %.1f | %d-%d | %d-%d | %d | %d |",
              substr(r$sitename, 1, 30), r$region, r$lat, r$lon,
              round(r$char_young), round(r$char_old),
              round(r$poll_young), round(r$poll_old),
              r$char_n, r$poll_n))
}

# =============================================================================
# Phase 2: Exceedance Analysis on feasible sites
# =============================================================================
log("")
log("## Phase 2: Exceedance Analysis")
log("")

exc_results <- list()

for (i in seq_len(nrow(feasible))) {
  site <- feasible[i, ]
  sid <- as.character(site$siteid)

  if (!sid %in% names(all_samples)) next

  cat(sprintf("\n  === Exceedance [%d/%d]: %s (site %s) ===\n",
              i, nrow(feasible), site$sitename, sid))

  data <- all_samples[[sid]]
  char_samp <- data$charcoal
  poll_samp <- data$pollen

  # --- CHARCOAL EXCEEDANCE ---
  # Aggregate charcoal by age (sum all size fractions)
  char_agg <- char_samp %>%
    filter(!is.na(age) & !is.na(value)) %>%
    group_by(age) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(age)

  char_exceed_pct <- NA
  char_signal <- "INSUFFICIENT"
  char_baseline_n <- 0
  char_test_n <- 0

  baseline_char <- char_agg %>% filter(age >= 5000)
  test_char <- char_agg %>% filter(age < 5000 & age >= 0)

  if (nrow(baseline_char) >= 5 && nrow(test_char) >= 3) {
    bl_mean <- mean(baseline_char$value, na.rm = TRUE)
    bl_sd <- sd(baseline_char$value, na.rm = TRUE)
    threshold <- bl_mean + 2 * bl_sd

    n_exceed <- sum(test_char$value > threshold, na.rm = TRUE)
    char_exceed_pct <- round(n_exceed / nrow(test_char) * 100, 1)
    char_signal <- if (char_exceed_pct > 10) "EXCEEDANCE" else "NO_SIGNAL"
    char_baseline_n <- nrow(baseline_char)
    char_test_n <- nrow(test_char)

    cat(sprintf("    Charcoal: baseline=%d samples (mean=%.1f, sd=%.1f, threshold=%.1f)\n",
                nrow(baseline_char), bl_mean, bl_sd, threshold))
    cat(sprintf("    Test period=%d samples, exceed=%d (%.1f%%)\n",
                nrow(test_char), n_exceed, char_exceed_pct))
  } else {
    cat(sprintf("    Charcoal: INSUFFICIENT (baseline=%d, test=%d)\n",
                nrow(baseline_char), nrow(test_char)))
  }

  # --- POLLEN EXCEEDANCE (AP ratio decrease) ---
  poll_exceed_pct <- NA
  poll_signal <- "INSUFFICIENT"
  poll_indicator <- "N/A"
  poll_baseline_n <- 0
  poll_test_n <- 0

  # Compute AP ratio per sample age
  ap_data <- poll_samp %>%
    filter(!is.na(age) & !is.na(value) & ecologicalgroup %in% c("TRSH", "UPHE")) %>%
    group_by(age) %>%
    summarise(
      ap = sum(value[ecologicalgroup == "TRSH"], na.rm = TRUE),
      nap = sum(value[ecologicalgroup == "UPHE"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter((ap + nap) > 0) %>%
    mutate(ap_ratio = ap / (ap + nap)) %>%
    arrange(age)

  baseline_ap <- ap_data %>% filter(age >= 5000)
  test_ap <- ap_data %>% filter(age < 5000 & age >= 0)

  if (nrow(baseline_ap) >= 5 && nrow(test_ap) >= 3) {
    bl_mean_ap <- mean(baseline_ap$ap_ratio, na.rm = TRUE)
    bl_sd_ap <- sd(baseline_ap$ap_ratio, na.rm = TRUE)
    # AP decrease = deforestation
    threshold_ap <- bl_mean_ap - 2 * bl_sd_ap

    n_exceed_ap <- sum(test_ap$ap_ratio < threshold_ap, na.rm = TRUE)
    poll_exceed_pct <- round(n_exceed_ap / nrow(test_ap) * 100, 1)
    poll_signal <- if (poll_exceed_pct > 10) "EXCEEDANCE" else "NO_SIGNAL"
    poll_indicator <- sprintf("AP_decrease (bl=%.2f, thr=%.2f)", bl_mean_ap, threshold_ap)
    poll_baseline_n <- nrow(baseline_ap)
    poll_test_n <- nrow(test_ap)

    cat(sprintf("    Pollen AP: baseline=%d (mean=%.3f, sd=%.3f, threshold=%.3f)\n",
                nrow(baseline_ap), bl_mean_ap, bl_sd_ap, threshold_ap))
    cat(sprintf("    Test period=%d, exceed=%d (%.1f%%)\n",
                nrow(test_ap), n_exceed_ap, poll_exceed_pct))
  } else {
    cat(sprintf("    Pollen AP: INSUFFICIENT (baseline=%d, test=%d)\n",
                nrow(baseline_ap), nrow(test_ap)))
  }

  # Check for crop taxa
  crop_taxa <- c("Cerealia", "Cerealia-type", "Secale", "Triticum", "Hordeum",
                 "Zea", "Zea mays", "Plantago lanceolata", "Plantago",
                 "Cannabis", "Oryza", "Rumex")
  found_crops <- unique(poll_samp$variablename[poll_samp$variablename %in% crop_taxa])
  if (length(found_crops) > 0) {
    poll_indicator <- paste(poll_indicator, "+", paste(found_crops, collapse="/"))
    cat(sprintf("    Crop taxa found: %s\n", paste(found_crops, collapse=", ")))
  }

  exc_results[[length(exc_results) + 1]] <- data.frame(
    siteid = as.integer(sid),
    sitename = site$sitename,
    lat = site$lat, lon = site$lon, region = site$region,
    char_baseline_n = char_baseline_n, char_test_n = char_test_n,
    char_exceed_pct = char_exceed_pct,
    char_signal = char_signal,
    poll_baseline_n = poll_baseline_n, poll_test_n = poll_test_n,
    poll_exceed_pct = poll_exceed_pct,
    poll_signal = poll_signal,
    poll_indicator = poll_indicator,
    stringsAsFactors = FALSE
  )
}

if (length(exc_results) > 0) {
  exc_df <- bind_rows(exc_results)
  saveRDS(exc_df, "/home/ayu/archeco/shared/cache/same_site_exceedance.rds")
  write.csv(exc_df, "/home/ayu/archeco/shared/same_site_exceedance.csv", row.names = FALSE)

  log(paste0("Exceedance analysis completed for: ", nrow(exc_df), " sites"))
  log("")
  log("### Results")
  log("")
  log("| Site | Region | Char BL/Test | Char Exceed% | Char Signal | Poll BL/Test | Poll Exceed% | Poll Signal |")
  log("|------|--------|-------------|-------------|-------------|-------------|-------------|-------------|")
  for (i in seq_len(nrow(exc_df))) {
    r <- exc_df[i, ]
    log(sprintf("| %s | %s | %d/%d | %s | %s | %d/%d | %s | %s |",
                substr(r$sitename, 1, 25), r$region,
                r$char_baseline_n, r$char_test_n,
                if(!is.na(r$char_exceed_pct)) paste0(r$char_exceed_pct, "%") else "N/A",
                r$char_signal,
                r$poll_baseline_n, r$poll_test_n,
                if(!is.na(r$poll_exceed_pct)) paste0(r$poll_exceed_pct, "%") else "N/A",
                r$poll_signal))
  }

  # =============================================================================
  # 2x2 Contingency Table
  # =============================================================================
  log("")
  log("## 2x2 Contingency Table")
  log("")

  valid <- exc_df %>% filter(char_signal != "INSUFFICIENT" & poll_signal != "INSUFFICIENT")
  log(paste0("Sites with both valid signals: ", nrow(valid)))

  if (nrow(valid) >= 2) {
    tab <- table(
      Charcoal = valid$char_signal,
      Pollen = valid$poll_signal
    )
    log("")
    log("```")
    log(paste(capture.output(print(tab)), collapse = "\n"))
    log("```")

    if (all(dim(tab) >= 2) && nrow(valid) >= 4) {
      ft <- fisher.test(tab)
      log("")
      log(paste0("Fisher's exact test p = ", format(ft$p.value, digits = 4)))
      log(paste0("Odds ratio = ", format(ft$estimate, digits = 3)))
    }

    # Concordance
    n_concordant <- sum(valid$char_signal == valid$poll_signal)
    log("")
    log(paste0("Concordance rate: ", n_concordant, "/", nrow(valid),
               " (", round(n_concordant / nrow(valid) * 100, 1), "%)"))

    # Interpretation
    log("")
    log("### Interpretation")
    log("")
    if (nrow(valid) >= 5) {
      ee <- sum(valid$char_signal == "EXCEEDANCE" & valid$poll_signal == "EXCEEDANCE")
      en <- sum(valid$char_signal == "EXCEEDANCE" & valid$poll_signal == "NO_SIGNAL")
      ne <- sum(valid$char_signal == "NO_SIGNAL" & valid$poll_signal == "EXCEEDANCE")
      nn <- sum(valid$char_signal == "NO_SIGNAL" & valid$poll_signal == "NO_SIGNAL")

      log(sprintf("- Both exceedance (fire + deforestation): %d sites", ee))
      log(sprintf("- Charcoal only (fire management, no deforestation): %d sites", en))
      log(sprintf("- Pollen only (deforestation without fire): %d sites", ne))
      log(sprintf("- Neither (natural baseline maintained): %d sites", nn))
    }
  } else {
    log("Insufficient valid sites for contingency table")
  }
}

# =============================================================================
# Save report
# =============================================================================
log("")
log("---")
log(paste0("Generated: ", Sys.time()))

# Update main report
main_report <- readLines("/home/ayu/archeco/shared/same_core_pollen_charcoal.md")
full_report <- c(main_report, "", "---", "", output_lines)
writeLines(full_report, "/home/ayu/archeco/shared/same_core_pollen_charcoal.md")

cat("\n=== Report appended to shared/same_core_pollen_charcoal.md ===\n")
