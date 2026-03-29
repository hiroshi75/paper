#!/usr/bin/env Rscript
# =============================================================================
# Find Same-Site Pollen + Charcoal Pairs in Neotoma
# =============================================================================
# Sites where BOTH pollen and charcoal datasets exist at the SAME site
# (same siteid = same core/lake/bog) — ideal for cross-proxy comparison.
# =============================================================================

library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)

output_lines <- c()
log <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  output_lines <<- c(output_lines, msg)
}

# Helper: extract coordinates from GeoJSON geography string
parse_coords <- function(geo_str) {
  tryCatch({
    geo <- fromJSON(geo_str)
    if (geo$type == "Point") {
      return(c(lon = geo$coordinates[1], lat = geo$coordinates[2]))
    }
    # Polygon → centroid
    coords <- geo$coordinates[[1]]
    return(c(lon = mean(coords[, 1]), lat = mean(coords[, 2])))
  }, error = function(e) c(lon = NA, lat = NA))
}

# =============================================================================
# STEP 1: Query ALL charcoal datasets (smaller set ~400)
# =============================================================================
log("# Same-Site Pollen + Charcoal Analysis")
log("")
log("## Step 1: Query Neotoma for charcoal datasets")
log("")

query_datasets <- function(datasettype, limit = 5000, offset = 0) {
  url <- "https://api.neotomadb.org/v2.0/data/datasets"
  all_items <- list()

  repeat {
    cat("  Querying", datasettype, "offset=", offset, "...\n")
    resp <- tryCatch(
      GET(url, query = list(datasettype = datasettype, limit = limit, offset = offset),
          timeout(120)),
      error = function(e) { cat("  Error:", e$message, "\n"); NULL }
    )
    if (is.null(resp) || status_code(resp) != 200) break

    data <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
    items <- data$data
    if (is.null(items) || !is.data.frame(items) || nrow(items) == 0) break

    all_items[[length(all_items) + 1]] <- items
    cat("  Got", nrow(items), "records\n")
    if (nrow(items) < limit) break
    offset <- offset + limit
  }

  if (length(all_items) > 0) bind_rows(all_items) else data.frame()
}

char_raw <- query_datasets("charcoal")
log(paste0("Charcoal dataset records: ", nrow(char_raw)))

# Extract site info with coordinates
extract_sites <- function(raw) {
  rows <- list()
  for (i in seq_len(nrow(raw))) {
    sid <- raw$site.siteid[i]
    sname <- raw$site.sitename[i]
    geo <- raw$site.geography[i]
    coords <- parse_coords(geo)

    # Extract dataset IDs from nested list
    ds_list <- raw$site.datasets[[i]]
    if (is.data.frame(ds_list) && nrow(ds_list) > 0) {
      for (j in seq_len(nrow(ds_list))) {
        rows[[length(rows) + 1]] <- data.frame(
          siteid = sid, sitename = sname,
          lat = coords["lat"], lon = coords["lon"],
          datasetid = ds_list$datasetid[j],
          datasettype = ds_list$datasettype[j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  bind_rows(rows)
}

char_sites <- extract_sites(char_raw)
char_sites <- char_sites %>% filter(datasettype == "charcoal")
cat("Charcoal datasets extracted:", nrow(char_sites), "\n")
cat("Unique charcoal sites:", length(unique(char_sites$siteid)), "\n")

# =============================================================================
# STEP 2: Query pollen datasets (larger, need pagination)
# =============================================================================
log("")
log("## Step 2: Query pollen datasets and find overlap")
log("")

pollen_raw <- query_datasets("pollen")
pollen_sites <- extract_sites(pollen_raw)
pollen_sites <- pollen_sites %>% filter(grepl("pollen", datasettype, ignore.case = TRUE))
cat("Pollen datasets extracted:", nrow(pollen_sites), "\n")
cat("Unique pollen sites:", length(unique(pollen_sites$siteid)), "\n")

# Need to get more pollen sites - the API truncated at 5000 dataset records
# but that covers many unique sites. Also query with offset
if (nrow(pollen_raw) >= 5000) {
  cat("Pollen truncated at 5000, fetching more...\n")
  pollen_raw2 <- query_datasets("pollen", limit = 5000, offset = 5000)
  if (nrow(pollen_raw2) > 0) {
    pollen_sites2 <- extract_sites(pollen_raw2)
    pollen_sites2 <- pollen_sites2 %>% filter(grepl("pollen", datasettype, ignore.case = TRUE))
    pollen_sites <- bind_rows(pollen_sites, pollen_sites2)
    cat("After 2nd batch - pollen datasets:", nrow(pollen_sites), "\n")
    cat("Unique pollen sites:", length(unique(pollen_sites$siteid)), "\n")
  }
}

# =============================================================================
# STEP 3: Find intersection
# =============================================================================
log("")
log("## Step 3: Same-site pairs")
log("")

shared_ids <- intersect(unique(char_sites$siteid), unique(pollen_sites$siteid))
log(paste0("Unique charcoal sites: ", length(unique(char_sites$siteid))))
log(paste0("Unique pollen sites: ", length(unique(pollen_sites$siteid))))
log(paste0("**Same-site pairs: ", length(shared_ids), "**"))

# Build summary table
site_summary <- char_sites %>%
  filter(siteid %in% shared_ids) %>%
  group_by(siteid, sitename, lat, lon) %>%
  summarise(n_char = n(), char_dsids = paste(datasetid, collapse = ","), .groups = "drop") %>%
  left_join(
    pollen_sites %>%
      filter(siteid %in% shared_ids) %>%
      group_by(siteid) %>%
      summarise(n_poll = n(), poll_dsids = paste(datasetid, collapse = ","), .groups = "drop"),
    by = "siteid"
  )

# Assign regions
site_summary <- site_summary %>%
  mutate(region = case_when(
    lon >= -170 & lon < -30 & lat >= 0 ~ "North America",
    lon >= -170 & lon < -30 & lat < 0 ~ "South/Central America",
    lon >= -30 & lon < 45 & lat >= 35 ~ "Europe",
    lon >= -30 & lon < 45 & lat >= 0 & lat < 35 ~ "Africa/Mediterranean",
    lon >= -30 & lon < 45 & lat < 0 ~ "Africa",
    lon >= 45 ~ "Asia/Oceania",
    TRUE ~ "Other"
  ))

log("")
log("### Geographic Distribution")
log("")
region_tab <- site_summary %>% count(region) %>% arrange(desc(n))
log("| Region | Count |")
log("|--------|-------|")
for (i in seq_len(nrow(region_tab))) {
  log(sprintf("| %s | %d |", region_tab$region[i], region_tab$n[i]))
}

# Full site table
log("")
log("### All Same-Site Pairs")
log("")
log("| Site ID | Site Name | Lat | Lon | Region | Char DS | Pollen DS |")
log("|---------|-----------|-----|-----|--------|---------|-----------|")
for (i in seq_len(nrow(site_summary))) {
  r <- site_summary[i, ]
  log(sprintf("| %d | %s | %.2f | %.2f | %s | %d | %d |",
              r$siteid, substr(r$sitename, 1, 35),
              r$lat, r$lon, r$region, r$n_char, r$n_poll))
}

saveRDS(site_summary, "/home/ayu/archeco/shared/cache/same_site_pollen_charcoal.rds")
write.csv(site_summary, "/home/ayu/archeco/shared/same_site_pollen_charcoal.csv", row.names = FALSE)

# =============================================================================
# STEP 4: Download data and check temporal overlap (first 50)
# =============================================================================
log("")
log("## Step 4: Temporal Overlap Assessment (first 50 sites)")
log("")

download_limit <- min(50, nrow(site_summary))
detailed_results <- list()

for (i in seq_len(download_limit)) {
  site <- site_summary[i, ]
  sid <- site$siteid
  cat(sprintf("  [%d/%d] Site %d: %s\n", i, download_limit, sid, site$sitename))

  # Get first charcoal and pollen dataset ID
  char_dsid <- as.integer(strsplit(site$char_dsids, ",")[[1]][1])
  poll_dsid <- as.integer(strsplit(site$poll_dsids, ",")[[1]][1])

  get_samples <- function(dsid) {
    cache_f <- sprintf("/home/ayu/archeco/shared/cache/ds_%d.rds", dsid)
    if (file.exists(cache_f)) return(readRDS(cache_f))

    url <- paste0("https://api.neotomadb.org/v2.0/data/downloads/", dsid)
    resp <- tryCatch(GET(url, timeout(60)), error = function(e) NULL)
    if (is.null(resp) || status_code(resp) != 200) return(NULL)

    data <- tryCatch(
      fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE),
      error = function(e) NULL
    )
    if (!is.null(data)) saveRDS(data$data, cache_f)
    data$data
  }

  char_data <- get_samples(char_dsid)
  poll_data <- get_samples(poll_dsid)

  # Extract age range from download data
  extract_age_range <- function(dl) {
    tryCatch({
      # Age range may be in site.collectionunit.dataset.agerange
      ar <- NULL
      cols <- names(dl)
      ar_col <- grep("agerange", cols, value = TRUE)
      if (length(ar_col) > 0) {
        ar_data <- dl[[ar_col[1]]]
        if (is.list(ar_data) && length(ar_data) > 0) {
          ar_df <- ar_data[[1]]
          if (is.data.frame(ar_df) && nrow(ar_df) > 0) {
            return(list(
              min_age = ar_df$ageold[1],
              max_age = ar_df$ageyoung[1],
              n_samples = NA
            ))
          }
        }
      }

      # Try to count samples
      samp_col <- grep("samples", cols, value = TRUE)
      if (length(samp_col) > 0) {
        samp <- dl[[samp_col[1]]]
        if (is.list(samp) && length(samp) > 0) {
          samp_df <- samp[[1]]
          if (is.data.frame(samp_df) && "age" %in% names(samp_df)) {
            ages <- samp_df$age
            return(list(
              min_age = min(ages, na.rm = TRUE),
              max_age = max(ages, na.rm = TRUE),
              n_samples = length(unique(ages))
            ))
          }
        }
      }
      NULL
    }, error = function(e) NULL)
  }

  char_ages <- extract_age_range(char_data)
  poll_ages <- extract_age_range(poll_data)

  # Calculate overlap
  has_overlap <- FALSE
  has_baseline <- FALSE
  overlap_min <- NA
  overlap_max <- NA

  if (!is.null(char_ages) && !is.null(poll_ages)) {
    c_young <- min(char_ages$min_age, char_ages$max_age, na.rm = TRUE)
    c_old <- max(char_ages$min_age, char_ages$max_age, na.rm = TRUE)
    p_young <- min(poll_ages$min_age, poll_ages$max_age, na.rm = TRUE)
    p_old <- max(poll_ages$min_age, poll_ages$max_age, na.rm = TRUE)

    overlap_min <- max(c_young, p_young)
    overlap_max <- min(c_old, p_old)
    has_overlap <- !is.na(overlap_min) && !is.na(overlap_max) && overlap_min < overlap_max
    has_baseline <- has_overlap && overlap_max >= 5000
  }

  detailed_results[[length(detailed_results) + 1]] <- data.frame(
    siteid = sid, sitename = site$sitename,
    lat = site$lat, lon = site$lon, region = site$region,
    char_dsid = char_dsid, poll_dsid = poll_dsid,
    char_young = if (!is.null(char_ages)) min(char_ages$min_age, char_ages$max_age, na.rm=T) else NA,
    char_old = if (!is.null(char_ages)) max(char_ages$min_age, char_ages$max_age, na.rm=T) else NA,
    poll_young = if (!is.null(poll_ages)) min(poll_ages$min_age, poll_ages$max_age, na.rm=T) else NA,
    poll_old = if (!is.null(poll_ages)) max(poll_ages$min_age, poll_ages$max_age, na.rm=T) else NA,
    has_overlap = has_overlap,
    overlap_min = overlap_min, overlap_max = overlap_max,
    has_baseline_5k = has_baseline,
    stringsAsFactors = FALSE
  )

  cat(sprintf("    Char: %s, Pollen: %s, Overlap: %s, Baseline>5k: %s\n",
              if(!is.null(char_ages)) paste(char_ages$min_age, "-", char_ages$max_age) else "N/A",
              if(!is.null(poll_ages)) paste(poll_ages$min_age, "-", poll_ages$max_age) else "N/A",
              if(has_overlap) "YES" else "NO",
              if(has_baseline) "YES" else "NO"))

  Sys.sleep(0.3)
}

detailed_df <- bind_rows(detailed_results)
saveRDS(detailed_df, "/home/ayu/archeco/shared/cache/same_site_detailed.rds")
write.csv(detailed_df, "/home/ayu/archeco/shared/same_site_detailed.csv", row.names = FALSE)

log(paste0("Downloaded metadata for: ", nrow(detailed_df), " sites"))
log(paste0("With temporal overlap: ", sum(detailed_df$has_overlap, na.rm = TRUE)))
log(paste0("With baseline >5000 BP: ", sum(detailed_df$has_baseline_5k, na.rm = TRUE)))

log("")
log("### Site Temporal Details")
log("")
log("| Site | Region | Char Ages | Pollen Ages | Overlap | Baseline>5k |")
log("|------|--------|-----------|-------------|---------|-------------|")
for (i in seq_len(nrow(detailed_df))) {
  r <- detailed_df[i, ]
  log(sprintf("| %s | %s | %s-%s | %s-%s | %s | %s |",
              substr(r$sitename, 1, 25), r$region,
              if(!is.na(r$char_young)) round(r$char_young) else "?",
              if(!is.na(r$char_old)) round(r$char_old) else "?",
              if(!is.na(r$poll_young)) round(r$poll_young) else "?",
              if(!is.na(r$poll_old)) round(r$poll_old) else "?",
              if(!is.na(r$has_overlap) && r$has_overlap) "YES" else "NO",
              if(!is.na(r$has_baseline_5k) && r$has_baseline_5k) "YES" else "NO"))
}

# =============================================================================
# STEP 5: Exceedance analysis on feasible sites
# =============================================================================
log("")
log("## Step 5: Exceedance Analysis")
log("")

feasible <- detailed_df %>% filter(has_baseline_5k == TRUE)
log(paste0("Feasible sites for exceedance: ", nrow(feasible)))

if (nrow(feasible) > 0) {
  # For each feasible site, download full sample data
  exc_results <- list()

  for (i in seq_len(nrow(feasible))) {
    site <- feasible[i, ]
    sid <- site$siteid
    cat(sprintf("  Exceedance [%d/%d]: site %d (%s)\n", i, nrow(feasible), sid, site$sitename))

    # Download full samples via API
    get_full_samples <- function(dsid) {
      url <- paste0("https://api.neotomadb.org/v2.0/data/downloads/", dsid)
      resp <- tryCatch(GET(url, timeout(60)), error = function(e) NULL)
      if (is.null(resp) || status_code(resp) != 200) return(NULL)

      data <- tryCatch(
        fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE),
        error = function(e) NULL
      )
      if (is.null(data)) return(NULL)

      # Extract samples
      tryCatch({
        dl <- data$data
        samp_col <- grep("samples$", names(dl), value = TRUE)
        if (length(samp_col) == 0) return(NULL)
        samp <- dl[[samp_col[length(samp_col)]]]  # Take the last (deepest) samples
        if (is.list(samp) && length(samp) > 0) {
          samp_df <- samp[[1]]
          if (is.data.frame(samp_df)) return(samp_df)
        }
        NULL
      }, error = function(e) NULL)
    }

    char_samp <- get_full_samples(site$char_dsid)
    poll_samp <- get_full_samples(site$poll_dsid)

    # --- Charcoal exceedance ---
    char_exceed_pct <- NA
    char_signal <- "INSUFFICIENT"

    if (!is.null(char_samp) && "age" %in% names(char_samp) && "value" %in% names(char_samp)) {
      # Aggregate by age (sum all charcoal variables)
      char_agg <- char_samp %>%
        group_by(age) %>%
        summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

      baseline <- char_agg %>% filter(age >= 5000)
      test_period <- char_agg %>% filter(age < 5000 & age >= 0)

      if (nrow(baseline) >= 5 && nrow(test_period) >= 3) {
        bl_mean <- mean(baseline$value, na.rm = TRUE)
        bl_sd <- sd(baseline$value, na.rm = TRUE)
        threshold <- bl_mean + 2 * bl_sd

        char_exceed_pct <- sum(test_period$value > threshold, na.rm = TRUE) / nrow(test_period) * 100
        char_signal <- if (char_exceed_pct > 10) "EXCEEDANCE" else "NO_SIGNAL"
      }
    }

    # --- Pollen exceedance (AP ratio decrease) ---
    poll_exceed_pct <- NA
    poll_signal <- "INSUFFICIENT"
    poll_indicator <- "N/A"

    if (!is.null(poll_samp) && "age" %in% names(poll_samp)) {
      # Check for ecological group info
      if ("ecologicalgroup" %in% names(poll_samp) && "value" %in% names(poll_samp)) {
        ap_data <- poll_samp %>%
          filter(ecologicalgroup %in% c("TRSH", "UPHE")) %>%
          group_by(age) %>%
          summarise(
            ap = sum(value[ecologicalgroup == "TRSH"], na.rm = TRUE),
            nap = sum(value[ecologicalgroup == "UPHE"], na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(ap_ratio = ap / (ap + nap))

        baseline_ap <- ap_data %>% filter(age >= 5000)
        test_ap <- ap_data %>% filter(age < 5000 & age >= 0)

        if (nrow(baseline_ap) >= 5 && nrow(test_ap) >= 3) {
          bl_mean <- mean(baseline_ap$ap_ratio, na.rm = TRUE)
          bl_sd <- sd(baseline_ap$ap_ratio, na.rm = TRUE)
          threshold <- bl_mean - 2 * bl_sd

          poll_exceed_pct <- sum(test_ap$ap_ratio < threshold, na.rm = TRUE) / nrow(test_ap) * 100
          poll_signal <- if (poll_exceed_pct > 10) "EXCEEDANCE" else "NO_SIGNAL"
          poll_indicator <- "AP_decrease"
        }
      }

      # Check for crop taxa
      if ("variablename" %in% names(poll_samp)) {
        crops <- c("Cerealia", "Secale", "Triticum", "Hordeum", "Zea",
                    "Plantago lanceolata", "Plantago", "Cannabis", "Oryza")
        found_crops <- unique(poll_samp$variablename[poll_samp$variablename %in% crops])
        if (length(found_crops) > 0) {
          poll_indicator <- paste(poll_indicator, "+", paste(found_crops, collapse="/"))
        }
      }
    }

    exc_results[[length(exc_results) + 1]] <- data.frame(
      siteid = sid, sitename = site$sitename,
      lat = site$lat, lon = site$lon, region = site$region,
      char_exceed_pct = round(char_exceed_pct, 1),
      char_signal = char_signal,
      poll_exceed_pct = round(poll_exceed_pct, 1),
      poll_signal = poll_signal,
      poll_indicator = poll_indicator,
      stringsAsFactors = FALSE
    )

    Sys.sleep(0.3)
  }

  if (length(exc_results) > 0) {
    exc_df <- bind_rows(exc_results)
    saveRDS(exc_df, "/home/ayu/archeco/shared/cache/same_site_exceedance.rds")
    write.csv(exc_df, "/home/ayu/archeco/shared/same_site_exceedance.csv", row.names = FALSE)

    log("")
    log("### Exceedance Results")
    log("")
    log("| Site | Region | Char Exceed% | Char Signal | Pollen Exceed% | Pollen Signal | Indicator |")
    log("|------|--------|-------------|-------------|----------------|---------------|-----------|")
    for (i in seq_len(nrow(exc_df))) {
      r <- exc_df[i, ]
      log(sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                  substr(r$sitename, 1, 22), r$region,
                  if(!is.na(r$char_exceed_pct)) paste0(r$char_exceed_pct, "%") else "N/A",
                  r$char_signal,
                  if(!is.na(r$poll_exceed_pct)) paste0(r$poll_exceed_pct, "%") else "N/A",
                  r$poll_signal, r$poll_indicator))
    }

    # 2x2 contingency table
    log("")
    log("### 2x2 Contingency Table")
    log("")

    valid <- exc_df %>% filter(char_signal != "INSUFFICIENT" & poll_signal != "INSUFFICIENT")

    if (nrow(valid) >= 4) {
      tab <- table(
        Charcoal = valid$char_signal,
        Pollen = valid$poll_signal
      )
      log("```")
      log(paste(capture.output(print(tab)), collapse = "\n"))
      log("```")
      log(paste0("Total sites in 2x2: ", nrow(valid)))

      if (all(dim(tab) >= 2)) {
        ft <- fisher.test(tab)
        log(paste0("Fisher's exact test p = ", format(ft$p.value, digits = 4)))
        log(paste0("Odds ratio = ", format(ft$estimate, digits = 3)))
      }
    } else {
      log(paste0("Insufficient for 2x2 (need >=4 with both signals): ", nrow(valid), " available"))
    }
  }
}

# =============================================================================
# Save report
# =============================================================================
log("")
log("---")
log(paste0("Generated: ", Sys.time()))

writeLines(output_lines, "/home/ayu/archeco/shared/same_core_pollen_charcoal.md")
cat("\n=== Report saved to shared/same_core_pollen_charcoal.md ===\n")
cat("=== CSV saved to shared/same_site_pollen_charcoal.csv ===\n")
cat("=== Done ===\n")
