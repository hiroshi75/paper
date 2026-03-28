#!/usr/bin/env Rscript
# ============================================================================
# Eastern North America Charcoal Analysis
# Independent proxy validation for Paper 7 (reviewer response)
# ============================================================================
#
# Purpose: Query Neotoma for charcoal datasets in ENA, compare with pollen
# sites, and test whether charcoal shows fire management signals during
# the Pre-Columbian agricultural period.
#
# Data source: Neotoma Paleoecology Database (charcoal datasets)
# Region: Eastern North America (lat 25-50, lon -95 to -65)
# ============================================================================

library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)

# --- 1. Query Neotoma API for charcoal datasets ---
cat("Querying Neotoma API for charcoal datasets in ENA...\n")
url <- "https://api.neotomadb.org/v2.0/data/datasets"
resp <- GET(url, query = list(
  datasettype = "charcoal",
  loc = "POLYGON((-95 25,-65 25,-65 50,-95 50,-95 25))",
  limit = 500
))
data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
sites_df <- data$data$site

# Parse coordinates from GeoJSON geography field
coords <- do.call(rbind, lapply(1:nrow(sites_df), function(i) {
  geo <- tryCatch(fromJSON(sites_df$geography[i]), error = function(e) NULL)
  if (!is.null(geo) && !is.null(geo$coordinates)) {
    data.frame(lon = geo$coordinates[1], lat = geo$coordinates[2])
  } else {
    data.frame(lon = NA, lat = NA)
  }
}))
sites_df$lat <- coords$lat
sites_df$lon <- coords$lon
sites_df$datasetid <- sapply(sites_df$datasets, function(x) {
  if (is.data.frame(x) && "datasetid" %in% names(x)) x$datasetid[1] else NA
})

# Filter to ENA bounding box
ena <- sites_df[!is.na(sites_df$lat) & sites_df$lat >= 25 & sites_df$lat <= 50 &
                  sites_df$lon >= -95 & sites_df$lon <= -65, ]
cat("ENA charcoal sites:", nrow(ena), "\n")

# --- 2. Spatial overlap with pollen sites ---
pollen_sites <- readRDS("shared/cache/pollen_sites_locs.rds")

haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat / 2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
  2 * R * asin(sqrt(a))
}

overlap_results <- do.call(rbind, lapply(1:nrow(ena), function(i) {
  dists <- haversine(ena$lat[i], ena$lon[i], pollen_sites$lat, pollen_sites$long)
  nearest <- which.min(dists)
  data.frame(
    char_siteid = ena$siteid[i], char_sitename = ena$sitename[i],
    char_datasetid = ena$datasetid[i],
    char_lat = ena$lat[i], char_lon = ena$lon[i],
    pollen_siteid = pollen_sites$siteid[nearest],
    pollen_sitename = pollen_sites$sitename[nearest],
    dist_km = min(dists)
  )
}))

cat("Within 10km:", sum(overlap_results$dist_km <= 10), "\n")
cat("Within 50km:", sum(overlap_results$dist_km <= 50), "\n")

# --- 3. Download charcoal data ---
target_datasets <- unique(ena$datasetid)

download_charcoal <- function(dsid) {
  url <- paste0("https://api.neotomadb.org/v2.0/data/downloads/", dsid)
  resp <- tryCatch(GET(url, timeout(30)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NULL)

  d <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = FALSE)
  if (d$status != "success") return(NULL)

  site <- d$data$site
  samples <- tryCatch(site$collectionunit$dataset$samples[[1]], error = function(e) NULL)
  if (is.null(samples) || nrow(samples) == 0) return(NULL)

  rows <- list()
  for (i in 1:nrow(samples)) {
    ages_df <- samples$ages[[i]]
    datum_df <- samples$datum[[i]]

    age <- NA; agetype <- NA
    if (!is.null(ages_df) && is.data.frame(ages_df) && nrow(ages_df) > 0 && !all(is.na(ages_df$age))) {
      age <- ages_df$age[1]
      agetype <- ages_df$agetype[1]
    }

    if (!is.null(datum_df) && is.data.frame(datum_df)) {
      char_rows <- datum_df[datum_df$ecologicalgroup == "CHAR", ]
      if (nrow(char_rows) > 0) {
        for (j in 1:nrow(char_rows)) {
          rows[[length(rows) + 1]] <- data.frame(
            datasetid = dsid, siteid = site$siteid, sitename = site$sitename,
            sampleid = samples$sampleid[i], depth = samples$depth[i],
            age = age, agetype = agetype,
            value = char_rows$value[j], units = char_rows$units[j],
            element = char_rows$element[j], variablename = char_rows$variablename[j],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (length(rows) > 0) do.call(rbind, rows) else NULL
}

all_charcoal <- do.call(rbind, lapply(target_datasets, function(ds) {
  cat("  Downloading dataset", ds, "...\n")
  download_charcoal(ds)
}))

# --- 4. Period analysis ---
aged <- all_charcoal[!is.na(all_charcoal$age), ]
char_by_sample <- aged %>%
  group_by(datasetid, siteid, sitename, sampleid, depth, age, agetype) %>%
  summarize(total_charcoal = sum(value, na.rm = TRUE), .groups = "drop")

char_by_sample$period <- cut(char_by_sample$age,
  breaks = c(-Inf, 0, 200, 458, 2000, 5000, Inf),
  labels = c("Modern (<0 BP)", "European clearing (0-200 BP)",
             "Post-1492 (200-458 BP)", "Pre-Columbian agric (458-2000 BP)",
             "Pre-agriculture (2000-5000 BP)", "Deep past (>5000 BP)"),
  right = TRUE
)

# Summary statistics by period
period_stats <- char_by_sample %>%
  group_by(period) %>%
  summarize(
    n_samples = n(), n_sites = n_distinct(siteid),
    mean_charcoal = round(mean(total_charcoal), 1),
    median_charcoal = round(median(total_charcoal), 1),
    .groups = "drop"
  )
print(period_stats)

# Per-site agricultural vs pre-agricultural comparison
for (sid in unique(char_by_sample$siteid)) {
  site_data <- char_by_sample[char_by_sample$siteid == sid, ]
  sname <- unique(site_data$sitename)
  pre_ag <- site_data[site_data$age > 2000 & site_data$age <= 5000, ]
  agric <- site_data[site_data$age > 458 & site_data$age <= 2000, ]

  if (nrow(pre_ag) >= 3 && nrow(agric) >= 3) {
    ratio <- mean(agric$total_charcoal) / mean(pre_ag$total_charcoal)
    tt <- t.test(agric$total_charcoal, pre_ag$total_charcoal)
    cat(sprintf("%s: ratio=%.2f, p=%.4f\n", sname, ratio, tt$p.value))
  }
}
