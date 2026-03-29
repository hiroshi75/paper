#!/usr/bin/env Rscript
# =============================================================================
# Pollen-Charcoal Cross-Proxy Analysis
# Question: At sites where pollen shows NO agricultural exceedance (H3),
# does charcoal show exceedance? → fire management detectable by charcoal alone?
# =============================================================================

library(dplyr)
library(tidyr)

cat("============================================================\n")
cat("Pollen-Charcoal Cross-Proxy Analysis\n")
cat("============================================================\n\n")

# =============================================================================
# STEP 1: Load all datasets
# =============================================================================
cat("--- Step 1: Loading datasets ---\n\n")

# --- Charcoal data ---
# Pilot global charcoal (raw samples - need to derive exceedance + get coords)
pilot_char_raw <- readRDS("shared/cache/pilot_charcoal_global.rds")
pilot_char_exc <- readRDS("shared/cache/pilot_charcoal_exceedance_results.rds")

# ENA charcoal sites with coordinates
ena_char_sites <- readRDS("shared/cache/ena_charcoal_sites_df.rds")
ena_char_filtered <- readRDS("shared/cache/ena_charcoal_filtered.rds")
ena_char_all <- readRDS("shared/cache/ena_charcoal_all_data.rds")

cat("Pilot charcoal exceedance results:", nrow(pilot_char_exc), "sites\n")
cat("ENA charcoal sites (df):", nrow(ena_char_sites), "sites\n")
cat("ENA charcoal filtered:", nrow(ena_char_filtered), "sites\n")

# --- European pollen data ---
pollen_uk <- read.csv("shared/pilot_signal_phase_uk.csv", stringsAsFactors = FALSE)
pollen_scand <- read.csv("shared/signal_phase_scandinavia.csv", stringsAsFactors = FALSE)
pollen_ce <- read.csv("shared/signal_phase_central_europe.csv", stringsAsFactors = FALSE)

# --- ENA pollen data ---
pollen_ena <- read.csv("shared/full_ena_exceedance_data.csv", stringsAsFactors = FALSE)

cat("Pollen UK:", nrow(pollen_uk), "datasets\n")
cat("Pollen Scandinavia:", nrow(pollen_scand), "datasets\n")
cat("Pollen Central Europe:", nrow(pollen_ce), "datasets\n")
cat("Pollen ENA:", nrow(pollen_ena), "sites\n")

# =============================================================================
# Standardize pollen datasets with H1/H3 classification
# =============================================================================
cat("\n--- Standardizing pollen classifications ---\n")

# For European data: classify based on available columns
# "persistent" = TRUE and signal detected → H1 (agricultural/persistent change)
# signal detected but not persistent → intermediate
# no signal → H3 (natural)

# UK: no 'persistent' column, use pct_exceed_after
# High exceedance (>50%) with onset in Neolithic range → H1-like
standardize_uk <- function(df) {
  df %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    mutate(
      has_signal = !is.na(signal_onset_age),
      # Classify: if signal onset exists and persistent exceedance >50%, call it H1
      classification = case_when(
        !has_signal ~ "H3_natural",
        pct_exceed_after >= 50 ~ "H1_agricultural",
        TRUE ~ "H3_natural"
      ),
      region = "UK_Ireland"
    ) %>%
    select(datasetid, sitename, lat, lon, classification, signal_onset_age, region)
}

standardize_scand <- function(df) {
  df %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    mutate(
      has_signal = !is.na(signal_onset_age),
      classification = case_when(
        !has_signal ~ "H3_natural",
        persistent == TRUE ~ "H1_agricultural",
        TRUE ~ "H3_natural"
      ),
      region = "Scandinavia"
    ) %>%
    select(datasetid, sitename, lat, lon, classification, signal_onset_age, region)
}

standardize_ce <- function(df) {
  df %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    mutate(
      has_signal = !is.na(signal_onset_age),
      classification = case_when(
        !has_signal ~ "H3_natural",
        persistent == TRUE ~ "H1_agricultural",
        TRUE ~ "H3_natural"
      ),
      region = "Central_Europe"
    ) %>%
    select(datasetid, sitename, lat, lon, classification, signal_onset_age, region)
}

standardize_ena <- function(df) {
  df %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    mutate(
      region = "ENA",
      datasetid = siteid  # ENA uses siteid as identifier
    ) %>%
    select(datasetid, sitename, lat, lon, classification,
           signal_onset_age = onset_age_bp, region)
}

pollen_std_uk <- standardize_uk(pollen_uk)
pollen_std_scand <- standardize_scand(pollen_scand)
pollen_std_ce <- standardize_ce(pollen_ce)
pollen_std_ena <- standardize_ena(pollen_ena)

# Combine all pollen
pollen_all <- bind_rows(pollen_std_uk, pollen_std_scand, pollen_std_ce, pollen_std_ena)

cat("Total pollen sites with coordinates:\n")
print(table(pollen_all$region, pollen_all$classification))
cat("Total:", nrow(pollen_all), "\n")

# =============================================================================
# Build unified charcoal site table with coordinates
# =============================================================================
cat("\n--- Building charcoal site table ---\n")

# Source 1: ENA charcoal sites (have lat/lon from sites_df)
ena_char_coords <- ena_char_sites %>%
  select(siteid, sitename, lat, lon = long) %>%
  filter(!is.na(lat), !is.na(lon))

# Source 2: ENA charcoal filtered (more sites with coords)
ena_char_filt_coords <- ena_char_filtered %>%
  select(siteid, sitename, lat, lon) %>%
  filter(!is.na(lat), !is.na(lon))

# Source 3: Pilot charcoal exceedance (has region but no coords)
# These are from Neotoma - need to find coords from siteid
# The pilot_char_exc has siteids that may match Neotoma sites

# Combine ENA charcoal coords
char_coords <- bind_rows(ena_char_coords, ena_char_filt_coords) %>%
  distinct(siteid, .keep_all = TRUE)

cat("ENA charcoal sites with coords:", nrow(char_coords), "\n")

# For pilot charcoal sites, try to get coords from Neotoma site names
# that match Scandinavia pollen sites (same Neotoma DB)
pilot_sites <- pilot_char_exc %>%
  select(siteid, sitename, region)

# Check if any pilot sites match pollen sites by name
cat("Pilot charcoal sites:\n")
print(pilot_sites[, c("siteid", "sitename", "region")])

# Try to match pilot charcoal sites to Scandinavian pollen sites by name
# (many pilot European sites are in Scandinavia)
scand_coords <- pollen_scand %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  select(sitename, lat, lon) %>%
  distinct(sitename, .keep_all = TRUE)

# Check for name matches
pilot_in_scand <- pilot_sites %>%
  inner_join(scand_coords, by = "sitename")

cat("\nPilot charcoal sites found in Scandinavian pollen data:", nrow(pilot_in_scand), "\n")
if (nrow(pilot_in_scand) > 0) print(pilot_in_scand)

# Also try matching against all pollen data
all_pollen_coords <- pollen_all %>%
  select(sitename, lat, lon) %>%
  distinct(sitename, .keep_all = TRUE)

pilot_in_pollen <- pilot_sites %>%
  inner_join(all_pollen_coords, by = "sitename")

cat("Pilot charcoal sites found in ANY pollen data:", nrow(pilot_in_pollen), "\n")
if (nrow(pilot_in_pollen) > 0) print(pilot_in_pollen)

# For pilot sites not in pollen data, try to get coords via Neotoma API
# But first, let's check: the pilot exceedance already tells us
# which sites have exceedance. We need coords to do spatial matching.
# Let's try to query Neotoma for the remaining sites.

# Actually, many pilot charcoal sites ARE pollen sites (same Neotoma site).
# The names like "Barheivatn", "Brurskardtjørni", "Dalmutladdo" appear in
# Scandinavian pollen data. These are CO-LOCATED by definition (same core).

# Let's check this systematically
cat("\n--- Checking same-core overlap (charcoal site = pollen site) ---\n")

# Pilot charcoal sites that appear in pollen datasets
same_core_scand <- pilot_char_exc %>%
  filter(region == "Europe") %>%
  inner_join(
    pollen_scand %>% select(sitename, datasetid, lat, lon, signal_onset_age, persistent),
    by = "sitename"
  )

cat("Same-core matches (pilot charcoal in Scandinavia pollen):", nrow(same_core_scand), "\n")
if (nrow(same_core_scand) > 0) {
  cat("Matched sites:\n")
  print(same_core_scand %>% select(siteid, sitename, lat, lon, region))
}

# For North America pilot sites
same_core_ena <- pilot_char_exc %>%
  filter(region == "North America") %>%
  inner_join(
    pollen_ena %>% select(sitename, siteid, lat, lon, classification, onset_age_bp),
    by = "sitename"
  )

cat("\nSame-core matches (pilot charcoal in ENA pollen):", nrow(same_core_ena), "\n")

# =============================================================================
# STEP 2: Find co-located sites (within 10km and 50km)
# =============================================================================
cat("\n\n--- Step 2: Finding co-located sites ---\n")

# Haversine distance function
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371  # Earth radius in km
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1 * pi/180) * cos(lat2 * pi/180) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R * c
}

# Strategy: Use ALL available charcoal site coordinates
# 1. Same-core matches from pilot (charcoal + pollen in same Neotoma site)
# 2. ENA charcoal sites matched to ENA pollen sites by proximity
# 3. Existing overlap data from cache

# Build comprehensive charcoal exceedance status
# For pilot charcoal: we have exceedance results directly
pilot_exc_status <- pilot_char_exc %>%
  mutate(
    char_has_exceedance = pct_exceed > 0,
    char_pct_exceed = pct_exceed,
    char_source = "pilot_global"
  )

# For ENA charcoal: need to compute exceedance from raw data
# Use the ena_charcoal_all_data with a simple baseline exceedance approach
cat("\n--- Computing ENA charcoal exceedance ---\n")

ena_char_exc_list <- list()
for (sid in unique(ena_char_all$siteid)) {
  site_data <- ena_char_all %>%
    filter(siteid == sid, !is.na(age), !is.na(value)) %>%
    arrange(desc(age))

  if (nrow(site_data) < 10) next

  # Use pre-5000 BP as baseline, post-5000 BP as agricultural period
  # (rough Neolithic boundary for ENA)
  baseline <- site_data %>% filter(age > 5000)
  agric <- site_data %>% filter(age <= 5000)

  if (nrow(baseline) < 3 || nrow(agric) < 3) next

  bl_mean <- mean(baseline$value, na.rm = TRUE)
  bl_sd <- sd(baseline$value, na.rm = TRUE)
  threshold <- bl_mean + 2 * bl_sd

  n_exc <- sum(agric$value > threshold, na.rm = TRUE)
  pct_exc <- n_exc / nrow(agric) * 100

  ena_char_exc_list[[as.character(sid)]] <- data.frame(
    siteid = sid,
    sitename = site_data$sitename[1],
    n_baseline = nrow(baseline),
    n_agric = nrow(agric),
    baseline_mean = bl_mean,
    baseline_sd = bl_sd,
    threshold = threshold,
    n_exceed = n_exc,
    pct_exceed = pct_exc,
    char_has_exceedance = pct_exc > 0,
    char_source = "ena_charcoal"
  )
}

ena_char_exc <- bind_rows(ena_char_exc_list)
cat("ENA charcoal exceedance computed for", nrow(ena_char_exc), "sites\n")
if (nrow(ena_char_exc) > 0) print(ena_char_exc)

# Add coordinates to ENA charcoal exceedance
ena_char_exc_geo <- ena_char_exc %>%
  left_join(char_coords %>% select(siteid, lat, lon), by = "siteid") %>%
  filter(!is.na(lat), !is.na(lon))

cat("ENA charcoal sites with coords and exceedance:", nrow(ena_char_exc_geo), "\n")

# =============================================================================
# Now find nearest pollen site for each charcoal site
# =============================================================================
cat("\n--- Finding nearest pollen sites ---\n")

# Combine all charcoal sites with coordinates and exceedance status
# Source A: Same-core matches (pilot charcoal in Scandinavian pollen)
if (nrow(same_core_scand) > 0) {
  same_core_pairs <- same_core_scand %>%
    mutate(
      char_siteid = siteid,
      char_sitename = sitename,
      char_lat = lat,
      char_lon = lon,
      char_has_exceedance = pct_exceed > 0,
      char_pct_exceed = pct_exceed,
      pollen_sitename = sitename,
      pollen_classification = case_when(
        is.na(signal_onset_age) ~ "H3_natural",
        persistent == TRUE ~ "H1_agricultural",
        TRUE ~ "H3_natural"
      ),
      pollen_signal_onset = signal_onset_age,
      dist_km = 0,
      match_type = "same_core"
    ) %>%
    select(char_siteid, char_sitename, char_lat, char_lon,
           char_has_exceedance, char_pct_exceed,
           pollen_sitename, pollen_classification, pollen_signal_onset,
           dist_km, match_type)

  cat("Same-core pairs:", nrow(same_core_pairs), "\n")
} else {
  same_core_pairs <- data.frame()
}

# Source B: Proximity-based matching for ENA charcoal → ENA pollen
proximity_pairs <- data.frame()

if (nrow(ena_char_exc_geo) > 0) {
  for (i in 1:nrow(ena_char_exc_geo)) {
    char_site <- ena_char_exc_geo[i, ]

    # Calculate distance to all ENA pollen sites
    dists <- haversine_km(
      char_site$lat, char_site$lon,
      pollen_std_ena$lat, pollen_std_ena$lon
    )

    nearest_idx <- which.min(dists)
    nearest_dist <- dists[nearest_idx]
    nearest_pollen <- pollen_std_ena[nearest_idx, ]

    proximity_pairs <- bind_rows(proximity_pairs, data.frame(
      char_siteid = char_site$siteid,
      char_sitename = char_site$sitename,
      char_lat = char_site$lat,
      char_lon = char_site$lon,
      char_has_exceedance = char_site$char_has_exceedance,
      char_pct_exceed = char_site$pct_exceed,
      pollen_sitename = nearest_pollen$sitename,
      pollen_classification = nearest_pollen$classification,
      pollen_signal_onset = nearest_pollen$signal_onset_age,
      dist_km = nearest_dist,
      match_type = "proximity_ena"
    ))
  }
}

# Source C: Proximity matching for pilot charcoal sites that have coords
# but aren't same-core matches
# For European pilot sites NOT in same-core, match to nearest European pollen
pilot_europe_unmatched <- pilot_char_exc %>%
  filter(region == "Europe") %>%
  filter(!(sitename %in% same_core_scand$sitename))

if (nrow(pilot_europe_unmatched) > 0) {
  # Try to get coords from pollen data
  for (i in 1:nrow(pilot_europe_unmatched)) {
    sname <- pilot_europe_unmatched$sitename[i]
    # Check if this site name appears anywhere in pollen data
    match_pollen <- pollen_all %>% filter(sitename == sname)
    if (nrow(match_pollen) > 0) {
      proximity_pairs <- bind_rows(proximity_pairs, data.frame(
        char_siteid = pilot_europe_unmatched$siteid[i],
        char_sitename = sname,
        char_lat = match_pollen$lat[1],
        char_lon = match_pollen$lon[1],
        char_has_exceedance = pilot_europe_unmatched$pct_exceed[i] > 0,
        char_pct_exceed = pilot_europe_unmatched$pct_exceed[i],
        pollen_sitename = match_pollen$sitename[1],
        pollen_classification = match_pollen$classification[1],
        pollen_signal_onset = match_pollen$signal_onset_age[1],
        dist_km = 0,
        match_type = "same_core_name_match"
      ))
    }
  }
}

# Also use the pre-computed overlap from cache
cached_overlap <- readRDS("shared/cache/charcoal_pollen_overlap.rds")
cached_overlap_full <- readRDS("shared/cache/charcoal_pollen_overlap_full.rds")

cat("Cached overlap pairs:", nrow(cached_overlap), "\n")
cat("Cached overlap full pairs:", nrow(cached_overlap_full), "\n")

# For cached overlap, add exceedance status
# These are ENA sites - check against our ENA charcoal exceedance
cached_with_exc <- cached_overlap %>%
  left_join(ena_char_exc %>% select(siteid, pct_exceed, char_has_exceedance),
            by = c("char_siteid" = "siteid"))

# Match pollen classification
cached_with_class <- cached_with_exc %>%
  left_join(pollen_std_ena %>% select(sitename, classification, signal_onset_age),
            by = c("pollen_sitename" = "sitename"))

# Combine all pairs
all_pairs <- bind_rows(
  same_core_pairs,
  proximity_pairs
)

# Deduplicate by charcoal site
all_pairs <- all_pairs %>%
  arrange(dist_km) %>%
  distinct(char_siteid, .keep_all = TRUE)

cat("\n=== Co-location Results ===\n")
cat("Total unique charcoal-pollen pairs:", nrow(all_pairs), "\n")
cat("  Within 0 km (same core):", sum(all_pairs$dist_km == 0), "\n")
cat("  Within 10 km:", sum(all_pairs$dist_km <= 10), "\n")
cat("  Within 50 km:", sum(all_pairs$dist_km <= 50), "\n")
cat("  Within 100 km:", sum(all_pairs$dist_km <= 100), "\n")

cat("\nAll pairs:\n")
print(all_pairs %>% select(char_sitename, pollen_sitename, dist_km,
                            char_has_exceedance, pollen_classification, match_type))

# =============================================================================
# STEP 3: Cross-exceedance analysis (2x2 table)
# =============================================================================
cat("\n\n--- Step 3: Cross-exceedance 2x2 table ---\n")

# Use pairs within 50km for the primary analysis
pairs_50km <- all_pairs %>% filter(dist_km <= 50)

cat("Pairs within 50km:", nrow(pairs_50km), "\n\n")

if (nrow(pairs_50km) > 0) {
  # Build 2x2 table
  cross_tab <- table(
    Pollen = ifelse(pairs_50km$pollen_classification == "H1_agricultural",
                    "H1_agricultural", "H3_natural_or_none"),
    Charcoal = ifelse(pairs_50km$char_has_exceedance,
                      "Exceedance", "No_exceedance")
  )

  cat("2x2 Cross-exceedance table (pairs within 50km):\n")
  print(cross_tab)

  cat("\n")

  # Cell labels
  if ("H1_agricultural" %in% rownames(cross_tab) && "Exceedance" %in% colnames(cross_tab)) {
    A <- cross_tab["H1_agricultural", "Exceedance"]
  } else A <- 0
  if ("H1_agricultural" %in% rownames(cross_tab) && "No_exceedance" %in% colnames(cross_tab)) {
    B <- cross_tab["H1_agricultural", "No_exceedance"]
  } else B <- 0
  if ("H3_natural_or_none" %in% rownames(cross_tab) && "Exceedance" %in% colnames(cross_tab)) {
    C <- cross_tab["H3_natural_or_none", "Exceedance"]
  } else C <- 0
  if ("H3_natural_or_none" %in% rownames(cross_tab) && "No_exceedance" %in% colnames(cross_tab)) {
    D <- cross_tab["H3_natural_or_none", "No_exceedance"]
  } else D <- 0

  cat("Cell A (H1 pollen + charcoal exceedance):", A, "\n")
  cat("Cell B (H1 pollen + no charcoal exceedance):", B, "\n")
  cat("Cell C (H3 pollen + charcoal exceedance) *** KEY:", C, "\n")
  cat("Cell D (H3 pollen + no charcoal exceedance):", D, "\n")

  # Also do the analysis for ALL pairs regardless of distance
  cat("\n\n--- 2x2 table for ALL pairs (any distance) ---\n")
  cross_tab_all <- table(
    Pollen = ifelse(all_pairs$pollen_classification == "H1_agricultural",
                    "H1_agricultural", "H3_natural_or_none"),
    Charcoal = ifelse(all_pairs$char_has_exceedance,
                      "Exceedance", "No_exceedance")
  )
  print(cross_tab_all)

  # Detail the Cell C sites
  cell_c_sites <- all_pairs %>%
    filter(pollen_classification != "H1_agricultural",
           char_has_exceedance == TRUE)

  cat("\n\n--- Cell C detail: Pollen H3/none BUT charcoal exceedance ---\n")
  cat("Number of Cell C sites:", nrow(cell_c_sites), "\n")
  if (nrow(cell_c_sites) > 0) {
    print(cell_c_sites %>% select(char_sitename, pollen_sitename, dist_km,
                                   char_pct_exceed, pollen_classification))
  }
}

# =============================================================================
# STEP 4: Temporal comparison at co-located sites
# =============================================================================
cat("\n\n--- Step 4: Temporal comparison at co-located sites ---\n")

# For Cell C sites: when does charcoal exceedance begin?
# We need the raw charcoal time series
if (nrow(cell_c_sites) > 0) {
  cat("Analyzing temporal patterns for Cell C sites:\n\n")

  for (i in 1:nrow(cell_c_sites)) {
    sid <- cell_c_sites$char_siteid[i]
    sname <- cell_c_sites$char_sitename[i]

    cat("Site:", sname, "(siteid:", sid, ")\n")

    # Check if in pilot charcoal raw data
    site_raw <- pilot_char_raw %>% filter(siteid == sid) %>% arrange(desc(age))

    if (nrow(site_raw) == 0) {
      # Check ENA charcoal
      site_raw <- ena_char_all %>% filter(siteid == sid) %>% arrange(desc(age))
    }

    if (nrow(site_raw) > 0) {
      # Find when values first exceed baseline
      bl <- site_raw %>% filter(age > 5000)
      if (nrow(bl) >= 3) {
        bl_mean <- mean(bl$value, na.rm = TRUE)
        bl_sd <- sd(bl$value, na.rm = TRUE)
        threshold <- bl_mean + 2 * bl_sd

        exceeding <- site_raw %>% filter(value > threshold, age <= 7000) %>% arrange(desc(age))

        if (nrow(exceeding) > 0) {
          first_exc_age <- max(exceeding$age, na.rm = TRUE)
          cat("  First charcoal exceedance:", round(first_exc_age), "cal BP\n")

          # Neolithic in Europe ~6000-4000 BP, in ENA ~3000-1000 BP
          if (first_exc_age > 6000) {
            cat("  → Pre-Neolithic: likely natural fire regime\n")
          } else if (first_exc_age > 3000) {
            cat("  → Neolithic period: possible fire management\n")
          } else {
            cat("  → Post-Neolithic: late fire management or land clearing\n")
          }
        }
      }
    } else {
      cat("  No raw time series available\n")
    }

    # Exceedance from pilot results
    pilot_match <- pilot_char_exc %>% filter(siteid == sid)
    if (nrow(pilot_match) > 0) {
      cat("  Charcoal exceedance:", round(pilot_match$pct_exceed[1], 1), "%\n")
      cat("  Ratio (agric/baseline):", round(pilot_match$ratio[1], 1), "x\n")
    }

    cat("  Pollen classification:", cell_c_sites$pollen_classification[i], "\n\n")
  }
}

# =============================================================================
# STEP 5: Complementary detection test
# =============================================================================
cat("\n--- Step 5: Complementary detection test ---\n\n")

# Calculate: what fraction of total detected human impact would be
# MISSED by pollen alone but RECOVERED by adding charcoal?

total_pairs <- nrow(all_pairs)
pollen_detected <- sum(all_pairs$pollen_classification == "H1_agricultural")
charcoal_detected <- sum(all_pairs$char_has_exceedance)
both_detected <- sum(all_pairs$pollen_classification == "H1_agricultural" &
                       all_pairs$char_has_exceedance)
pollen_only <- sum(all_pairs$pollen_classification == "H1_agricultural" &
                     !all_pairs$char_has_exceedance)
charcoal_only <- sum(all_pairs$pollen_classification != "H1_agricultural" &
                       all_pairs$char_has_exceedance)
neither <- sum(all_pairs$pollen_classification != "H1_agricultural" &
                 !all_pairs$char_has_exceedance)

cat("Total co-located pairs:", total_pairs, "\n")
cat("Pollen detects human impact (H1):", pollen_detected, "\n")
cat("Charcoal detects exceedance:", charcoal_detected, "\n")
cat("Both detect:", both_detected, "\n")
cat("Pollen only:", pollen_only, "\n")
cat("Charcoal only (Cell C):", charcoal_only, "\n")
cat("Neither:", neither, "\n\n")

# Union of detected sites
union_detected <- both_detected + pollen_only + charcoal_only

if (union_detected > 0) {
  pollen_alone_coverage <- (both_detected + pollen_only) / union_detected * 100
  added_by_charcoal <- charcoal_only / union_detected * 100

  cat("=== Complementary Detection Summary ===\n")
  cat("Sites with ANY detected impact:", union_detected, "\n")
  cat("Coverage by pollen alone:", round(pollen_alone_coverage, 1), "%\n")
  cat("Additional coverage from charcoal:", round(added_by_charcoal, 1), "%\n")
  cat("Fraction MISSED by pollen, RECOVERED by charcoal:",
      round(added_by_charcoal, 1), "%\n\n")

  if (charcoal_only > 0) {
    cat("*** FINDING: Adding charcoal recovers", charcoal_only,
        "site(s) (", round(added_by_charcoal, 1),
        "%) of total detected human impact that pollen alone would miss.\n")
  } else {
    cat("No sites where charcoal detects impact that pollen misses.\n")
  }
} else {
  cat("No impact detected at any co-located site by either proxy.\n")
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================
cat("\n\n============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n\n")

cat("Data inventory:\n")
cat("  European pollen sites:", nrow(pollen_std_uk) + nrow(pollen_std_scand) + nrow(pollen_std_ce), "\n")
cat("  ENA pollen sites:", nrow(pollen_std_ena), "\n")
cat("  Charcoal sites (pilot global):", nrow(pilot_char_exc), "\n")
cat("  Charcoal sites (ENA):", nrow(ena_char_exc), "\n")

cat("\nCo-location:\n")
cat("  Same-core matches:", sum(all_pairs$dist_km == 0), "\n")
cat("  Within 10km:", sum(all_pairs$dist_km <= 10), "\n")
cat("  Within 50km:", sum(all_pairs$dist_km <= 50), "\n")

cat("\n2x2 Table (all co-located pairs):\n")
cat("                          Charcoal+   Charcoal-\n")
cat(sprintf("  Pollen H1 (agricultural)    %d           %d\n",
            sum(all_pairs$pollen_classification == "H1_agricultural" & all_pairs$char_has_exceedance),
            sum(all_pairs$pollen_classification == "H1_agricultural" & !all_pairs$char_has_exceedance)))
cat(sprintf("  Pollen H3 (natural/none)    %d           %d\n",
            sum(all_pairs$pollen_classification != "H1_agricultural" & all_pairs$char_has_exceedance),
            sum(all_pairs$pollen_classification != "H1_agricultural" & !all_pairs$char_has_exceedance)))

cat("\n============================================================\n")
