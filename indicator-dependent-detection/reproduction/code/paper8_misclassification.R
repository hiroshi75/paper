#!/usr/bin/env Rscript
# =============================================================================
# Paper 8: Misclassification Analysis (Analysis 9.4)
# =============================================================================
# Quantifies the "epistemic error rate" — what fraction of ENA sites are
# misclassified as "no agricultural impact" when European-only indicators
# are applied.
#
# Author: ArchEco research agent
# Date: 2026-03-28
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# --- Load indicator sensitivity results ---
results_file <- "/home/ayu/archeco/shared/cache/paper8_indicator_sensitivity_results.rds"
OUT_DIR <- "/home/ayu/archeco/shared"

if (!file.exists(results_file)) {
  stop("Run paper8_indicator_sensitivity.R first!")
}

all_results <- readRDS(results_file)
CONTACT_YEAR <- 458  # 1492 CE in cal BP

cat("========================================================\n")
cat("PAPER 8: MISCLASSIFICATION ANALYSIS (Analysis 9.4)\n")
cat("========================================================\n\n")

# --- Extract results by set ---
set_a <- all_results$A_European   # European indicators
set_b <- all_results$B_Indigenous # Indigenous (EAC)
set_c <- all_results$C_Disturbance # Disturbance
set_d <- all_results$D_Combined   # Combined

n_total <- nrow(set_a)
cat(sprintf("Total ENA sites analyzed: %d\n\n", n_total))

# ====================================================================
# STEP 1: European-only classification
# ====================================================================
cat("=== STEP 1: European-only Classification ===\n")

# Sites where European indicators say "no impact"
# Two definitions:
# (a) Strict: taxa not even present, or present but no exceedance
# (b) Lenient: testable but no exceedance
euro_no_impact_strict <- set_a %>% filter(!exceedance)
euro_no_impact_testable <- set_a %>% filter(testable & !exceedance)
euro_signal <- set_a %>% filter(exceedance)

n_euro_no_impact <- nrow(euro_no_impact_strict)
n_euro_signal <- nrow(euro_signal)

cat(sprintf("European indicators say 'signal detected': %d (%.1f%%)\n",
            n_euro_signal, n_euro_signal / n_total * 100))
cat(sprintf("European indicators say 'no impact': %d (%.1f%%)\n",
            n_euro_no_impact, n_euro_no_impact / n_total * 100))

# Note: Some "European" detections are from native Plantago species
# Not true agricultural signals — this inflates Set A
euro_precontact_only <- set_a %>% filter(exceedance & pre_contact & !post_contact)
cat(sprintf("  (of which pre-contact only — likely native Plantago: %d)\n",
            nrow(euro_precontact_only)))

# ====================================================================
# STEP 2: How many "no impact" sites actually show indigenous signal?
# ====================================================================
cat("\n=== STEP 2: Indigenous Signal Among 'No Impact' Sites ===\n")

# Join European and Indigenous results
combined <- set_a %>%
  select(siteid, sitename, lat, long,
         euro_exceedance = exceedance, euro_testable = testable,
         euro_onset = onset_age_bp) %>%
  left_join(
    set_b %>% select(siteid,
                     indig_exceedance = exceedance, indig_testable = testable,
                     indig_onset = onset_age_bp, indig_pre_contact = pre_contact),
    by = "siteid"
  ) %>%
  left_join(
    set_c %>% select(siteid,
                     disturb_exceedance = exceedance, disturb_testable = testable,
                     disturb_onset = onset_age_bp, disturb_pre_contact = pre_contact),
    by = "siteid"
  ) %>%
  left_join(
    set_d %>% select(siteid,
                     combined_exceedance = exceedance, combined_testable = testable,
                     combined_onset = onset_age_bp),
    by = "siteid"
  )

# Sites European says "no impact" but indigenous says "signal"
misclassified_indig <- combined %>%
  filter(!euro_exceedance & indig_exceedance)

n_misclass_indig <- nrow(misclassified_indig)
misclass_rate_indig <- n_misclass_indig / n_total * 100

cat(sprintf("Sites classified 'no impact' by European: %d\n", n_euro_no_impact))
cat(sprintf("Of those, showing Indigenous (EAC) exceedance: %d\n", n_misclass_indig))
cat(sprintf("Misclassification rate (as %% of all sites): %.1f%%\n", misclass_rate_indig))
cat(sprintf("Misclassification rate (as %% of 'no impact' sites): %.1f%%\n",
            n_misclass_indig / n_euro_no_impact * 100))

# Pre-contact indigenous signal (stronger evidence of actual agriculture)
misclass_precontact <- misclassified_indig %>%
  filter(indig_pre_contact == TRUE)
cat(sprintf("Of misclassified, with pre-contact indigenous signal: %d\n",
            nrow(misclass_precontact)))

# Sites European says "no impact" but disturbance says "signal"
misclassified_disturb <- combined %>%
  filter(!euro_exceedance & disturb_exceedance)
n_misclass_disturb <- nrow(misclassified_disturb)
cat(sprintf("\nSites 'no impact' by European but disturbance exceedance: %d (%.1f%%)\n",
            n_misclass_disturb, n_misclass_disturb / n_total * 100))

# Sites European misses but ANY other set detects
misclassified_any <- combined %>%
  filter(!euro_exceedance & (indig_exceedance | disturb_exceedance))
n_misclass_any <- nrow(misclassified_any)
cat(sprintf("Sites 'no impact' by European but detected by ANY other: %d (%.1f%%)\n",
            n_misclass_any, n_misclass_any / n_total * 100))

# ====================================================================
# STEP 3: Full Cross-Classification
# ====================================================================
cat("\n=== STEP 3: Cross-Classification Matrix ===\n\n")

# Create cross-classification
combined <- combined %>%
  mutate(
    euro_class = ifelse(euro_exceedance, "Detected", "No impact"),
    indig_class = ifelse(indig_exceedance, "Detected",
                         ifelse(is.na(indig_testable) | !indig_testable,
                                "Not testable", "No impact")),
    disturb_class = ifelse(disturb_exceedance, "Detected",
                           ifelse(is.na(disturb_testable) | !disturb_testable,
                                  "Not testable", "No impact")),
    combined_class = ifelse(combined_exceedance, "Detected",
                            ifelse(is.na(combined_testable) | !combined_testable,
                                   "Not testable", "No impact"))
  )

cat("European x Indigenous cross-tabulation:\n")
print(table(European = combined$euro_class, Indigenous = combined$indig_class))

cat("\nEuropean x Disturbance cross-tabulation:\n")
print(table(European = combined$euro_class, Disturbance = combined$disturb_class))

cat("\nEuropean x Combined cross-tabulation:\n")
print(table(European = combined$euro_class, Combined = combined$combined_class))

# ====================================================================
# STEP 4: Detection Matrix (Paper's key figure)
# ====================================================================
cat("\n=== STEP 4: Detection Matrix ===\n\n")

# Classify each site by what each framework detects
detection_categories <- combined %>%
  mutate(
    category = case_when(
      euro_exceedance & indig_exceedance ~ "Both detect",
      euro_exceedance & !indig_exceedance ~ "European only",
      !euro_exceedance & indig_exceedance ~ "Indigenous only (MISSED by European)",
      !euro_exceedance & !indig_exceedance ~ "Neither detects"
    )
  )

cat("Detection categories (European vs Indigenous):\n")
cat_table <- table(detection_categories$category)
for (nm in names(cat_table)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", nm, cat_table[nm],
              cat_table[nm] / n_total * 100))
}

# ====================================================================
# STEP 5: Extrapolation to Global Context
# ====================================================================
cat("\n=== STEP 5: Global Extrapolation ===\n\n")

# ENA empirical rates
cat("EMPIRICAL (ENA):\n")
cat(sprintf("  Sites with European-only framework: %d 'no impact'\n", n_euro_no_impact))
cat(sprintf("  Of those, with indigenous agricultural signal: %d\n", n_misclass_indig))
cat(sprintf("  Misclassification rate: %.1f%%\n\n", misclass_rate_indig))

cat("ANALOGOUS GLOBAL SYSTEMS (cautious extrapolation):\n\n")

global_systems <- data.frame(
  Region = c("Eastern North America",
             "Sub-Saharan Africa",
             "South/Southeast Asia",
             "Amazonia",
             "Oceania/Pacific"),
  Agricultural_System = c("EAC (Chenopodium, Amaranthus, sunflower, squash)",
                           "Yam, millet, sorghum, teff, enset cultivation",
                           "Rice in forest, taro, yam, breadfruit",
                           "Manioc, managed forest, terra preta",
                           "Taro, yam, breadfruit, coconut"),
  European_Indicator_Applicable = c("No", "No", "Partially (rice)",
                                     "No", "No"),
  Indigenous_Indicators_Developed = c("Yes (this study)",
                                       "Limited (pearl millet pollen)",
                                       "Partially (Oryza pollen)",
                                       "Very limited",
                                       "Limited"),
  Estimated_Blind_Spot = c(sprintf("%.1f%% (measured)", misclass_rate_indig),
                            "Unknown but plausibly >20%",
                            "Unknown, complex",
                            "Unknown but plausibly high",
                            "Unknown"),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(global_systems)) {
  gs <- global_systems[i, ]
  cat(sprintf("  %s:\n", gs$Region))
  cat(sprintf("    System: %s\n", gs$Agricultural_System))
  cat(sprintf("    European indicators applicable: %s\n", gs$European_Indicator_Applicable))
  cat(sprintf("    Indigenous indicators developed: %s\n", gs$Indigenous_Indicators_Developed))
  cat(sprintf("    Estimated blind spot: %s\n\n", gs$Estimated_Blind_Spot))
}

# ====================================================================
# STEP 6: Generate Markdown Report
# ====================================================================
cat("--- Generating report ---\n")

md <- c(
  "# Paper 8: Misclassification Analysis",
  "",
  sprintf("**Date**: %s", Sys.Date()),
  sprintf("**Total ENA sites**: %d", n_total),
  "",
  "## The Epistemic Error Rate",
  "",
  "If a researcher applies ONLY European agricultural indicators (Set A: *Plantago lanceolata*,",
  "*Rumex*, Cerealia-type, *Secale*) to Eastern North American pollen records, what fraction",
  "of agricultural signals are missed?",
  "",
  "## Step 1: European-Only Classification",
  "",
  sprintf("- European indicators detect exceedance at: **%d sites (%.1f%%)**", n_euro_signal, n_euro_signal/n_total*100),
  sprintf("- European indicators say 'no impact': **%d sites (%.1f%%)**", n_euro_no_impact, n_euro_no_impact/n_total*100),
  sprintf("- Note: %d European 'detections' are pre-contact only (likely native *Plantago*, not agriculture)",
          nrow(euro_precontact_only)),
  "",
  "## Step 2: Indigenous Signal Hidden by European Framework",
  "",
  sprintf("Of the %d sites classified 'no impact' by European indicators:", n_euro_no_impact),
  "",
  sprintf("- **%d sites (%.1f%%)** show Amaranthaceae/Chenopodiaceae exceedance (Set B)",
          n_misclass_indig, n_misclass_indig / n_euro_no_impact * 100),
  sprintf("- **%d sites (%.1f%%)** show disturbance exceedance (Set C)",
          n_misclass_disturb, n_misclass_disturb / n_euro_no_impact * 100),
  sprintf("- **%d sites (%.1f%%)** show signal under ANY non-European indicator",
          n_misclass_any, n_misclass_any / n_euro_no_impact * 100),
  "",
  sprintf("### Misclassification rate: %.1f%% of all ENA sites", misclass_rate_indig),
  "",
  sprintf("**%d out of %d ENA sites** have indigenous agricultural signals that the European",
          n_misclass_indig, n_total),
  "framework completely misses. These sites would be classified as 'pristine' or",
  "'no agricultural impact' in any study using only European indicator species.",
  "",
  "## Step 3: Detection Matrix",
  "",
  "### European vs Indigenous Indicators",
  "",
  "| | Indigenous: Detected | Indigenous: No impact | Indigenous: Not testable |",
  "|---|---|---|---|"
)

# Build detection matrix from cross-tab
xtab <- table(European = combined$euro_class, Indigenous = combined$indig_class)
for (rn in rownames(xtab)) {
  vals <- paste(sapply(colnames(xtab), function(cn) {
    v <- xtab[rn, cn]
    sprintf("%d", v)
  }), collapse = " | ")
  md <- c(md, sprintf("| European: %s | %s |", rn, vals))
}

md <- c(md, "",
  "### Interpretation of Detection Categories",
  "",
  "| Category | N sites | % | Meaning |",
  "|---|---|---|---|"
)

for (nm in names(cat_table)) {
  meaning <- switch(nm,
    "Both detect" = "True signal visible to both frameworks",
    "European only" = "Post-contact European agriculture (Plantago/Rumex rise)",
    "Indigenous only (MISSED by European)" = "**FALSE NEGATIVE**: EAC agriculture invisible to European framework",
    "Neither detects" = "No detectable agricultural signal under either framework"
  )
  md <- c(md, sprintf("| %s | %d | %.1f%% | %s |", nm, cat_table[nm],
                        cat_table[nm]/n_total*100, meaning))
}

md <- c(md, "",
  "## Step 4: Conceptual Detection Matrix",
  "",
  "This matrix shows the epistemic consequences of indicator choice:",
  "",
  "| Actual State | European Indicators Say | Indigenous Indicators Say | Combined Say |",
  "|---|---|---|---|",
  "| Pre-contact EAC agriculture | 'No impact' (FALSE NEGATIVE) | 'Signal detected' | 'Signal detected' |",
  "| Post-contact European agriculture | 'Signal detected' | May or may not detect | 'Signal detected' |",
  "| Post-contact forest clearance | 'Signal detected' (Rumex/Plantago) | Varies | 'Signal detected' |",
  "| No agriculture (pristine) | 'No impact' (TRUE NEGATIVE) | 'No impact' | 'No impact' |",
  "| Natural vegetation change | 'No impact' (TRUE NEGATIVE) | 'No impact' | Varies |",
  "",
  "## Step 5: Global Extrapolation",
  "",
  "### The Argument",
  "",
  sprintf("If %.1f%% of ENA sites are misclassified under European indicators alone,", misclass_rate_indig),
  "similar indicator-dependent blindness plausibly exists wherever European frameworks",
  "are applied to non-European agricultural systems.",
  "",
  "### Analogous Global Systems",
  "",
  "| Region | Agricultural System | European Indicators Applicable? | Indigenous Indicators Developed? | Estimated Blind Spot |",
  "|---|---|---|---|---|"
)

for (i in 1:nrow(global_systems)) {
  gs <- global_systems[i, ]
  md <- c(md, sprintf("| %s | %s | %s | %s | %s |",
                        gs$Region, gs$Agricultural_System,
                        gs$European_Indicator_Applicable,
                        gs$Indigenous_Indicators_Developed,
                        gs$Estimated_Blind_Spot))
}

md <- c(md, "",
  "### Cautious Statement for Paper",
  "",
  sprintf("\"At minimum, %.1f%% of ENA pollen sites (%d of %d) carry indigenous agricultural",
          misclass_rate_indig, n_misclass_indig, n_total),
  "signals that are completely invisible to the European indicator framework.",
  "This demonstrates that indicator selection is not a neutral methodological choice",
  "but an epistemic filter that determines which agricultural systems can be detected.",
  "Similar rates of misclassification are plausible wherever European pollen indicators",
  "are applied to regions with non-European domesticates — including sub-Saharan Africa",
  "(yam, millet, sorghum), South/Southeast Asia (rice in forest systems), Amazonia",
  "(manioc, managed forest), and Oceania (taro, breadfruit).\"",
  "",
  "## Methodological Notes",
  "",
  "- Baseline: >3000 cal BP (pre-maize era in ENA)",
  "- Exceedance threshold: mean + 2 SD of baseline period",
  "- Pre-contact: >458 cal BP (before 1492 CE)",
  "- Post-contact: <=458 cal BP (after 1492 CE)",
  sprintf("- Note: European Set A detects %d sites (%.1f%%), not 0%%, because broad", n_euro_signal, n_euro_signal/n_total*100),
  "  *Plantago* genus includes native North American species. Strict *P. lanceolata*-only",
  "  detection would yield a lower rate. This is itself an indicator specificity issue.",
  ""
)

report_path <- file.path(OUT_DIR, "paper8_misclassification.md")
writeLines(md, report_path)
cat(sprintf("Report saved: %s\n", report_path))

# Save combined data
write.csv(combined, file.path(OUT_DIR, "paper8_cross_classification.csv"), row.names = FALSE)
cat(sprintf("Cross-classification CSV saved\n"))

cat("\n=== ANALYSIS 9.4 COMPLETE ===\n")
