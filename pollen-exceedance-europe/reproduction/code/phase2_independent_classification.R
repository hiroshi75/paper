#!/usr/bin/env Rscript
# Phase 2: Independent (Time-Based) Classification to Break Circularity
#
# Problem: Original classification uses composition (pastoral/arable %) to define
# H1 vs H3, then tests whether H1 differs in composition — tautological.
#
# Solution: Classify sites by TIMING (proximity to regional Neolithic arrival),
# then test whether time-classified groups differ in COMPOSITION.
# This completely breaks the circularity.
#
# Neolithic arrival dates (cal BP):
#   UK = 6000, Scandinavia = 5800, Central Europe = 7500
# H1_time: signal_onset_age within ±1000 yr of Neolithic arrival
# H3_time: all other sites

library(dplyr)
library(tidyr)

# ── 1. Load and harmonise data ───────────────────────────────────────────────

uk <- read.csv("/home/ayu/archeco/shared/signal_decomposition_uk_corrected.csv",
               stringsAsFactors = FALSE) %>%
  mutate(region = "UK",
         pastoral_pct = pastoral_arable_pct,
         climate_pct  = climate_tree_pct,
         succession_pct = succession_tree_pct,
         bc_mag = total_bc,
         orig_class = classification)

scand <- read.csv("/home/ayu/archeco/shared/signal_decomposition_scandinavia.csv",
                  stringsAsFactors = FALSE) %>%
  mutate(region = "Scandinavia",
         pastoral_pct = bc_anthro_pct,
         climate_pct  = bc_climate_tree_pct,
         succession_pct = bc_succession_tree_pct,
         bc_mag = bc_total,
         orig_class = classification)

ce <- read.csv("/home/ayu/archeco/shared/signal_decomposition_central_europe.csv",
               stringsAsFactors = FALSE) %>%
  mutate(region = "CE",
         pastoral_pct = pastoral_arable_pct,
         climate_pct  = climate_tree_pct,
         succession_pct = succession_tree_pct,
         bc_mag = total_bc,
         orig_class = classification)

# Common columns
cols <- c("datasetid", "sitename", "signal_onset_age", "pastoral_pct",
          "climate_pct", "succession_pct", "bc_mag", "orig_class", "region")

all_data <- bind_rows(
  uk[, cols],
  scand[, cols],
  ce[, cols]
)

cat("Total sites loaded:", nrow(all_data), "\n")
cat("  UK:", sum(all_data$region == "UK"),
    " Scand:", sum(all_data$region == "Scandinavia"),
    " CE:", sum(all_data$region == "CE"), "\n\n")

# ── 2. Time-based classification ─────────────────────────────────────────────

neolithic <- c(UK = 6000, Scandinavia = 5800, CE = 7500)
window <- 1000  # ±1000 yr

all_data <- all_data %>%
  mutate(
    neo_arrival = neolithic[region],
    onset_offset = abs(signal_onset_age - neo_arrival),
    time_class = ifelse(onset_offset <= window, "H1_time", "H3_time")
  )

cat("=== TIME-BASED CLASSIFICATION ===\n")
cat("Window: Neolithic arrival ±", window, "years\n\n")

# Summary by region
class_summary <- all_data %>%
  group_by(region, time_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = time_class, values_from = n, values_fill = 0)
print(as.data.frame(class_summary))
cat("\n")

# Cross-tabulation: time_class vs orig_class
cat("Cross-tabulation: time_class vs original classification\n")
print(table(all_data$time_class, all_data$orig_class))
cat("\n")

# ── 3. Mann-Whitney tests: H1_time vs H3_time ────────────────────────────────

run_mw_test <- function(data, var_name, label) {
  h1 <- data[[var_name]][data$time_class == "H1_time"]
  h3 <- data[[var_name]][data$time_class == "H3_time"]

  if (length(h1) < 3 | length(h3) < 3) {
    return(data.frame(variable = label, n_H1 = length(h1), n_H3 = length(h3),
                      median_H1 = median(h1, na.rm = TRUE),
                      median_H3 = median(h3, na.rm = TRUE),
                      W = NA, p_value = NA, effect_r = NA,
                      stringsAsFactors = FALSE))
  }

  test <- wilcox.test(h1, h3, exact = FALSE)
  # Effect size r = Z / sqrt(N)
  z <- qnorm(test$p.value / 2)
  n_total <- length(h1) + length(h3)
  r <- abs(z) / sqrt(n_total)

  data.frame(
    variable = label,
    n_H1 = length(h1), n_H3 = length(h3),
    median_H1 = round(median(h1, na.rm = TRUE), 2),
    median_H3 = round(median(h3, na.rm = TRUE), 2),
    W = test$statistic,
    p_value = round(test$p.value, 4),
    effect_r = round(r, 3),
    stringsAsFactors = FALSE
  )
}

variables <- c("pastoral_pct", "climate_pct", "succession_pct", "bc_mag")
labels    <- c("Pastoral/Arable BC%", "Climate-tree BC%", "Succession-tree BC%", "Total BC magnitude")

# --- 3a. Pooled (all regions) ---
cat("=== MANN-WHITNEY TESTS: H1_time vs H3_time (POOLED) ===\n")
pooled_results <- do.call(rbind, mapply(run_mw_test,
                                         var_name = variables,
                                         label = labels,
                                         MoreArgs = list(data = all_data),
                                         SIMPLIFY = FALSE))
print(pooled_results, row.names = FALSE)
cat("\n")

# --- 3b. By region ---
cat("=== MANN-WHITNEY TESTS BY REGION ===\n")
for (reg in c("UK", "Scandinavia", "CE")) {
  cat("\n--- ", reg, " ---\n")
  reg_data <- all_data %>% filter(region == reg)
  cat("H1_time:", sum(reg_data$time_class == "H1_time"),
      " H3_time:", sum(reg_data$time_class == "H3_time"), "\n")
  reg_results <- do.call(rbind, mapply(run_mw_test,
                                        var_name = variables,
                                        label = labels,
                                        MoreArgs = list(data = reg_data),
                                        SIMPLIFY = FALSE))
  print(reg_results, row.names = FALSE)
}
cat("\n")

# ── 4. Archaeological distance gradient (Spearman) ───────────────────────────

cat("=== ARCHAEOLOGICAL DISTANCE GRADIENT (Spearman) ===\n")
cat("x = |signal_onset_age - neolithic_arrival|, y = pastoral_pct\n\n")

# Pooled
sp_pooled <- cor.test(all_data$onset_offset, all_data$pastoral_pct,
                       method = "spearman", exact = FALSE)
cat("POOLED: rho =", round(sp_pooled$estimate, 3),
    ", p =", format.pval(sp_pooled$p.value, digits = 4),
    ", n =", nrow(all_data), "\n")

# By region
for (reg in c("UK", "Scandinavia", "CE")) {
  rd <- all_data %>% filter(region == reg)
  sp <- cor.test(rd$onset_offset, rd$pastoral_pct,
                  method = "spearman", exact = FALSE)
  cat(reg, ": rho =", round(sp$estimate, 3),
      ", p =", format.pval(sp$p.value, digits = 4),
      ", n =", nrow(rd), "\n")
}
cat("\n")

# Directionality check: is the correlation NEGATIVE (closer = more pastoral)?
cat("Note: NEGATIVE rho means sites closer to Neolithic arrival have HIGHER pastoral %\n")
cat("(smaller offset → higher pastoral_pct)\n\n")

# ── 5. Sensitivity: varying the time window ──────────────────────────────────

cat("=== SENSITIVITY: VARYING TIME WINDOW ===\n")
cat("Testing windows from ±500 to ±2000 yr\n\n")

windows <- c(500, 750, 1000, 1250, 1500, 2000)
sens_results <- data.frame()

for (w in windows) {
  tmp <- all_data %>%
    mutate(time_class_w = ifelse(onset_offset <= w, "H1_time", "H3_time"))

  h1 <- tmp$pastoral_pct[tmp$time_class_w == "H1_time"]
  h3 <- tmp$pastoral_pct[tmp$time_class_w == "H3_time"]

  if (length(h1) >= 3 & length(h3) >= 3) {
    test <- wilcox.test(h1, h3, exact = FALSE)
    z <- qnorm(test$p.value / 2)
    r <- abs(z) / sqrt(length(h1) + length(h3))
    sens_results <- rbind(sens_results, data.frame(
      window = w, n_H1 = length(h1), n_H3 = length(h3),
      med_H1 = round(median(h1, na.rm = TRUE), 2),
      med_H3 = round(median(h3, na.rm = TRUE), 2),
      p_value = round(test$p.value, 4),
      effect_r = round(r, 3)
    ))
  }
}
print(sens_results, row.names = FALSE)
cat("\n")

# ── 6. Agreement between time-based and original classification ──────────────

cat("=== AGREEMENT: TIME-BASED vs ORIGINAL CLASSIFICATION ===\n")
agreement <- all_data %>%
  mutate(
    orig_H1 = grepl("H1", orig_class),
    time_H1 = time_class == "H1_time",
    agree = orig_H1 == time_H1
  )
cat("Agreement rate:", round(mean(agreement$agree) * 100, 1), "%\n")
cat("Cohen's kappa (if both binary):\n")

# Simple kappa
a <- sum(agreement$orig_H1 & agreement$time_H1)
b <- sum(agreement$orig_H1 & !agreement$time_H1)
c_val <- sum(!agreement$orig_H1 & agreement$time_H1)
d <- sum(!agreement$orig_H1 & !agreement$time_H1)
n <- nrow(agreement)
po <- (a + d) / n
pe <- ((a + b) * (a + c_val) + (c_val + d) * (b + d)) / n^2
kappa <- (po - pe) / (1 - pe)
cat("  kappa =", round(kappa, 3), "\n")
cat("  Observed agreement (po) =", round(po, 3), "\n")
cat("  Expected agreement (pe) =", round(pe, 3), "\n\n")

# Confusion matrix
cat("Confusion matrix:\n")
cat("                    Time-H1  Time-H3\n")
cat(sprintf("  Orig-H1 (anthro)  %4d     %4d\n", a, b))
cat(sprintf("  Orig-H3 (natural) %4d     %4d\n", c_val, d))
cat("\n")

# ── 7. Onset age distributions ───────────────────────────────────────────────

cat("=== ONSET AGE DISTRIBUTIONS BY REGION ===\n")
for (reg in c("UK", "Scandinavia", "CE")) {
  rd <- all_data %>% filter(region == reg)
  cat(reg, "(Neolithic =", neolithic[reg], "BP):\n")
  cat("  Mean onset:", round(mean(rd$signal_onset_age, na.rm = TRUE)),
      "  Median:", round(median(rd$signal_onset_age, na.rm = TRUE)),
      "  Range:", range(rd$signal_onset_age, na.rm = TRUE), "\n")
}
cat("\n")

cat("=== ANALYSIS COMPLETE ===\n")
