#!/usr/bin/env Rscript
# phase2_mixed_effects.R
# Proper statistical methods for combined regional H1 vs H3 analysis
# Reviewer response: mixed-effects model with region as random effect
#
# Analyses:
# 1. Pooled Mann-Whitney U test with effect size
# 2. Mixed-effects model (lme4/lmerTest) with Satterthwaite df
# 3. Meta-analytic approach (per-region effect sizes, inverse-variance weighting)
# 4. Same for total_bc

# ── 0. Setup ──────────────────────────────────────────────────────────────────

# Install packages if needed
required_pkgs <- c("lme4", "lmerTest", "effectsize", "metafor")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(lme4)
library(lmerTest)  # overrides lmer to provide Satterthwaite p-values
suppressPackageStartupMessages(library(metafor))

cat("=== Phase 2: Mixed-Effects Analysis for H1 vs H3 ===\n\n")

# ── 1. Load and harmonise data ────────────────────────────────────────────────

uk <- read.csv("/home/ayu/archeco/shared/signal_decomposition_uk_corrected.csv",
               stringsAsFactors = FALSE)
scand <- read.csv("/home/ayu/archeco/shared/signal_decomposition_scandinavia.csv",
                  stringsAsFactors = FALSE)
ce <- read.csv("/home/ayu/archeco/shared/signal_decomposition_central_europe.csv",
               stringsAsFactors = FALSE)

# Harmonise column names
# UK: pastoral_arable_pct, total_bc, classification
# Scandinavia: bc_anthro_pct, bc_total, classification
# CE: pastoral_arable_pct, total_bc, classification

harmonise <- function(df, region_name) {
  # Map to common names
  if ("bc_anthro_pct" %in% names(df)) {
    df$pastoral_arable_pct <- df$bc_anthro_pct
  }
  if ("bc_total" %in% names(df)) {
    df$total_bc <- df$bc_total
  }

  df$region <- region_name
  df[, c("sitename", "datasetid", "pastoral_arable_pct", "total_bc",
         "classification", "region")]
}

uk_h <- harmonise(uk, "UK")
scand_h <- harmonise(scand, "Scandinavia")
ce_h <- harmonise(ce, "Central_Europe")

dat <- rbind(uk_h, scand_h, ce_h)

# Filter to H1 and H3 only
dat <- dat[dat$classification %in% c("H1_anthropogenic", "H3_natural"), ]
dat$is_H1 <- ifelse(dat$classification == "H1_anthropogenic", 1, 0)
dat$classification <- factor(dat$classification, levels = c("H3_natural", "H1_anthropogenic"))
dat$region <- factor(dat$region)

cat("Sample sizes:\n")
print(table(dat$region, dat$classification))
cat("\n")

# ── Sink results to file ─────────────────────────────────────────────────────

results_file <- "/home/ayu/archeco/shared/phase2_mixed_effects.md"
sink_lines <- character(0)
add <- function(...) {
  results_file  # keep in scope
  sink_lines <<- c(sink_lines, paste0(...))
}

add("# Phase 2: Mixed-Effects Analysis — H1 vs H3 Pastoral BC%")
add("")
add("## Sample Sizes")
add("")
add("| Region | H3_natural | H1_anthropogenic |")
add("|--------|-----------|-----------------|")
tbl <- table(dat$region, dat$classification)
for (r in rownames(tbl)) {
  add(sprintf("| %s | %d | %d |", r, tbl[r, "H3_natural"], tbl[r, "H1_anthropogenic"]))
}
add(sprintf("| **Total** | **%d** | **%d** |", sum(tbl[, "H3_natural"]), sum(tbl[, "H1_anthropogenic"])))
add("")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS A: pastoral_arable_pct (H1 vs H3)
# ══════════════════════════════════════════════════════════════════════════════

add("---")
add("")
add("## A. Pastoral/Arable BC% (pastoral_arable_pct)")
add("")

# ── A1. Pooled Mann-Whitney ───────────────────────────────────────────────────

h1_vals <- dat$pastoral_arable_pct[dat$classification == "H1_anthropogenic"]
h3_vals <- dat$pastoral_arable_pct[dat$classification == "H3_natural"]

wt <- wilcox.test(h1_vals, h3_vals, exact = FALSE, correct = TRUE)

# Rank-biserial correlation (= Cliff's delta for two groups)
# r_rb = 1 - 2U / (n1 * n2)
n1 <- length(h1_vals)
n2 <- length(h3_vals)
U <- wt$statistic
r_rb <- 1 - (2 * U) / (n1 * n2)

# Also compute Cliff's delta directly
cliff_delta <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  dom <- outer(x, y, function(a, b) sign(a - b))
  sum(dom) / (nx * ny)
}
cd <- cliff_delta(h1_vals, h3_vals)

cat("A1. Pooled Mann-Whitney:\n")
cat(sprintf("  U = %.1f, p = %.6f\n", wt$statistic, wt$p.value))
cat(sprintf("  Cliff's delta = %.4f, rank-biserial r = %.4f\n", cd, -r_rb))
cat(sprintf("  H1 median = %.2f, H3 median = %.2f\n\n", median(h1_vals), median(h3_vals)))

add("### A1. Pooled Mann-Whitney U Test")
add("")
add(sprintf("- **U** = %.1f", wt$statistic))
add(sprintf("- **p-value** = %.6f", wt$p.value))
add(sprintf("- **Cliff's delta** = %.4f (probability of H1 > H3 minus probability of H3 > H1)", cd))
add(sprintf("- **Rank-biserial r** = %.4f", -r_rb))
add(sprintf("- H1 median = %.2f%%, H3 median = %.2f%%", median(h1_vals), median(h3_vals)))
add(sprintf("- H1 mean = %.2f%%, H3 mean = %.2f%%", mean(h1_vals), mean(h3_vals)))
add(sprintf("- n(H1) = %d, n(H3) = %d", n1, n2))
add("")
add("**Note**: This test ignores regional structure (treats all observations as independent).")
add("")

# ── A2. Mixed-Effects Model ───────────────────────────────────────────────────

# Check distribution
add("### A2. Mixed-Effects Model")
add("")

# Log-transform check
add("#### Distribution check")
shap_raw <- shapiro.test(dat$pastoral_arable_pct)
dat$log_pastoral <- log1p(dat$pastoral_arable_pct)  # log(1+x) for zero-safe
shap_log <- shapiro.test(dat$log_pastoral)

add(sprintf("- Raw pastoral_arable_pct: Shapiro-Wilk W = %.4f, p = %.4f", shap_raw$statistic, shap_raw$p.value))
add(sprintf("- log(1+pastoral_arable_pct): Shapiro-Wilk W = %.4f, p = %.4f", shap_log$statistic, shap_log$p.value))

# Fit both raw and log-transformed
cat("A2. Mixed-Effects Model:\n")

# Raw model
m_raw <- lmer(pastoral_arable_pct ~ classification + (1 | region), data = dat)
s_raw <- summary(m_raw)

# Log-transformed model
m_log <- lmer(log_pastoral ~ classification + (1 | region), data = dat)
s_log <- summary(m_log)

# Choose model based on residual normality
res_raw <- residuals(m_raw)
res_log <- residuals(m_log)
shap_res_raw <- shapiro.test(res_raw)
shap_res_log <- shapiro.test(res_log)

add(sprintf("- Raw model residuals: Shapiro-Wilk p = %.4f", shap_res_raw$p.value))
add(sprintf("- Log model residuals: Shapiro-Wilk p = %.4f", shap_res_log$p.value))
add("")

# Use log if it gives better residual normality, otherwise raw
if (shap_res_log$p.value > shap_res_raw$p.value) {
  best_model <- m_log
  best_summary <- s_log
  transform_used <- "log(1 + pastoral_arable_pct)"
  add("**Selected model**: log-transformed (better residual normality)")
} else {
  best_model <- m_raw
  best_summary <- s_raw
  transform_used <- "pastoral_arable_pct (untransformed)"
  add("**Selected model**: untransformed (better residual normality)")
}
add("")

# Report both models anyway
for (mod_label in c("raw", "log")) {
  if (mod_label == "raw") {
    s <- s_raw
    resp <- "pastoral_arable_pct"
  } else {
    s <- s_log
    resp <- "log(1 + pastoral_arable_pct)"
  }

  coef_tbl <- coef(s)
  fe <- coef_tbl["classificationH1_anthropogenic", ]

  add(sprintf("#### Model: %s ~ classification + (1|region)", resp))
  add("")
  add("| Parameter | Estimate | SE | t-value | df | p-value |")
  add("|-----------|----------|-----|---------|-----|---------|")
  add(sprintf("| Intercept (H3) | %.4f | %.4f | %.3f | %.1f | %.6f |",
              coef_tbl["(Intercept)", "Estimate"],
              coef_tbl["(Intercept)", "Std. Error"],
              coef_tbl["(Intercept)", "t value"],
              coef_tbl["(Intercept)", "df"],
              coef_tbl["(Intercept)", "Pr(>|t|)"]))
  add(sprintf("| H1 vs H3 | %.4f | %.4f | %.3f | %.1f | %.6f |",
              fe["Estimate"], fe["Std. Error"], fe["t value"], fe["df"], fe["Pr(>|t|)"]))
  add("")

  # Random effects
  vc <- as.data.frame(VarCorr(if (mod_label == "raw") m_raw else m_log))
  add(sprintf("- Random effect (region): SD = %.4f", vc$sdcor[vc$grp == "region"]))
  add(sprintf("- Residual SD = %.4f", vc$sdcor[vc$grp == "Residual"]))

  # ICC
  var_region <- vc$vcov[vc$grp == "region"]
  var_resid <- vc$vcov[vc$grp == "Residual"]
  icc <- var_region / (var_region + var_resid)
  add(sprintf("- ICC (region) = %.4f", icc))
  add("")

  cat(sprintf("  %s model: beta = %.4f, SE = %.4f, t = %.3f, p = %.6f\n",
              mod_label, fe["Estimate"], fe["Std. Error"], fe["t value"], fe["Pr(>|t|)"]))
}
cat("\n")

# ── A3. Meta-Analytic Approach ────────────────────────────────────────────────

add("### A3. Meta-Analytic Approach (per-region effect sizes)")
add("")

regions <- levels(dat$region)
meta_results <- data.frame(
  region = character(),
  n_h1 = integer(),
  n_h3 = integer(),
  cliff_d = numeric(),
  r_rb = numeric(),
  stringsAsFactors = FALSE
)

for (reg in regions) {
  d_reg <- dat[dat$region == reg, ]
  h1_r <- d_reg$pastoral_arable_pct[d_reg$classification == "H1_anthropogenic"]
  h3_r <- d_reg$pastoral_arable_pct[d_reg$classification == "H3_natural"]

  if (length(h1_r) < 2 || length(h3_r) < 2) {
    cat(sprintf("  Skipping %s: insufficient data (n_H1=%d, n_H3=%d)\n", reg, length(h1_r), length(h3_r)))
    next
  }

  cd_r <- cliff_delta(h1_r, h3_r)

  meta_results <- rbind(meta_results, data.frame(
    region = reg,
    n_h1 = length(h1_r),
    n_h3 = length(h3_r),
    cliff_d = cd_r,
    stringsAsFactors = FALSE
  ))
}

# SE of Cliff's delta: SE ≈ sqrt((2*(m+n) + 1) / (6*m*n)) for approximate
# More precise: Long et al. (2003) formula
meta_results$n_total <- meta_results$n_h1 + meta_results$n_h3
meta_results$se <- sqrt(
  (2 * meta_results$n_total + 1) / (6 * meta_results$n_h1 * meta_results$n_h3)
)

# Convert Cliff's delta to r for meta-analysis (Cliff's delta = 2*r_rb - 1... actually r_rb = delta)
# Use Fisher's z transformation for combining
# Actually, let's use metafor with raw effect sizes
# Cliff's delta can be treated as a standardised mean difference for combining

add("| Region | n(H1) | n(H3) | Cliff's delta | SE | 95% CI |")
add("|--------|-------|-------|--------------|------|--------|")
for (i in seq_len(nrow(meta_results))) {
  ci_lo <- meta_results$cliff_d[i] - 1.96 * meta_results$se[i]
  ci_hi <- meta_results$cliff_d[i] + 1.96 * meta_results$se[i]
  add(sprintf("| %s | %d | %d | %.4f | %.4f | [%.3f, %.3f] |",
              meta_results$region[i], meta_results$n_h1[i], meta_results$n_h3[i],
              meta_results$cliff_d[i], meta_results$se[i], ci_lo, ci_hi))
}
add("")

# Inverse-variance weighted meta-analysis
if (nrow(meta_results) >= 2) {
  # Using metafor for proper random-effects meta-analysis
  ma <- tryCatch({
    rma(yi = cliff_d, sei = se, data = meta_results, method = "REML")
  }, error = function(e) {
    # If REML fails, try fixed-effects
    rma(yi = cliff_d, sei = se, data = meta_results, method = "FE")
  })

  cat("A3. Meta-analytic combined effect:\n")
  cat(sprintf("  Combined Cliff's delta = %.4f [%.4f, %.4f]\n", ma$beta, ma$ci.lb, ma$ci.ub))
  cat(sprintf("  z = %.3f, p = %.6f\n", ma$zval, ma$pval))
  cat(sprintf("  I² = %.1f%%, Q = %.3f (p = %.4f)\n\n", ma$I2, ma$QE, ma$QEp))

  add(sprintf("**Combined effect (random-effects meta-analysis, REML)**:"))
  add(sprintf("- Combined Cliff's delta = **%.4f** [95%% CI: %.4f, %.4f]",
              as.numeric(ma$beta), ma$ci.lb, ma$ci.ub))
  add(sprintf("- z = %.3f, **p = %.6f**", ma$zval, ma$pval))
  add(sprintf("- Heterogeneity: I² = %.1f%%, Q = %.3f (p = %.4f)", ma$I2, ma$QE, ma$QEp))
  add(sprintf("- tau² = %.6f", ma$tau2))
  add("")

  # Also fixed-effects for comparison
  ma_fe <- rma(yi = cliff_d, sei = se, data = meta_results, method = "FE")
  add(sprintf("**Fixed-effects comparison**: delta = %.4f [%.4f, %.4f], z = %.3f, p = %.6f",
              as.numeric(ma_fe$beta), ma_fe$ci.lb, ma_fe$ci.ub, ma_fe$zval, ma_fe$pval))
  add("")
}

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS B: total_bc (H1 vs H3)
# ══════════════════════════════════════════════════════════════════════════════

add("---")
add("")
add("## B. Total Bray-Curtis Dissimilarity (total_bc)")
add("")

# ── B1. Pooled Mann-Whitney ───────────────────────────────────────────────────

h1_bc <- dat$total_bc[dat$classification == "H1_anthropogenic"]
h3_bc <- dat$total_bc[dat$classification == "H3_natural"]

wt_bc <- wilcox.test(h1_bc, h3_bc, exact = FALSE, correct = TRUE)
cd_bc <- cliff_delta(h1_bc, h3_bc)

add("### B1. Pooled Mann-Whitney U Test")
add("")
add(sprintf("- **U** = %.1f, **p-value** = %.6f", wt_bc$statistic, wt_bc$p.value))
add(sprintf("- **Cliff's delta** = %.4f", cd_bc))
add(sprintf("- H1 median = %.4f, H3 median = %.4f", median(h1_bc), median(h3_bc)))
add(sprintf("- H1 mean = %.4f, H3 mean = %.4f", mean(h1_bc), mean(h3_bc)))
add("")

cat("B1. Pooled Mann-Whitney (total_bc):\n")
cat(sprintf("  U = %.1f, p = %.6f, Cliff's delta = %.4f\n\n", wt_bc$statistic, wt_bc$p.value, cd_bc))

# ── B2. Mixed-Effects Model ──────────────────────────────────────────────────

add("### B2. Mixed-Effects Model")
add("")

m_bc_raw <- lmer(total_bc ~ classification + (1 | region), data = dat)
s_bc_raw <- summary(m_bc_raw)

dat$log_bc <- log(dat$total_bc)  # total_bc should always be > 0
m_bc_log <- lmer(log_bc ~ classification + (1 | region), data = dat)
s_bc_log <- summary(m_bc_log)

# Check residuals
shap_bc_raw <- shapiro.test(residuals(m_bc_raw))
shap_bc_log <- shapiro.test(residuals(m_bc_log))

add(sprintf("- Raw model residuals: Shapiro-Wilk p = %.4f", shap_bc_raw$p.value))
add(sprintf("- Log model residuals: Shapiro-Wilk p = %.4f", shap_bc_log$p.value))
add("")

for (mod_label in c("raw", "log")) {
  if (mod_label == "raw") {
    s <- s_bc_raw
    m <- m_bc_raw
    resp <- "total_bc"
  } else {
    s <- s_bc_log
    m <- m_bc_log
    resp <- "log(total_bc)"
  }

  coef_tbl <- coef(s)
  fe <- coef_tbl["classificationH1_anthropogenic", ]

  add(sprintf("#### Model: %s ~ classification + (1|region)", resp))
  add("")
  add("| Parameter | Estimate | SE | t-value | df | p-value |")
  add("|-----------|----------|-----|---------|-----|---------|")
  add(sprintf("| Intercept (H3) | %.4f | %.4f | %.3f | %.1f | %.6f |",
              coef_tbl["(Intercept)", "Estimate"],
              coef_tbl["(Intercept)", "Std. Error"],
              coef_tbl["(Intercept)", "t value"],
              coef_tbl["(Intercept)", "df"],
              coef_tbl["(Intercept)", "Pr(>|t|)"]))
  add(sprintf("| H1 vs H3 | %.4f | %.4f | %.3f | %.1f | %.6f |",
              fe["Estimate"], fe["Std. Error"], fe["t value"], fe["df"], fe["Pr(>|t|)"]))
  add("")

  vc <- as.data.frame(VarCorr(m))
  add(sprintf("- Random effect (region): SD = %.4f", vc$sdcor[vc$grp == "region"]))
  add(sprintf("- Residual SD = %.4f", vc$sdcor[vc$grp == "Residual"]))
  var_region <- vc$vcov[vc$grp == "region"]
  var_resid <- vc$vcov[vc$grp == "Residual"]
  icc <- var_region / (var_region + var_resid)
  add(sprintf("- ICC (region) = %.4f", icc))
  add("")

  cat(sprintf("  %s model (total_bc): beta = %.4f, SE = %.4f, t = %.3f, p = %.6f\n",
              mod_label, fe["Estimate"], fe["Std. Error"], fe["t value"], fe["Pr(>|t|)"]))
}

# ── B3. Meta-Analytic ─────────────────────────────────────────────────────────

add("### B3. Meta-Analytic Approach (total_bc)")
add("")

meta_bc <- data.frame(
  region = character(), n_h1 = integer(), n_h3 = integer(),
  cliff_d = numeric(), stringsAsFactors = FALSE
)

for (reg in regions) {
  d_reg <- dat[dat$region == reg, ]
  h1_r <- d_reg$total_bc[d_reg$classification == "H1_anthropogenic"]
  h3_r <- d_reg$total_bc[d_reg$classification == "H3_natural"]

  if (length(h1_r) < 2 || length(h3_r) < 2) next

  meta_bc <- rbind(meta_bc, data.frame(
    region = reg, n_h1 = length(h1_r), n_h3 = length(h3_r),
    cliff_d = cliff_delta(h1_r, h3_r), stringsAsFactors = FALSE
  ))
}

meta_bc$se <- sqrt((2 * (meta_bc$n_h1 + meta_bc$n_h3) + 1) / (6 * meta_bc$n_h1 * meta_bc$n_h3))

add("| Region | n(H1) | n(H3) | Cliff's delta | SE |")
add("|--------|-------|-------|--------------|------|")
for (i in seq_len(nrow(meta_bc))) {
  add(sprintf("| %s | %d | %d | %.4f | %.4f |",
              meta_bc$region[i], meta_bc$n_h1[i], meta_bc$n_h3[i],
              meta_bc$cliff_d[i], meta_bc$se[i]))
}
add("")

if (nrow(meta_bc) >= 2) {
  ma_bc <- tryCatch({
    rma(yi = cliff_d, sei = se, data = meta_bc, method = "REML")
  }, error = function(e) {
    rma(yi = cliff_d, sei = se, data = meta_bc, method = "FE")
  })

  add(sprintf("**Combined effect**: delta = %.4f [%.4f, %.4f], z = %.3f, p = %.6f, I² = %.1f%%",
              as.numeric(ma_bc$beta), ma_bc$ci.lb, ma_bc$ci.ub, ma_bc$zval, ma_bc$pval, ma_bc$I2))
  add("")
}

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY & RECOMMENDATION
# ══════════════════════════════════════════════════════════════════════════════

add("---")
add("")
add("## Summary: Method Comparison")
add("")
add("### Pastoral/Arable BC%")
add("")
add("| Method | Statistic | p-value | Effect size |")
add("|--------|-----------|---------|-------------|")
add(sprintf("| Pooled Mann-Whitney | U = %.1f | %.6f | Cliff's d = %.3f |",
            wt$statistic, wt$p.value, cd))

# Get lmer p-values for summary
s_raw_coef <- coef(s_raw)["classificationH1_anthropogenic", ]
s_log_coef <- coef(s_log)["classificationH1_anthropogenic", ]
add(sprintf("| LMM (raw) | t = %.3f | %.6f | beta = %.3f |",
            s_raw_coef["t value"], s_raw_coef["Pr(>|t|)"], s_raw_coef["Estimate"]))
add(sprintf("| LMM (log) | t = %.3f | %.6f | beta = %.3f |",
            s_log_coef["t value"], s_log_coef["Pr(>|t|)"], s_log_coef["Estimate"]))
if (exists("ma")) {
  add(sprintf("| Meta-analysis (RE) | z = %.3f | %.6f | delta = %.3f |",
              ma$zval, ma$pval, as.numeric(ma$beta)))
}
add("")

add("### Total BC")
add("")
add("| Method | Statistic | p-value | Effect size |")
add("|--------|-----------|---------|-------------|")
add(sprintf("| Pooled Mann-Whitney | U = %.1f | %.6f | Cliff's d = %.3f |",
            wt_bc$statistic, wt_bc$p.value, cd_bc))
s_bc_raw_coef <- coef(s_bc_raw)["classificationH1_anthropogenic", ]
s_bc_log_coef <- coef(s_bc_log)["classificationH1_anthropogenic", ]
add(sprintf("| LMM (raw) | t = %.3f | %.6f | beta = %.3f |",
            s_bc_raw_coef["t value"], s_bc_raw_coef["Pr(>|t|)"], s_bc_raw_coef["Estimate"]))
add(sprintf("| LMM (log) | t = %.3f | %.6f | beta = %.3f |",
            s_bc_log_coef["t value"], s_bc_log_coef["Pr(>|t|)"], s_bc_log_coef["Estimate"]))
if (exists("ma_bc")) {
  add(sprintf("| Meta-analysis (RE) | z = %.3f | %.6f | delta = %.3f |",
              ma_bc$zval, ma_bc$pval, as.numeric(ma_bc$beta)))
}
add("")

add("## Recommendation for Paper")
add("")
add("The **mixed-effects model (lmerTest with Satterthwaite df)** is the most")
add("appropriate method for the paper because:")
add("")
add("1. It directly addresses the reviewer's concern about regional structure")
add("2. It properly accounts for non-independence of observations within regions")
add("3. The ICC quantifies how much variance is explained by region")
add("4. Satterthwaite degrees of freedom provide accurate p-values for small")
add("   numbers of random-effect levels (k=3 regions)")
add("")
add("The meta-analytic approach provides a useful complement:")
add("- Shows consistency of the effect across regions (I² statistic)")
add("- Per-region effect sizes are interpretable for the reader")
add("")
add("Report in paper as:")
add("```")
add("We tested the difference in pastoral/arable BC% between H1 and H3 sites")
add("using a linear mixed-effects model with classification as a fixed effect")
add("and region as a random intercept (lme4/lmerTest; Satterthwaite df).")
add("```")

# Write results
writeLines(sink_lines, results_file)
cat(sprintf("\nResults written to: %s\n", results_file))
cat("Done.\n")
