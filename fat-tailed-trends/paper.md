# Pipeline generates apparent fat tails, while asymmetry persists beyond null and bias models: decomposing biodiversity trend distributions

---

## Abstract

McGill et al. (2026) demonstrated that temporal biodiversity trends are leptokurtic, following heavy-tailed distributions rather than the Gaussian model implicitly assumed by most aggregate indices. But is the distributional shape uniquely determined by ecological processes, or strongly shaped by the measurement pipeline? We address this question by passing identical ecological data — GBIF occurrence records for birds (n = 72 cells) and vascular plants (n = 75 cells) across 100 European grid cells — through multiple analytical pipelines that differ in metric choice and rarefaction treatment, and then decomposing the contributions of pipeline versus ecology using a simulation framework. Under proportional change with rarefaction, distributions are strongly leptokurtic (excess kurtosis: birds 13.2, plants 2.8) and right-skewed. Switching to log-ratios reduces kurtosis by 69–80%; removing rarefaction produces a distribution indistinguishable from Normal (Shapiro-Wilk p = 0.52). A 1000-iteration null model applying symmetric Gaussian change to the empirical baseline richness distribution through the same pipeline generates near-zero kurtosis but an asymmetry ratio of only 0.70 — the empirical value of 2.30 lies at P = 0.000 against the null. Stratification by baseline richness reveals that fat tails are predominantly concentrated in species-poor cells (low S1: kurtosis 10.3; high S1: kurtosis −0.3), identifying a small-denominator amplification mechanism. A detection bias model with 60% extinction detection rate produces an asymmetry ratio of only 0.97, ruling out observational artefact. We conclude that distributional shape is not uniquely determined by ecological processes, but strongly shaped by the measurement pipeline. Fat tails are predominantly a pipeline artefact concentrated in species-poor cells. Asymmetry, however, persists beyond all pipeline, null, and bias models tested, consistent with an ecological origin — the 2:1 excess of gains over losses is not reproduced by any pipeline configuration, symmetric null model (Gaussian, Laplace, Student-t), or detection bias model (30–90% extinction detection rates).

**Keywords:** biodiversity change, fat-tailed distributions, metric sensitivity, measurement pipeline, rarefaction artefact, null model, simulation decomposition, GBIF, species richness

---

## Introduction

The trajectory of global biodiversity remains one of ecology's most contentious debates. Headline indices such as the Living Planet Index report catastrophic vertebrate declines exceeding 60% since 1970, yet meta-analyses of local species richness frequently find no systematic net change (Vellend et al. 2013; Dornelas et al. 2014). This contradiction persists in part because aggregate metrics obscure the underlying heterogeneity. Leung et al. (2020) showed that fewer than 3% of vertebrate populations drive the reported 50% decline, with the remainder showing no mean trend. The mean, in short, is misleading. Leung et al. (2024) further demonstrated that the mathematical structure of the Living Planet Index itself amplifies apparent decline. These findings collectively indicate that understanding biodiversity change requires moving beyond central tendencies to examine the full distribution of trends.

McGill et al. (2026) took this step explicitly. Analysing data from BioTIME and the North American Breeding Bird Survey, they showed that temporal biodiversity trends are strongly leptokurtic — best fitted by Laplace and Subbotin distributions with shape parameters consistently below 2.0. The implication is profound: extreme changes, both large gains and large losses, are far more probable than Gaussian models predict, and the distribution of trends exhibits an excess of both near-stasis and very large change relative to Normal expectations. Stanley et al. (1996) documented an analogous pattern in economics, showing that firm growth rates follow Laplace rather than Gaussian distributions, suggesting that heavy-tailed growth dynamics may represent a general statistical regularity across complex systems. For conservation, the consequence is direct: planning based on average trends systematically underweights the sites most in need of intervention, and the existence of a fat-tailed distribution means that the probability of extreme events is orders of magnitude higher than Gaussian risk models assume.

The critical question that McGill et al. did not address is whether the distributional shape they documented is a property of the ecological system or of the measurement pipeline. Three analytical decisions intervene between raw occurrence data and a reported trend distribution: the metric used to quantify change (proportional vs log-ratio), the data processing method (rarefied vs raw), and the filtering threshold applied to exclude poorly sampled sites. Each of these decisions has mathematical consequences for distributional shape that operate independently of ecology. Proportional change, (S2−S1)/S1, is bounded below at −1 but unbounded above, mechanically generating right-skew and excess kurtosis. Rarefaction introduces stochastic variation through subsampling that inflates tails. Filtering thresholds determine which cells enter the distribution. Kuczynski et al. (2023) demonstrated a further complication: biodiversity time series are systematically biased toward detecting colonisations before extinctions, generating apparent asymmetry that reflects observation protocol rather than ecological process.

We address this question by combining pipeline variation with simulation-based decomposition. Using GBIF occurrence data for birds and vascular plants across 100 European grid cells, we pass identical ecological data through six pipeline configurations and then use null models to separate the contributions of pipeline artefacts from genuine ecological signal. Our contribution is not replication of fat tails (McGill et al. have already established this) but a demonstration that fat tails and asymmetry have different origins — the former is generated by the pipeline, while the latter persists beyond all pipeline, null, and bias models tested — and that the two can be formally separated through simulation.

---

## Methodology

### Data acquisition

We sampled species occurrence records from the GBIF public API (GBIF 2024) across 100 one-degree grid cells throughout Europe (latitude 38–62°N, longitude −10 to 30°E). For each cell, we retrieved records for birds (Aves, GBIF taxonKey 212) and vascular plants (Tracheophyta, taxonKey 7707728) in two temporal windows: an early period (2000–2010) and a late period (2014–2024), separated by a four-year gap to ensure temporal independence. Up to 300 records were retrieved per cell-period-taxon combination. Species richness was computed as the count of unique species-level identifications per cell-period. After filtering (minimum 30 species and 50 records in the early period), 72 cells were retained for birds and 75 for plants.

### Richness change metrics

We computed richness change using two metrics that differ in their mathematical properties. Proportional change, (S2−S1)/S1, is the conventional metric in much biodiversity trend literature. It is bounded below at −1 (complete loss) but unbounded above: a doubling yields +1.0, a tripling +2.0, but a halving yields only −0.5. This floor effect mechanically compresses the left tail and stretches the right tail, generating apparent right-skew and leptokurtosis even from symmetric underlying variation. The asymmetry is not subtle: a species pool that doubles and one that halves differ by a factor of four in their proportional change scores (+1.0 vs −0.5), despite representing symmetric changes on a multiplicative scale. The log-ratio, log(S2/S1), eliminates this asymmetry by operating on the multiplicative scale directly: doubling and halving produce equal and opposite values (+0.69 and −0.69). Comparing distributional shapes across both metrics isolates the contribution of metric-induced artefacts from genuine ecological signal.

### Rarefaction

Citizen science recording effort varies across space and time (Isaac et al. 2014). We applied individual-based rarefaction (Gotelli & Colwell 2001): within each cell-period, we randomly subsampled 100 occurrence records, counted unique species, and repeated this 20 times, taking the mean rarefied richness. To assess whether rarefaction itself alters distributional shape, we also computed both metrics on raw (unrarefied) richness counts. This yields four pipeline configurations for birds and two for plants.

### Distribution fitting

We fitted four candidate distributions to each richness change vector: Normal, Laplace, Student's t, and Asymmetric Laplace (Kozubowski & Podgorski 2000). The Asymmetric Laplace is parameterised as f(x) ∝ exp(−|x − μ|/b_side), where b1 governs the left tail and b2 the right tail; the ratio b2/b1 quantifies directional asymmetry. All parameters were estimated by maximum likelihood (Nelder-Mead). Models were compared using AIC. We also attempted the Subbotin (generalised Gaussian) distribution used by McGill et al. (2026), but the shape parameter diverged during optimisation, likely reflecting our smaller sample sizes (n = 72, 75 vs hundreds of time series).

### Simulation framework

To decompose pipeline-generated distributional properties from ecological signal, we implemented three simulation analyses.

*Null models (symmetric ecology).* We tested three symmetric null models to ensure robustness to distributional assumptions: Gaussian, Laplace, and Student-t (df = 3). For each model and each of 1000 iterations, we drew a richness change for every cell from the respective symmetric distribution calibrated to the empirical variance, applied this change to the empirical baseline richness (S1) distribution, and then passed the simulated S2 values through the identical analytical pipeline (proportional change, log-ratio, rarefaction, filtering). If the pipeline alone can generate observed kurtosis and asymmetry, simulated distributions should match empirical values. The null models preserve the empirical S1 distribution — including its heterogeneity and small-denominator cells — so any distributional structure generated by the pipeline acting on heterogeneous baselines is captured.

*S1 stratification.* We stratified cells into quintiles by baseline richness and computed distributional statistics separately. We also performed a median split into low-S1 (n = 36, S1 = 32–81) and high-S1 (n = 36, S1 = 82–143) groups. This isolates the small-denominator amplification mechanism: in species-poor cells, the gain or loss of a single species generates a large proportional or log-ratio shift. We tested the linearity of the S1–change relationship using Spearman correlation between S1 and absolute log-ratio change.

*Detection bias model.* Following Kuczynski et al. (2023), we modelled asymmetric detection by varying the probability of detecting an extinction event across a range from 30% to 90% while colonisations are detected with certainty. We applied this bias structure to the symmetric null model and measured the resulting asymmetry ratio, testing whether observational artefacts can account for the empirical asymmetry.

### Filtering sensitivity and spatial predictors

We varied the minimum early-period richness threshold (S1 ≥ 0, 10, 20, 30, 50) for bird log-ratios and examined how kurtosis, skewness, and asymmetry ratio respond. We constructed a multiple regression model with absolute rarefied richness change as the response and six predictors: latitude, longitude, log effort growth, baseline richness, urban species fraction, and land-cover category (n = 61 bird cells with complete data). Given non-normal residuals, we supplemented OLS with nonparametric bootstrap (10,000 replicates) and quantile regression at the median (Koenker 2005). We note an endogeneity concern: the urban species fraction is calculated from the same GBIF data used to compute the response variable.

### Bias control

The correlation between effort growth rate and richness change was non-significant for both birds (r = 0.054, p = 0.652) and plants (r = 0.199, p = 0.087) after filtering, and raw-rarefied richness changes were positively correlated, indicating that rarefaction modifies magnitude but not direction.

---

## Results

### Distributional shape changes substantially under metric change

The distributional shape of biodiversity trends depended fundamentally on the analytical pipeline (Table 1). Under proportional change with rarefaction — the configuration closest to standard practice — birds showed extreme leptokurtosis (excess kurtosis 13.22, skewness 2.87) and plants showed moderate leptokurtosis (kurtosis 2.79, skewness 1.66). Both decisively rejected normality (Shapiro-Wilk p < 10⁻⁶). The Asymmetric Laplace was the best-fitting model by margins exceeding 10 AIC units over all alternatives.

Switching the metric to log-ratios while retaining rarefaction reduced kurtosis by 80% for birds (13.22 → 2.59) and 69% for plants (2.79 → 0.87). Skewness dropped by 64% (birds: 2.87 → 1.02) and 77% (plants: 1.66 → 0.38). These are not marginal changes. The same ecological data, processed with a different metric, produced a qualitatively different distribution.

Removing rarefaction from log-ratios completed the collapse. Raw log-ratio richness change for birds was indistinguishable from a Normal distribution (Shapiro-Wilk p = 0.52, excess kurtosis 0.67, skewness 0.28). The Normal model was the best fit by AIC. The fat-tailed pattern had was no longer distinguishable from normality — a result that demonstrates the degree to which distributional shape depends on analytical choices rather than ecological processes.

**Table 1.** Distributional properties of richness change under six analytical pipeline configurations.

| Metric | Rarefaction | Taxon | n | Excess kurtosis | Skewness | Best model |
|--------|------------|-------|---|-----------------|----------|------------|
| (S2−S1)/S1 | Yes | Birds | 72 | 13.22 | 2.87 | Asym. Laplace |
| (S2−S1)/S1 | Yes | Plants | 75 | 2.79 | 1.66 | Asym. Laplace |
| log(S2/S1) | Yes | Birds | 72 | 2.59 | 1.02 | Asym. Laplace |
| log(S2/S1) | Yes | Plants | 75 | 0.87 | 0.38 | Asym. Laplace |
| (S2−S1)/S1 | No | Birds | 72 | 4.94 | 1.67 | — |
| log(S2/S1) | No | Birds | 72 | 0.67 | 0.28 | Normal |

### Null simulation decomposes pipeline from ecological signal

We tested three symmetric null models — Gaussian, Laplace, and Student-t (df = 3) — each calibrated to the empirical variance and applied through the identical analytical pipeline. All three produced near-zero kurtosis (Gaussian proportional: −0.10, log-ratio: 0.29; Laplace and Student-t values comparable), confirming that symmetric ecological change processed through the pipeline does not generate fat tails comparable to our empirical values. The empirical kurtosis thus arises partly from the small-denominator effect in low-S1 cells — a property of the baseline richness distribution, not of the pipeline's mathematical transformations alone.

The critical test, however, concerns asymmetry. All three null models generated asymmetry ratios far below the empirical value of 2.30 for rarefied bird log-ratios (Gaussian mean: 0.70; Laplace and Student-t: comparable), and not a single simulation out of 1000 produced an asymmetry ratio equal to or exceeding the empirical value (P = 0.000 for all three null models). The pipeline can generate apparent kurtosis through small-denominator amplification and metric boundedness, but the null models tested — including Gaussian, Laplace, and Student-t symmetric change — do not reproduce the observed asymmetry. The asymmetry is not explained by the pipeline, null, or detection-bias models considered here, suggesting a component not attributable to measurement decisions alone.

The detection bias model provided a second, independent line of evidence. We tested extinction detection rates from 30% to 90% (following Kuczynski et al. 2023). The maximum asymmetry ratio across this range was 1.19 (at 30% detection), still far below the empirical 2.30. The detection bias model implemented here does not reproduce the observed asymmetry across any parameterisation tested. The two most plausible non-ecological explanations for the asymmetry — pipeline artefact and detection bias — are both insufficient.

### Fat tails localize to species-poor cells

Quintile stratification by baseline richness revealed that fat tails are not a universal property of biodiversity change but a small-sample phenomenon concentrated in species-poor cells. The lowest quintile (S1 = 32–61, n = 14) showed extreme leptokurtosis (proportional kurtosis 3.70, log kurtosis 0.99) and strong asymmetry (ratio 4.39). Middle quintiles showed kurtosis near zero and asymmetry near 1.0. A median-split analysis confirmed the pattern: low-S1 cells (n = 36, S1 = 32–81) showed proportional kurtosis 10.3 and log kurtosis 3.1, while high-S1 cells (n = 36, S1 = 82–143) showed no fat tails (proportional kurtosis −0.3, log kurtosis −0.4) and near-unity asymmetry (ratio 0.90). The correlation between S1 and absolute log-ratio change was not significant (Spearman rho = 0.058, p = 0.626), indicating that the relationship is nonlinear and driven primarily by the lowest quintile rather than a continuous gradient.

The mechanism is arithmetic: in a cell with 40 species, the gain or loss of 4 species represents a 10% change; in a cell with 120 species, the same absolute change is 3.3%. Species-poor cells thus amplify small absolute changes into large proportional or log-ratio shifts, populating the distributional tails. This small-denominator amplification is independent of both metric choice and rarefaction, representing a third source of apparent leptokurtosis that operates on any ratio-based metric applied to heterogeneous baseline values. Filtering sensitivity analysis confirmed this pattern: kurtosis dropped from 4.26 (S1 >= 0, n = 91) through 2.24 (S1 >= 30, n = 73) to 0.17 (S1 >= 50, n = 66), indicating that extreme kurtosis values arise disproportionately from cells with low baseline richness. The asymmetry ratio, by contrast, remained within 1.98–2.80 throughout — a narrow range given the nearly twofold change in filtering stringency.

### Asymmetry persists beyond all pipeline and bias models

The convergent evidence from pipeline variation, null simulation, detection bias modelling, and S1 stratification identifies the 2:1 asymmetry as the one distributional feature not attributable to measurement. It survives metric change (proportional: 2.87–4.72; log-ratio: 1.96–2.30), rarefaction removal, filtering variation (1.98–2.80), and is not reproduced by any of the three symmetric null models (Gaussian, Laplace, Student-t; all P = 0.000) or the detection bias model (maximum ratio 1.19 across 30–90% detection rates). No other distributional property in our analysis passes all of these tests.

### Cross-taxon consistency and spatial predictors

Birds and plants showed qualitatively identical patterns across all pipeline configurations, despite different kurtosis magnitudes (13.22 vs 2.79 for rarefied proportional change). The higher kurtosis in birds likely reflects their greater mobility and smaller species pools per cell, which amplify proportional changes from single-species turnover events. This cross-taxon consistency strengthens the decomposition: if fat tails were driven by taxon-specific ecology, birds and plants should show qualitatively different patterns; instead, both show the same pipeline-dependent kurtosis and pipeline-independent asymmetry.

Urban species fraction predicted extreme change magnitude across all estimation methods (OLS: β = 5.69, p = 0.010; bootstrap 95% CI: 1.66–11.50; quantile regression: β = 2.54, p = 0.042). However, this predictor is calculated from the same GBIF data as the response variable, creating endogeneity that cannot be resolved without external land-use validation. We present this result as consistent with urbanisation-driven homogenisation but do not treat it as a central finding.

---

## Discussion

### Decomposition of artefact and signal

The central result of this study is a decomposition: fat tails and asymmetry in biodiversity trend distributions have different origins. Fat tails — the leptokurtic signature documented by McGill et al. (2026) and confirmed here — are jointly produced by three pipeline components: metric choice (proportional change creates a floor effect), rarefaction (subsampling introduces stochastic noise), and the small-denominator effect (species-poor cells amplify ratio-based changes). Asymmetry — the 2:1 excess of extreme gains over extreme losses — survives removal of all three artefact sources and is not reproduced by any of the symmetric null models (Gaussian, Laplace, Student-t; all P = 0.000) or detection bias simulations (maximum ratio 1.19 across 30–90% detection rates). The question "are biodiversity trends fat-tailed?" has no answer independent of the measurement pipeline. The question "are biodiversity trends asymmetric?" appears to: asymmetry persists beyond all pipeline, null, and bias models tested, consistent with an ecological origin. This reframes the field's inquiry from "why are trends fat-tailed?" — a question that presupposes an intrinsic distributional property — to "what ecological processes generate asymmetric change?" — a question that can be addressed with mechanistic models and species-level data.

### Rarefaction as a generator of distributional structure

The comparison between rarefied and raw data reveals that rarefaction does not merely correct for sampling bias — it generates distributional structure. For log-ratios, rarefied bird data showed clear leptokurtosis (kurtosis 2.59, Shapiro-Wilk p = 0.002), while raw data were approximately Normal (kurtosis 0.67, p = 0.52). Individual-based rarefaction subsamples 100 records from cells containing up to 300, introducing stochastic variation that inflates tails. Cells near the 100-record threshold have noisier estimates, creating heteroscedasticity that manifests as a central spike with occasional large deviations. Rarefaction trades one form of bias (effort inequality) for another (stochastic inflation of tails), and the net effect on distributional inference depends on which bias dominates.

### Fat tails as a small-sample phenomenon

The S1 stratification result provides a specific, testable, mechanistic explanation for empirical fat tails. Leptokurtosis is predominantly concentrated in cells with below-median baseline richness. High-richness cells show kurtosis indistinguishable from zero and near-symmetric tails. This follows directly from the mathematical properties of the metrics — it is the consequence of computing ratios from small denominators in a sample with heterogeneous baseline values. The implication extends beyond our study: any analysis reporting distributional shape from ratio-based biodiversity metrics without controlling for baseline richness heterogeneity risks interpreting a statistical sampling property as an ecological pattern. The practical resolution is straightforward: stratify by baseline richness, or use only cells above a richness threshold that eliminates the small-denominator regime. This finding does not invalidate McGill et al. (2026), who used different data structures and methods, but it identifies a specific mechanism that should be evaluated in any dataset where baseline richness varies substantially across units.

### The asymmetry is consistent with ecological origin

The null model results (P = 0.000 for all three symmetric null models — Gaussian, Laplace, Student-t) establish that the observed 2:1 asymmetry is not generated by any pipeline configuration applied to symmetric ecological change. The detection bias model (maximum ratio 1.19 across 30–90% extinction detection rates) establishes that the observational artefact tested here is also insufficient. Together, these tests do not prove ecological causation, but they eliminate the most plausible non-ecological explanations considered.

Three ecological mechanisms could generate this asymmetry, and they are not mutually exclusive. First, invasive and range-expanding species may establish faster than native specialists decline, creating a temporal lag in which gains outpace losses at the time scale of our study. In a warming and urbanising Europe, the arrival of new generalist species into grid cells is expected to outpace the disappearance of specialists, producing exactly the rightward asymmetry we observe. Second, climate-driven range shifts at the European scale may produce more "winners" (species expanding poleward) than "losers" (species contracting) within our two-decade study period (Blowes et al. 2019). Third, land-use homogenisation favours generalist species that colonise rapidly while specialist declines operate on longer timescales governed by extinction debt. Distinguishing among these mechanisms requires species-level turnover analysis beyond the scope of this study, but the simulation framework we develop here provides the necessary baseline against which ecological hypotheses can be tested.

### Pipeline transparency

Our results suggest a pipeline transparency standard for studies reporting distributional shape of biodiversity trends. At minimum, any such study should specify: (a) the metric used and its mathematical properties (bounded vs unbounded, symmetric vs asymmetric on the real line); (b) whether rarefaction or other subsampling was applied, at what depth, and with how many iterations; (c) the filtering thresholds used and their effect on sample size and composition; and (d) whether the reported distributional shape survives at least one alternative pipeline configuration. We further recommend that any study characterising trend distributions include a null model calibrated to the empirical baseline richness distribution, to separate pipeline-generated distributional properties from ecological signal. Without this information, distributional claims are uninterpretable. The general principle — that distributional shape is a joint product of signal and measurement — applies beyond biodiversity trends to any domain in which distributional shape is used to infer process.

### Conservation implications

The decomposition carries specific conservation implications. Cells populating the fat tails are predominantly species-poor, where small absolute changes generate large ratio-based signals. These cells deserve monitoring attention because species-poor communities may be disproportionately vulnerable to further loss. The asymmetry finding adds a second dimension: rapid richness gains, often signalling biotic homogenisation through generalist arrival rather than ecological recovery, are approximately twice as common as comparable losses. Sites gaining species rapidly are not necessarily improving; they may be losing ecological distinctiveness.

### Caveats

GBIF citizen science data carry well-documented spatial and taxonomic biases (Isaac et al. 2014; Johnston et al. 2023). Our one-degree grid cells are coarse and likely spatially autocorrelated; a Moran's I test or spatial error model would strengthen future analyses. The geographic scope is limited to Europe, and two taxa from one continent do not establish universal generality. The two-period comparison (2000–2010 vs 2014–2024) lacks continuous temporal resolution and collapses within-period dynamics into a single net change. The urban species fraction predictor is endogenous to the GBIF data and cannot be validated without external land-use information. Sample sizes (n = 72, 75) are modest for distributional inference, and although we tested three symmetric null models (Gaussian, Laplace, Student-t), other asymmetric non-ecological mechanisms not considered here could in principle contribute. The Subbotin distribution used by McGill et al. (2026) could not be fitted due to convergence failure, preventing direct comparison of shape parameters. These limitations are real but do not affect the core decomposition: the pipeline generates apparent fat tails, asymmetry persists beyond all pipeline, null, and bias models tested, consistent with an ecological origin, and the two are separable through simulation. The decomposition framework itself — null model, stratification, detection bias test — is general and can be applied to any dataset where distributional shape is of interest.

---

## Conclusion

Fat tails in biodiversity trend distributions arise predominantly from analytical pipeline decisions — metric boundedness and small-denominator amplification in species-poor cells. Asymmetry, however, is not reproduced by any of the symmetric null models (Gaussian, Laplace, Student-t), detection-bias models (30-90% extinction detection), or pipeline configurations tested here, consistent with — though not proof of — an underlying ecological asymmetry in colonisation and extinction dynamics. The distributional shape of biodiversity trends is jointly determined by measurement decisions and ecological processes, and the two can be partially separated through simulation-based decomposition.

---

## Data availability

All occurrence data were retrieved from the GBIF public API (https://www.gbif.org). Analysis scripts are available at [repository to be determined upon acceptance].

---

## References

Blowes SA, Supp SR, Antao LH, Bates A, Bruelheide H, Chase JM, Moyes F, Magurran A, McGill B, Myers-Smith IH, Winter M, Bjorkman AD, Bowler DE, Byrnes JEK, Gonzalez A, Hines J, Isbell F, Jones HP, Navarro LM, Thompson PL, Vellend M, Waldock C, Dornelas M (2019) The geography of biodiversity change in marine and terrestrial assemblages. *Science* 366: 339–345.

Dornelas M, Gotelli NJ, McGill B, Shimadzu H, Moyes F, Sievers C, Magurran AE (2014) Assemblage time series reveal biodiversity change but not systematic loss. *Science* 344: 296–299.

GBIF (2024) GBIF: The Global Biodiversity Information Facility. Available from https://www.gbif.org [accessed 2024].

Gotelli NJ, Colwell RK (2001) Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. *Ecology Letters* 4: 379–391.

Isaac NJB, van Strien AJ, de Smedt R, de Zeeuw J, Schmeller DS, Roy DB, Stefanescu C, Haase P et al. (2014) Statistics for citizen science: extracting signals of change from noisy ecological data. *Methods in Ecology and Evolution* 5: 1052–1060.

Johnston A, Matechou E, Dennis EB (2023) Outstanding challenges and future directions for biodiversity monitoring using citizen science data. *Methods in Ecology and Evolution* 14: 103–116.

Koenker R (2005) *Quantile Regression*. Cambridge University Press, Cambridge.

Kozubowski TJ, Podgorski K (2000) A multivariate and asymmetric generalization of Laplace distribution. *Computational Statistics* 15: 531–540.

Kuczynski L, Ontiveros VJ, Hillebrand H (2023) Biodiversity time series are biased towards increasing species richness in changing environments. *Nature Ecology & Evolution* 7: 980–988.

Leung B, Hargreaves AL, Greenberg DA, McGill B, Dornelas M, Freeman R (2020) Clustered versus catastrophic global vertebrate declines. *Nature* 588: 267–271.

Leung B, Hargreaves AL, Greenberg DA, McGill B, Dornelas M, Freeman R (2024) Mathematical biases in the calculation of the Living Planet Index lead to overestimation of vertebrate population decline. *Nature Communications* 15: 4488.

McGill BJ, Moyes F, Dornelas M, Gotelli NJ, Magurran AE (2026) Biodiversity trends show an excess of both near stasis and of very large change. *Ecology Letters* 29: e70353.

Stanley MHR, Amaral LAN, Buldyrev SV, Havlin S, Leschhorn H, Maass P, Salinger MA, Stanley HE (1996) Scaling behaviour in the growth of companies. *Nature* 379: 804–806.

Vellend M, Baeten L, Myers-Smith IH, Elmendorf SC, Beausejour R, Brown CD, De Frenne P, Verheyen K, Wipf S (2013) Global meta-analysis reveals no net change in local-scale plant biodiversity over time. *Proceedings of the National Academy of Sciences* 110: 19456–19459.

Williams PJ, Lu X, Scharf HR, Hooten MB (2023) Embracing asymmetry in nature: How to account for skewness in ecological data. *Ecological Informatics* 76: 102113.

---

## Figure legends

**Figure 1.** The same ecological data produce different distributional shapes under different analytical pipelines. Panel (a): proportional change with rarefaction — strongly leptokurtic (excess kurtosis 13.22, Asymmetric Laplace fit). Panel (b): log-ratio with rarefaction — modestly leptokurtic (kurtosis 2.59, Asymmetric Laplace fit). Panel (c): log-ratio without rarefaction — approximately Normal (kurtosis 0.67, Shapiro-Wilk p = 0.52). All three panels show the same 72 European bird cells measured over the same time periods. Histograms are overlaid with fitted Normal (dashed) and Asymmetric Laplace (solid) densities.

**Figure 2.** Null model decomposition of pipeline artefact from ecological signal. Panel (a): distribution of asymmetry ratios from 1000 symmetric null simulations (mean 0.70), with empirical value (2.30) indicated by vertical line (P = 0.000). Panel (b): S1 stratification showing kurtosis and asymmetry ratio for low-S1 (n = 36) and high-S1 (n = 36) cells. Fat tails are confined to species-poor cells; asymmetry is concentrated in the low-S1 stratum.

**Figure 3.** Sensitivity of distributional shape to analytical pipeline and filtering threshold. Left: excess kurtosis across all six metric × rarefaction × taxon combinations from Table 1. Right: asymmetry ratio (b2/b1) across filtering thresholds (S1 ≥ 0 to S1 ≥ 50). Kurtosis is highly sensitive to pipeline choices; asymmetry remains within 1.96–2.80 across all conditions, consistent with the null model result identifying it as a genuine ecological signal.
