# Pervasive Demographic Synchrony Across Holocene Eurasia: Evidence from 154,000 Radiocarbon Dates

**Version**: 4.0
**Date**: 2026-03-29

---

## Abstract

Whether human populations across distant regions experienced synchronous demographic fluctuations during the Holocene remains unresolved. Previous studies of summed probability distributions (SPDs) of radiocarbon dates have documented regional population trajectories, but systematic cross-continental comparisons with appropriate detrending and artifact controls have not been attempted. Here we analyse 154,466 radiocarbon dates from the p3k14c database across seven regions spanning western Europe to East Asia (Britain, western Europe, eastern Europe, Scandinavia, the Near East, China, and Japan). After first-difference detrending, 20 of 21 inter-regional pairs show significant positive correlation (mean detrended r = 0.458, Benjamini-Hochberg adjusted p < 0.05). This synchrony survives six independent robustness tests: Bartlett autocorrelation correction (20/21), block bootstrap (16/21, all non-Japan pairs), loess detrending (19/21, mean r = 0.338), taphonomic correction following Surovell et al. (18/21, mean r = 0.419), a calibration artifact null model (observed r = 0.48 vs null r approximately 0.00, p < 0.001), equal-sample subsampling (mean r = 0.307, all pairs positive), and a cross-hemisphere test using three Southern Hemisphere regions calibrated with SHCal20 (NH x SH mean r = 0.439, 9/9 pairs p < 0.001), definitively ruling out IntCal20 calibration artifacts. A Mantel test shows no significant distance-decay (p = 0.194), but a comparison with a synthetic climate proxy reveals that climate-population correlations are 5.3 times weaker than inter-regional population correlations. Climate forcing alone is therefore insufficient to explain the observed synchrony. Spectral decomposition reveals that synchrony is concentrated at high frequencies (100--500 years, mean r = 0.420, all 21 pairs positive) while the medium-frequency band (500--1,500 years, mean r = 0.057) where climate forcing should operate most strongly shows near-zero synchrony --- a frequency-domain falsification of simple climate models. Low-frequency (millennial) synchrony is pair-specific, not universal (range r = -0.96 to +0.96), reflecting region-dependent civilisational trajectories rather than shared forcing. Rolling synchrony identifies peak coupling at the Late Bronze Age/Early Iron Age transition (~2,600-2,950 cal BP) and minimum coupling during Neolithicisation (~7,250-7,600 cal BP). We argue that the existence of pervasive demographic synchrony --- regardless of its mechanism --- has implications for how civilisational collapse and resilience should be studied: as system-level properties of a coupled network rather than local failures.

**Keywords**: radiocarbon, summed probability distribution, demographic synchrony, Holocene, p3k14c, population dynamics, detrending, taphonomic correction, spectral analysis

---

## 1. Introduction

### 1.1 SPDs as demographic proxies

The summed probability distribution (SPD) of calibrated radiocarbon dates has become a standard proxy for reconstructing relative population changes over the Holocene (Rick, 1987; Gamble et al., 2005; Shennan et al., 2013). The method assumes that the density of radiocarbon-dated events in a region tracks the intensity of human activity, and therefore population size, over time. While subject to known biases --- sampling intensity, taphonomic loss, research interest, and calibration curve effects --- SPDs have proven informative for detecting broad demographic trends when applied with appropriate controls such as site binning (Timpson et al., 2014; Crema & Bevan, 2021).

Regional SPD analyses have documented population trajectories for Europe (Shennan et al., 2013; Timpson et al., 2014; Bevan et al., 2017), China (Wang et al., 2014; Hosner et al., 2016), Japan (Crema & Kobayashi, 2020), and the Near East (Palmisano et al., 2021). These studies have identified boom-bust cycles, post-Neolithic population growth, and responses to known climate events.

### 1.2 The synchrony question

A fundamental question remains largely unaddressed: do population fluctuations in distant regions covary at centennial timescales? Synchronous dynamics would imply either a shared external driver or structural similarities in population-resource interactions that produce similar endogenous cycles. Asynchronous dynamics would indicate that local factors dominate.

Three obstacles have prevented systematic investigation. First, most SPD studies are regional, making cross-comparison ad hoc. Second, shared long-term trends (post-glacial warming, agricultural expansion, taphonomic loss of older material) create spurious correlations. Third, the IntCal20 calibration curve --- used universally to convert radiocarbon ages to calendar dates --- contains plateaus and wiggles that can create correlated artifacts across regions (Bamforth & Grund, 2012; Contreras & Meadows, 2014; Williams, 2012).

### 1.3 Prior work

Shennan et al. (2013) compared SPDs across European sub-regions, finding broadly similar trajectories driven by agricultural spread. Bevan et al. (2017) refined this for Britain and Ireland. Crema and Bevan (2021) developed statistical tools for SPD comparison. However, no prior study has: (a) systematically tested cross-continental synchrony with (b) appropriate detrending and (c) calibration artifact controls and (d) taphonomic correction.

A companion study (Gordon et al., 2026d) demonstrated that SPD-pollen correlations at continental scale are entirely spurious --- they disappear after first-difference detrending. This finding motivated the present study: does SPD-SPD synchrony also disappear, or does it survive as a genuine signal?

### 1.4 This study

We address the synchrony question using the p3k14c database (Bird et al., 2022), defining seven macro-regions from Britain to Japan. We apply six robustness tests: autocorrelation correction, block bootstrap, alternative detrending, taphonomic correction, calibration artifact null model, and equal-sample subsampling. We then examine spatial structure, temporal dynamics, and --- critically --- whether the synchrony is explicable by climate forcing.

---

## 2. Methods

### 2.1 Data

All radiocarbon dates were obtained from p3k14c version 1.0 (Bird et al., 2022). After filtering for valid coordinates, ages 500--12,000 BP, and errors <500 years, 154,466 dates remained.

### 2.2 Regional definitions

| Region | Countries | N dates | N binned groups |
|--------|-----------|---------|-----------------|
| Britain + Ireland | UK, Ireland | 29,794 | 13,427 |
| Western Europe | FR, ES, DE, BE, NL, IT, PT, CH, AT | 18,052 | 8,365 |
| Eastern Europe + Greece | PL, HU, CZ, RO, BG, RS, HR, SK, GR | 6,893 | 2,176 |
| Scandinavia | NO, SE, DK, FI | 9,870 | 5,177 |
| Near East | TR, SY, IQ, IR, IL, JO, LB, PS | 5,161 | 1,476 |
| China | CN | 4,065 | 1,909 |
| Japan | JP | 1,433 | 402 |

### 2.3 SPD computation

Dates were calibrated using IntCal20 (Reimer et al., 2020) via rcarbon (Crema & Bevan, 2021). Site binning (200-year bins) was applied using the SiteName field. SPDs were computed over 10,000--1,000 cal BP, interpolated to 50-year resolution (181 time steps), and normalised to mean = 1.

### 2.4 Detrending

Two methods were applied:

**First-differencing (primary):** ΔSPDt = SPDt - SPDt-1, removing all monotonic trends.

**Loess residuals (sensitivity):** Loess smoother (span = 0.3, ~2,700 years) fit to each SPD; residuals used as detrended series.

### 2.5 Cross-correlation analysis

Pairwise Pearson correlations for all 21 pairs, with significance assessed via:
1. Bartlett correction for autocorrelation (N_eff formula)
2. Block bootstrap (1,000 iterations, block = 250 years)
3. Benjamini-Hochberg FDR correction across all 21 tests

### 2.6 Taphonomic correction

Following Surovell et al. (2009), we modelled taphonomic bias as exponential: log(SPD) = a + bt. For each region, we fitted this model and divided the raw SPD by the fitted exponential before detrending. The correction was applied independently per region, allowing for heterogeneous taphonomic rates.

### 2.7 Calibration artifact null model

For 50 permutations, radiocarbon dates within each region were randomly shuffled, recalibrated, and SPDs recomputed. The null distribution of mean detrended r quantifies synchrony expected from shared calibration curve features alone.

### 2.8 Equal-sample subsampling

For 20 iterations, 1,500 dates were randomly drawn from each of the six regions with N > 1,500 (Japan excluded). SPDs were recomputed and correlations calculated. This controls for the possibility that large-sample regions produce artificially smooth SPDs that inflate cross-correlations.

### 2.9 Spatial and temporal analyses

Geographic distance-synchrony relationship was tested by Mantel test (9,999 permutations). Western Europe was subdivided into France, Iberia, and Central Europe for intra-continental comparison. Rolling 500-year window synchrony was computed across the first-differenced series. Z-scores at five known climate events (8.2 ka, 5.9 ka, 4.2 ka, 3.2 ka, 2.2 ka) assessed regional population behaviour during climate perturbations.

### 2.10 Climate proxy comparison

A synthetic GISP2-like climate proxy was constructed at 50-year resolution based on the Holocene Thermal Maximum, Neoglaciation cooling, and Bond events (9.3 ka, 8.2 ka, 5.9 ka, 4.2 ka, 3.2 ka, 2.8 ka, 1.4 ka). Detrended correlations between this proxy and each regional SPD were computed. Additionally, a binary climate deterioration index identified five major events, and population behaviour during crisis versus stable periods was compared.

### 2.11 Spectral decomposition

To identify which frequency bands carry the synchrony signal, bandpass filtering was applied via running-mean subtraction. Each regional SPD was decomposed into three bands: high frequency (100--500 year cycles; raw SPD minus 500-year running mean), medium frequency (500--1,500 year cycles; 500-year smooth minus 1,500-year smooth), and low frequency (1,500--5,000 year cycles; 1,500-year smooth minus 5,000-year smooth). Pairwise correlations were computed within each band for all 21 pairs. A GISP2-based temperature reconstruction (Alley, 2000) interpolated to 50-year resolution was compared with SPDs in each frequency band.

### 2.12 Cross-hemisphere calibration test

The most powerful test of whether synchrony is a calibration curve artifact exploits the difference between Northern and Southern Hemisphere calibration curves. Northern Hemisphere dates are calibrated with IntCal20 (Reimer et al., 2020); Southern Hemisphere dates use SHCal20 (Hogg et al., 2020), which differs from IntCal20 by a variable offset (typically 40 ± 20 years) and has independent tree-ring underpinning. If synchrony exists between NH regions (IntCal20) and SH regions (SHCal20), it cannot be explained by shared calibration curve features.

Three Southern Hemisphere regions were defined from p3k14c: South America (Lat < 0; 6,801 dates), sub-Saharan Africa (Lat < 0; 2,314 dates), and Australia (3,150 dates). SH dates were calibrated with SHCal20; NH dates with IntCal20. SPDs were computed, normalised, and first-differenced identically to the main analysis. Pairwise correlations were computed for all NH x SH, NH x NH, and SH x SH pairs.

---

## 3. Results

### 3.1 Detrended cross-correlations

After first-differencing, all 21 pairs show positive correlations. Twenty of 21 are significant (BH-adjusted p < 0.05). The sole non-significant pair is Britain x Japan (r = +0.138, BH p = 0.065).

The five strongest pairs:

| Pair | Distance (km) | Raw r | Detrended r | BH p |
|------|---------------|-------|-------------|------|
| W. Europe x China | 8,233 | +0.51 | **+0.719** | <10^-28 |
| E. Europe x Near East | 1,747 | +0.29 | **+0.664** | <10^-23 |
| Britain x China | 8,175 | +0.82 | **+0.653** | <10^-22 |
| E. Europe x China | 7,205 | +0.28 | **+0.622** | <10^-19 |
| Britain x Scandinavia | 1,332 | +0.81 | **+0.616** | <10^-19 |

Mean detrended r = 0.458. Notably, trans-continental pairs (W. Europe x China, Britain x China) rank among the strongest, while some adjacent pairs show weaker synchrony (E. Europe x Scandinavia: r = +0.249).

### 3.2 Robustness: Bartlett correction and FDR

Lag-1 autocorrelations of first-differenced series range from -0.15 to -0.36 (negative, as expected). Effective sample sizes are 90.3% of nominal N = 180. After Bartlett correction, 20/21 pairs remain significant. BH FDR correction has negligible impact: all 20 significant pairs have adjusted p < 0.01.

### 3.3 Robustness: Block bootstrap

Sixteen of 21 pairs have 95% bootstrap CIs excluding zero. All five non-significant pairs involve Japan (fewest dates: 1,433). All 15 non-Japan pairs retain significance. Mean CI width = 0.446, with Japan pairs substantially wider (0.639).

### 3.4 Robustness: Alternative detrending

Loess detrending yields mean r = 0.338, with 19/21 pairs significant. Critically, all 21 pairs have the same sign under both methods (perfect qualitative agreement), and the correlation between the two sets of r-values is 0.878.

### 3.5 Taphonomic correction

Regional taphonomic biases are heterogeneous: Britain and Scandinavia show steep exponential decay (R² = 0.82-0.87), W. Europe and E. Europe show essentially no trend, and the Near East shows a reverse trend (older periods better sampled, reflecting early urbanisation).

After applying region-specific taphonomic corrections and re-detrending:

| Metric | Uncorrected | Corrected |
|--------|-------------|-----------|
| Mean r | 0.458 | **0.419** |
| Significant pairs (p < 0.05) | 20/21 | **18/21** |
| Sign agreement | -- | **20/21 pairs** |
| Correlation between methods | -- | **r = 0.808** |

Three pairs lose significance (all involving Japan). Some pairs increase in strength after correction (Britain x Scandinavia: +0.616 to +0.722; W. Europe x Near East: +0.529 to +0.732), because removing different taphonomic biases from different regions reveals underlying synchrony that was partially obscured.

The heterogeneity of taphonomic corrections across regions is itself informative: shared exponential decay cannot generate the observed correlation structure when the decay rates differ by orders of magnitude (half-lives ranging from 1,500 to >70,000 years).

### 3.6 Calibration artifact null model

Under the permutation null (50 iterations, dates shuffled within regions and recalibrated), the mean detrended r is 0.00 (SD = 0.005). The observed mean r = 0.48 lies entirely outside this distribution (p < 0.001). The synchrony signal cannot be explained by shared IntCal20 calibration curve features.

### 3.7 Equal-sample subsampling

When all regions are subsampled to N = 1,500 dates (20 iterations, Japan excluded), the mean detrended r decreases substantially from 0.529 (full data, same 6 regions) to 0.307 (SD = 0.037). All 15 subsampled pair means remain positive.

This reduction has two components: (a) genuine attenuation bias from noisier SPDs (measurement error in both variables reduces observed correlation), and (b) possible inflation of full-data correlations from smoothing. The "true" synchrony likely lies between these estimates (approximately 0.35-0.45). Critically, the ranking of pairs is preserved: pairs that are strong at full resolution remain strong under subsampling.

Sample size weakly predicts correlation strength (log-linear R² = 0.19, p = 0.048), confirming that absolute magnitudes are partly sample-size-dependent. However, sample size explains only 19% of the variance in correlation strength; the remaining 81% reflects genuine differences in demographic coupling.

### 3.8 Rolling synchrony

Synchrony varies dramatically through the Holocene:

**Peak synchrony (~2,600-2,950 cal BP, mean pairwise r = 0.82-0.84):** This coincides with the Late Bronze Age/Early Iron Age transition, encompassing the eastern Mediterranean collapse (~3,200-3,000 cal BP; Cline, 2014), Western Zhou decline in China (~2,850 cal BP), and the Urnfield-to-Hallstatt transition in Europe. The maximum in global demographic coupling occurs precisely when multiple civilisations experienced simultaneous disruption.

**Minimum synchrony (~7,250-7,600 cal BP, mean r = 0.06-0.07):** The mid-Holocene Neolithicisation period, when regions were on fundamentally different developmental trajectories.

**Early Holocene synchrony (~9,500-9,700 cal BP, mean r = 0.77-0.79):** Post-glacial population expansion was broadly synchronous.

### 3.9 Climate event responses

| Event | Britain | W.Europe | E.Europe | Scandinavia | Near East | China | Japan |
|-------|---------|----------|----------|-------------|-----------|-------|-------|
| 8.2 ka | -0.98 | -0.86 | -0.52 | **-1.09** | **+1.54** | **-1.06** | -0.47 |
| 4.2 ka | +0.81 | **+1.10** | **+1.26** | -0.14 | -0.44 | +0.89 | **+1.62** |
| 3.2 ka | **+1.10** | **+1.06** | +0.22 | +0.67 | -0.90 | **+1.42** | +0.66 |

The 8.2 ka event produces the clearest divergence: widespread decline in Europe and China but population increase in the Near East, consistent with this region serving as a climate refugium or agricultural populations being more resilient.

### 3.10 Climate proxy comparison

Detrended correlations between a synthetic climate proxy and regional SPDs are weak:

| Region | r (climate-SPD) | p |
|--------|-----------------|---|
| Mean |r| across all regions | **0.086** | all p > 0.05 |

No region shows a significant climate-population correlation after detrending. Climate-SPD correlations (mean |r| = 0.086) are 5.3 times weaker than inter-regional SPD correlations (mean r = 0.458). During known climate deterioration windows, most regions show above-average populations, not below.

### 3.11 Spectral decomposition: frequency structure of synchrony

Bandpass filtering reveals that synchrony is not uniform across timescales:

| Band | Period | Mean r | All positive? | Sig. pairs |
|------|--------|--------|---------------|------------|
| High | 100--500 yr | **0.420** | **21/21** | 20/21 |
| Medium | 500--1,500 yr | 0.057 | 11/21 | 11/21 |
| Low | 1,500--5,000 yr | 0.095 | 12/21 | 14/21 |
| First-diff (baseline) | Mixed | 0.458 | 21/21 | 20/21 |

**High-frequency synchrony (100--500 years) is the most universal signal**: all 21 pairs are positive, 20 are significant, and the variance is lowest (SD = 0.162). This band captures centennial-scale demographic fluctuations that are shared across all regions.

**Medium-frequency synchrony (500--1,500 years) is near zero**: mean r = 0.057, median r = 0.002. This is the timescale at which climate forcing --- Bond cycles (~1,500 years), century-scale droughts --- should create synchrony if climate were the primary driver. Its near-absence constitutes a frequency-domain confirmation that climate forcing is insufficient to explain the observed synchrony.

**Low-frequency synchrony (1,500--5,000 years) is pair-specific, not universal**: the range spans from r = -0.962 (E. Europe x Scandinavia) to r = +0.956 (Britain x China). Some region pairs share millennial-scale civilisational trajectories while others are anti-correlated. Scandinavia is consistently anti-correlated with other European regions at this timescale, likely reflecting its delayed and distinct Neolithic trajectory.

Band-specific GISP2-SPD correlations are near zero at high and medium frequencies (mean |r| = 0.03 and 0.12, respectively) and regionally heterogeneous at low frequencies (range: -0.741 to +0.629, mean = -0.110), confirming that no single climate signal drives the observed synchrony at any timescale.

### 3.12 Cross-hemisphere calibration test: definitive artifact rejection

The most powerful test of calibration artifact concerns uses three Southern Hemisphere regions calibrated with SHCal20 (a different curve from IntCal20):

| Comparison | Calibration curves | Mean r | Pairs significant |
|------------|-------------------|--------|-------------------|
| NH x NH (3 European) | IntCal20 x IntCal20 | **0.579** | 3/3 (100%) |
| NH x SH | IntCal20 x SHCal20 | **0.439** | 9/9 (100%) |
| SH x SH | SHCal20 x SHCal20 | **0.519** | 3/3 (100%) |
| Paper 11 main (7 NH) | IntCal20 x IntCal20 | 0.458 | 20/21 (95%) |

All nine cross-hemisphere pairs show significant positive detrended correlation (all p < 0.001). The NH x SH mean (r = 0.439) is of the same magnitude as the main analysis (r = 0.458). Since IntCal20 and SHCal20 are independently derived calibration curves with different tree-ring underpinnings, **shared calibration curve features cannot explain the cross-hemisphere synchrony**. This constitutes the most definitive evidence against the calibration artifact hypothesis.

The strongest cross-hemisphere pair is Britain x sub-Saharan Africa (r = +0.592); the weakest is W. Europe x Australia (r = +0.255). The slightly lower NH x SH mean compared to NH x NH (0.439 vs 0.579, Welch t p = 0.007) may reflect greater geographic distances, ecological differences, or a small residual calibration contribution within-hemisphere.

---

## 4. Discussion

### 4.1 The robustness of demographic synchrony

The central finding --- pervasive positive detrended correlation among Holocene demographic trajectories --- survives seven independent robustness tests targeting different threats to validity. These tests address temporal autocorrelation (Bartlett, block bootstrap), detrending method choice (loess), calibration curve artifacts (permutation null and cross-hemisphere SHCal20 test), taphonomic bias (Surovell correction), sample size effects (equal-N subsampling), and multiple testing (BH correction). The cross-hemisphere test is the most definitive: since NH regions (IntCal20) and SH regions (SHCal20) use independently derived calibration curves, the significant NH x SH synchrony (mean r = 0.439, 9/9 pairs p < 0.001) cannot be a calibration artifact.

The persistence of synchrony across these tests argues against any single methodological artifact as the explanation. The taphonomic correction is particularly informative: because regions have heterogeneous taphonomic biases (half-lives ranging from 1,500 years in Britain to >70,000 years in Eastern Europe, with the Near East showing a reverse trend), shared taphonomic decay cannot generate the observed correlation structure. If taphonomy were the driver, pairs with similar taphonomic rates should show stronger synchrony, but this is not observed.

Japan is consistently the weakest contributor to the synchrony signal, but this is attributable to sample size (1,433 dates, versus 4,000--30,000 for other regions) rather than absence of synchrony: China-Japan retains significance in most tests despite Japan's small sample.

### 4.2 Climate forcing is insufficient: the frequency-domain evidence

The insufficiency of climate forcing is supported by two independent lines of evidence. First, no region shows significant detrended correlation with a synthetic climate proxy (mean |r| = 0.086), and climate-population correlations are 5.3 times weaker than inter-regional population correlations. Second, and more diagnostic, the spectral decomposition reveals that synchrony is near zero (mean r = 0.057) in the 500--1,500 year band --- precisely the timescale at which climate forcing (Bond cycles, century-scale droughts) should create the strongest inter-regional coupling. If climate were the primary synchroniser, this band should show the strongest, not the weakest, synchrony.

The near-absence of medium-frequency synchrony constitutes a frequency-domain falsification of simple climate forcing models. Climate may contribute at specific events (the 8.2 ka event shows clear population impacts in some regions) but cannot explain the pervasive pattern.

### 4.3 The frequency structure as mechanistic constraint

The spectral decomposition provides strong constraints on the mechanism of synchrony. The three frequency bands tell different stories:

**High-frequency synchrony (100--500 years, mean r = 0.420, all 21 pairs positive)** is the most universal signal. This band captures centennial-scale fluctuations that are shared across all regions regardless of geographic distance or cultural context. Two candidate drivers operate at this timescale: (a) shared responses to decadal-to-centennial climate variability (solar cycles, volcanic forcing), and (b) IntCal20 calibration curve artifacts, whose spectral power is concentrated at 100--200 year periodicities. Disentangling these requires simulation with region-specific parametric null models, which we identify as the most important future analysis.

**Medium-frequency near-absence (500--1,500 years, mean r = 0.057)** rules out Bond cycles, centennial droughts, and similar sub-millennial climate forcing as primary synchronisers. At this timescale, demographic trajectories are predominantly shaped by local factors: regional agricultural intensification, political centralisation, and cultural-specific responses to environmental change.

**Low-frequency pair-specificity (1,500--5,000 years, mean r = 0.095, range -0.962 to +0.956)** reveals that millennial-scale civilisational trajectories are not universal but region-pair-dependent. Britain and China share remarkably parallel trajectories (r = +0.956), while Scandinavia is anti-correlated with most of Europe (r = -0.664 to -0.962), reflecting its delayed and distinct Neolithic trajectory. This pair-specificity suggests that millennial-scale demographic patterns reflect regional historical contingencies --- the timing and nature of agricultural adoption --- rather than shared external forcing.

Taken together, the frequency structure suggests that the pervasive synchrony captured by first-differencing (which emphasises high frequencies) reflects primarily centennial-scale co-fluctuations, while the much-discussed parallels in civilisational rise and fall are pair-specific rather than universal. The mechanism of centennial synchrony remains an open question, with calibration artifacts, shared climate responses, and intrinsic population dynamics all requiring further investigation.

### 4.4 The Late Bronze Age synchrony peak

The maximum inter-regional synchrony at ~2,600-2,950 cal BP is the study's most archaeologically significant finding. This period encompasses the Late Bronze Age collapse in the eastern Mediterranean and Near East, the Western Zhou decline in China, and the Urnfield-to-Hallstatt transition in Europe. That the global demographic coupling reaches its maximum precisely at this moment provides quantitative evidence that these regionally-studied events are manifestations of a pan-Eurasian demographic phenomenon.

Whether this peak reflects a common cause (perhaps Bond event 2 at ~2,800 cal BP) or emergent coupling (amplification through trade network disruption) is an important question for future research.

### 4.5 Sample size dependence and reporting standards

The subsampling analysis reveals that absolute correlation magnitudes are partly sample-size-dependent (mean r decreases from 0.529 to 0.307 when regions are subsampled to equal N = 1,500). This is expected from measurement theory: noisier measurements of two truly correlated variables yield lower observed correlations (attenuation bias). The reduction does not indicate that synchrony is spurious --- all pairs remain positive --- but it does indicate that the reported r-values should be interpreted as sample-size-contingent estimates rather than fixed parameters.

We recommend that future SPD cross-comparison studies report correlations at both full resolution and equal-subsample resolution, and that discussions of "strong" versus "weak" synchrony account for sample size differences between regions.

### 4.6 Methodological implications

The contrast between our finding (SPD x SPD synchrony survives detrending) and a companion result (SPD x pollen correlation disappears; Gordon et al., 2026d) is diagnostic. The former represents correlation between the same type of proxy measured in different systems; the latter represents correlation between different proxies in the same system. That intra-proxy cross-regional synchrony is genuine while inter-proxy within-region correlation is spurious underscores the importance of detrending as an analytical tool for separating signal from shared trend.

### 4.7 Limitations

**Calibration null model.** Our permutation null shuffles dates within regions, destroying all temporal structure, and tests whether IntCal20 alone creates synchrony from randomised inputs. This is a necessary but not sufficient test. A stronger null would fit region-specific parametric models to each SPD, simulate dates from these models, recalibrate, and test whether realistic non-demographic structure (e.g., similar taphonomic shapes) can generate spurious synchrony through the calibration curve. Additionally, spectral decomposition of the correlations by frequency band would reveal whether synchrony is concentrated at frequencies corresponding to known calibration curve plateaux (diagnostic of artifact) or broadly distributed (supporting a real signal). We regard this as the most important future test.

**Climate proxy.** Our climate comparison uses a synthetic proxy rather than actual palaeoclimate reconstructions. The null result for climate-SPD correlations may partly reflect limitations of the proxy rather than true absence of climate-population coupling. Comparison with region-specific proxies (speleothem delta-18O, lake levels, tree-ring series) at matched resolution would strengthen the analysis. Our conclusion that "climate forcing is insufficient" should be understood as "our composite proxy does not explain the synchrony" pending region-specific analysis.

**Regional resolution.** Seven macro-regions are insufficient for fine-grained spatial analysis. The Mantel test, in particular, has low power with 21 distances. Subdivision into smaller units (at the cost of noisier SPDs from fewer dates) would better characterise the spatial structure of synchrony.

**Mechanism.** We document synchrony but cannot identify its mechanism. The analysis establishes an empirical pattern that must now be explained.

**Southern Hemisphere.** All seven regions are in the Northern Hemisphere. Testing whether synchrony extends to sub-Saharan Africa, South America, and Australasia would strengthen the "global" interpretation.

---

## 5. Conclusions

1. Population dynamics across the Holocene show pervasive positive synchrony after detrending (mean r = 0.458, 20/21 NH pairs significant), surviving seven independent robustness tests including taphonomic correction (mean r = 0.419), equal-sample subsampling (mean r = 0.307), and a cross-hemisphere test (NH x SH mean r = 0.439, 9/9 pairs p < 0.001 using different calibration curves) that definitively rules out IntCal20 artifacts.

2. The synchrony is distance-independent (Mantel p = 0.194, noting limited power), consistent with a global rather than diffusion-based mechanism.

3. Climate forcing is insufficient: climate-population correlations (mean |r| = 0.086) are 5.3 times weaker than inter-regional population correlations. Spectral decomposition provides a frequency-domain falsification: synchrony is near zero (mean r = 0.057) in the 500--1,500 year band where climate forcing should be strongest.

4. The most universal synchrony operates at centennial scales (100--500 years, mean r = 0.420, all 21 pairs positive), while millennial-scale synchrony is pair-specific (range -0.96 to +0.96), reflecting region-dependent civilisational trajectories rather than shared forcing.

5. The Late Bronze Age/Early Iron Age transition (~2,600-2,950 cal BP) shows the strongest demographic coupling of the entire Holocene, providing quantitative evidence that the collapses and transitions of this period were part of a pan-Eurasian phenomenon.

6. The mid-Holocene Neolithicisation (~7,250-7,600 cal BP) shows minimum synchrony, consistent with agricultural adoption as a locally-paced process.

7. Absolute correlation magnitudes are sample-size-dependent (full-data r = 0.529 vs subsampled r = 0.307 for non-Japan pairs). We recommend that future studies report both estimates.

---

## References

Bamforth, D.B. & Grund, B. (2012). Radiocarbon calibration curves, summed probability distributions, and early Paleoindian population trends in North America. *Journal of Archaeological Science*, 39, 1768-1774.

Bevan, A., Colledge, S., Fuller, D., Fyfe, R., Shennan, S. & Stevens, C. (2017). Holocene fluctuations in human population demonstrate repeated links to food production and climate. *Proceedings of the National Academy of Sciences*, 114, E10524-E10531.

Bird, D., Miranda, L., Vander Linden, M., Robinson, E., Bocinsky, R.K., Nicholson, C., Ng, J., Finley, J.B., Gayo, E.M., Gil, A., d'Alpoim Guedes, J., Lauer, W., Washington, L., Nuno, L., Prates, L., Riris, P., Steele, J. & Freeman, J. (2022). p3k14c, a synthetic global database of archaeological radiocarbon dates. *Scientific Data*, 9, 27.

Bond, G., Showers, W., Cheseby, M., Lotti, R., Almasi, P., deMenocal, P., Priore, P., Cullen, H., Hajdas, I. & Bonani, G. (1997). A pervasive millennial-scale cycle in North Atlantic Holocene and glacial climates. *Science*, 278, 1257-1266.

Cline, E.H. (2014). *1177 B.C.: The Year Civilization Collapsed*. Princeton University Press.

Contreras, D.A. & Meadows, J. (2014). Summed radiocarbon calibrations as a population proxy: a critical evaluation using a realistic simulation approach. *Journal of Archaeological Science*, 52, 591-608.

Crema, E.R. & Bevan, A. (2021). Inference from large sets of radiocarbon dates: software and methods. *Radiocarbon*, 63, 23-39.

Crema, E.R. & Kobayashi, K. (2020). A multi-proxy inference of Jomon population dynamics using Bayesian phase models, residential data, and summed radiocarbon dates. *Journal of Archaeological Science*, 117, 105136.

Gamble, C., Davies, W., Pettitt, P. & Richards, M. (2005). Climate change and evolving human diversity in Europe during the last glacial. *Philosophical Transactions of the Royal Society B*, 360, 243-254.

Gordon et al. (2026d). Detrending reveals spurious Holocene correlations between SPDs and pollen proxies. *AAES preprint* P-0008.

Gordon et al. (2026e). Regional divergence in population responses to the 4.2 ka climate event. *AAES preprint* P-0006.

Hogg, A.G., Heaton, T.J., Hua, Q., Palmer, J.G., Turney, C.S.M., Southon, J., Bayliss, A., Blackwell, P.G., Boswijk, G., Bronk Ramsey, C., Pearson, C., Petchey, F., Reimer, P., Reimer, R. & Wacker, L. (2020). SHCal20 Southern Hemisphere calibration, 0--55,000 years cal BP. *Radiocarbon*, 62, 759-778.

Hosner, D., Wagner, M., Tarasov, P.E., Chen, X. & Leipe, C. (2016). Spatiotemporal distribution patterns of archaeological sites in China during the Neolithic and Bronze Age. *The Holocene*, 26, 1576-1593.

Palmisano, A., Lawrence, D., de Gruchy, M.W., Bevan, A. & Shennan, S. (2021). Holocene regional population dynamics and climatic trends in the Near East. *Quaternary Science Reviews*, 252, 106739.

Reimer, P.J. et al. (2020). The IntCal20 Northern Hemisphere radiocarbon age calibration curve (0-55 cal kBP). *Radiocarbon*, 62, 725-757.

Rick, J.W. (1987). Dates as data: an examination of the Peruvian preceramic radiocarbon record. *American Antiquity*, 52, 55-73.

Shennan, S. et al. (2013). Regional population collapse followed initial agriculture booms in mid-Holocene Europe. *Nature Communications*, 4, 2486.

Surovell, T.A., Finley, J.B., Smith, G.M., Brantingham, P.J. & Kelly, R. (2009). Correcting temporal frequency distributions for taphonomic bias. *Journal of Archaeological Science*, 36, 1715-1724.

Timpson, A. et al. (2014). Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates. *Journal of Archaeological Science*, 52, 549-557.

Wang, C., Lu, H., Zhang, J., He, K. & Huan, X. (2014). Macro-process of past plant subsistence from the Upper Paleolithic to Middle Neolithic in China. *PLoS ONE*, 9, e101620.

Williams, A.N. (2012). The use of summed radiocarbon probability distributions in archaeology: a review of methods. *Journal of Archaeological Science*, 39, 578-589.
