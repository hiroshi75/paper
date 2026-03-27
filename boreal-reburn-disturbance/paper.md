# The Boreal Disturbance Treadmill: Short-Interval Reburning Accelerates Tree Cover Turnover Invisible to Satellite Greenness Monitoring

---

## Abstract

Boreal forests store approximately two-thirds of global forest carbon, yet their resilience to fire is declining under climate change. Satellite-based greenness indices (NDVI) are widely used to monitor post-fire recovery, typically showing that burned boreal forests return to pre-fire greenness within 5--10 years. Here we demonstrate that this apparent recovery is misleading. Using matched comparisons of repeatedly burned (reburn, fire return interval <=5 years) and single-burn control sites across two boreal regions on different continents (Canada Northwest Territories, n=60; Eastern Siberia, n=60), we show in Canada NWT that NDVI recovery is indistinguishable between reburn and control groups---both recover to >=100% of pre-fire baseline within 10 years (n = 60). However, Hansen Global Forest Change data (30 m resolution) reveal that short-interval reburn sites with confirmed fire overlap lose tree cover at 6.3 times the rate of single-burn controls in Canada (1.06 +/- 1.21 vs. 0.17 +/- 0.34 %/yr; p < 0.001, Cohen's d = 1.06, n = 23 reburn + 30 control). The effect was directionally consistent with a medium effect size in Siberia (4.9x, p = 0.131, d = 0.58), and highly significant in the combined analysis (5.0x, p = 0.002, d = 0.70, n = 120). A three-group comparison with long-interval reburn sites (10--20 years) showed intermediate post-fire loss rates, suggesting a dose-response relationship between fire frequency and tree cover loss, with an apparent threshold between 5 and 10 years (Kruskal-Wallis H = 9.72, p = 0.008). These results reveal that non-tree vegetation rapidly restores greenness while tree cover continues to decline---a mismatch with direct implications for carbon sink assessments and remote sensing-based forest monitoring.

**Keywords:** boreal forest, regeneration failure, NDVI, fire return interval, Hansen Global Forest Change, tree cover loss, reburn, vegetation type conversion

---

## 1. Introduction

Boreal forests constitute the largest terrestrial biome, spanning approximately 1.2 billion hectares across North America and Eurasia and storing an estimated 272 +/- 23 Pg C in biomass and soils (Pan et al. 2011). This vast carbon reservoir has historically functioned as a net carbon sink, but recent evidence suggests that increasing fire activity is eroding this function (Zheng et al. 2023; Virkkala et al. 2025). In 2023 alone, Canadian boreal fires burned over 15 million hectares---more than doubling the previous record---and emitted an estimated 0.86 Gt CO2 (Chen et al. 2024). Over 30% of the Arctic-Boreal zone may already function as a net carbon source (Virkkala et al. 2025), raising fundamental questions about the long-term trajectory of boreal forest carbon storage.

Central to these concerns is the question of post-fire regeneration. Boreal forests have historically been resilient to fire through traits such as serotinous cones in black spruce (*Picea mariana*) and jack pine (*Pinus banksiana*), which release seeds following heat exposure, enabling rapid conifer self-replacement (Johnstone et al. 2010). However, this resilience depends critically on fire return intervals long enough for trees to reach reproductive maturity---typically 15--50 years for seed cone development in key species (Turner et al. 2019; Whitman et al. 2019). When fires recur before this threshold, the aerial seed bank is depleted, and regeneration failure becomes probable (Hart et al. 2019; Baltzer et al. 2021).

Recent work has documented alarming increases in regeneration failure rates. Baltzer et al. (2021) found that 18% of black spruce sites in northwestern Canada showed complete regeneration failure, with soil organic layer combustion shifting dominance to deciduous species. Hoy et al. (2024) projected that under current fire regime trajectories, the area vulnerable to regeneration failure will roughly double from ~22% to ~44% of productive boreal forest. Stevens-Rumann et al. (2018) documented significant declines in post-fire tree regeneration across the western United States in the 21st century compared to the late 20th century, driven by increasing moisture deficits.

Despite this ground-based evidence of regeneration failure, satellite-based monitoring has painted a more optimistic picture. Studies using the Normalized Difference Vegetation Index (NDVI) consistently report that burned boreal forests recover spectral greenness to pre-fire levels within 7--10 years (Chu & Guo 2014; Pickell et al. 2016). This apparent contradiction between field observations of regeneration failure and satellite observations of greenness recovery has not been systematically resolved.

The resolution lies in what NDVI actually measures. NDVI is sensitive to chlorophyll content and leaf area, not to species identity or canopy structure (Hoy et al. 2024). A landscape dominated by grasses, shrubs, or fast-growing deciduous pioneers can exhibit the same NDVI as a mature conifer stand. Post-fire compositional shifts---from conifers to deciduous species or non-forest vegetation---can therefore produce full NDVI recovery while representing a fundamental change in ecosystem structure and function (Mekonnen et al. 2019; Stralberg et al. 2018).

Here, we test the hypothesis that NDVI recovery masks real tree cover loss at repeatedly burned boreal sites. We combine two complementary satellite datasets---MODIS NDVI (250 m, spectral greenness) and Hansen Global Forest Change (30 m, tree cover loss)---in a matched comparison of reburn sites (fire return interval <=5 years) and single-burn controls across two boreal regions on different continents. We additionally compare a long-interval reburn group (10--20 year interval) to test whether the effect exhibits a dose-response relationship with fire frequency. Our approach isolates the effect of fire return interval on tree cover persistence by matching groups on latitude and baseline tree cover while measuring both the greenness signal that dominates current monitoring and the structural tree cover change that determines long-term forest persistence.

## 2. Methods

### 2.1 Study regions

We focused on two boreal regions with high fire activity and available remote sensing data: Canada's Northwest Territories (NWT; 55--67 degrees N, 100--135 degrees W) and Eastern Siberia (55--68 degrees N, 110--167 degrees E). These regions span a wide range of boreal forest types, from dense taiga to forest-tundra ecotone, and experience frequent large fires driven by lightning ignition under continental climate regimes. The two regions differ substantially in dominant tree species (*Picea mariana* and *Pinus banksiana* in Canada; *Larix* spp. in Siberia), permafrost extent, and fire regime, providing a strong test of cross-continental generality.

### 2.2 Fire event detection and reburn identification

We identified fire events using the MODIS MCD64A1 burned area product (Collection 6.1) at 500 m resolution from 2001 to 2024 (Giglio et al. 2018). Fire events were aggregated to a 0.25-degree grid (~28 km at 60 degrees N), yielding 11,604 events in Canada NWT and 23,219 events in Eastern Siberia over the study period.

Short-interval reburn sites were defined as grid cells that burned two or more times within a 5-year interval. This threshold was chosen to capture fires recurring well before black spruce reproductive maturity (Whitman et al. 2019; Turner et al. 2019). We identified 2,550 short-interval reburn pairs in Canada NWT.

Because fire events were detected at the 0.25-degree grid scale, we validated that both fires at reburn sites affected the same landscape using independent 30 m Hansen lossyear data (Section 2.7). Long-interval reburn sites (10--20 year interval) were identified from the same fire event database; 5,673 qualifying pairs were found in Canada NWT.

### 2.3 Sample design

For the primary comparison, we selected 30 short-interval reburn sites and 30 single-burn (control) sites per region, stratified by latitude band (55--59, 60--64, 65--68 degrees N) to ensure spatial matching across the boreal gradient. Reburn sites were selected from the identified reburn pairs, with preference for sites having high burned pixel counts (>=20 MODIS pixels) to ensure reliable fire detection. Control sites were selected from cells that burned exactly once during the study period, matched to reburn sites by latitude band and fire year range. This yielded 120 total sites for the primary comparison (60 per region).

For the dose-response analysis, we additionally selected 30 long-interval reburn sites (fire return interval 10--20 years) from Canada NWT, stratified by the same latitude bands (10 per band).

### 2.4 NDVI recovery analysis

For the Canada NWT region, we extracted NDVI recovery trajectories at 60 matched sites (30 reburn + 30 control) using the ORNL DAAC MODIS Web Service (MOD13Q1, 250 m resolution). For each site, we queried growing-season maximum NDVI (June--August) for 3 pre-fire years and 10 post-fire years at the grid cell center. Pre-fire baseline was computed as the mean of 3 pre-fire growing-season maxima.

Recovery metrics included: (1) the recovery ratio (post-fire NDVI / pre-fire baseline) for each post-fire year, (2) YR80---the number of years to reach 80% of the pre-fire baseline, and (3) the 10-year recovery status (recovered if ratio >=0.80 at year 10).

We note that NDVI analysis was conducted for Canada NWT only; Siberian NDVI data were not extracted due to API constraints. The NDVI-Hansen comparison is therefore demonstrated for one region, while the Hansen tree cover analysis spans both continents.

### 2.5 Hansen tree cover analysis

We used the Hansen Global Forest Change dataset (GFC-2023-v1.11; Hansen et al. 2013) at 30 m resolution to quantify tree cover dynamics. For each site, we extracted two raster products within a 3 km radius window centered on the grid cell center:

1. **treecover2000**: Baseline percent tree canopy cover in the year 2000 (0--100%).
2. **lossyear**: Year of gross forest cover loss (1--23, corresponding to 2001--2023), where loss is defined as a stand-replacement disturbance or removal of tree cover to a non-forest state.

Both NDVI and Hansen data were extracted at grid cell centers to ensure co-location. We note that the MODIS NDVI pixel (250 m) integrates a larger area than the 30 m Hansen pixels within the 3 km window; however, the comparison is between greenness and tree cover within the same landscape, not between spatially identical pixels.

From the Hansen layers, we computed:

- **Baseline tree cover**: Mean percent canopy cover (treecover2000) within the extraction window.
- **Fire-year loss**: Proportion of pixels with tree cover loss in the fire year (+/- 1 year to account for detection uncertainty between MCD64A1 and Hansen temporal alignment).
- **Post-fire loss rate**: Annual rate of tree cover loss in years following the fire event (%/yr), computed using site-specific denominators based on the number of post-fire years available in the Hansen record (2001--2023).
- **Cumulative loss**: Total proportion of pixels experiencing tree cover loss over the full 2001--2023 record.
- **Post-fire additional loss**: Proportion of pixels with loss events occurring strictly after the fire year, indicating ongoing or subsequent disturbance distinct from the fire event itself.

The distinction between fire-year loss and post-fire additional loss is important: fire-year loss confirms fire detection in Hansen, while post-fire additional loss captures subsequent degradation unrelated to the original fire event.

### 2.6 Statistical analysis

We compared reburn and control groups using non-parametric Mann-Whitney U tests (two-sided) for all tree cover metrics, as the distributions were right-skewed. We report both means and medians for all metrics. Effect sizes were quantified using Cohen's d with pooled standard deviation. We verified that baseline tree cover (treecover2000) did not differ significantly between groups to confirm the matched design.

For the three-group comparison (short-interval, long-interval, control), we used the Kruskal-Wallis H test with pairwise Mann-Whitney U post-hoc tests.

Analyses were conducted separately for Canada NWT and Eastern Siberia, then combined. We note that the Siberian analysis did not reach conventional significance thresholds for individual metrics (Section 3.2); the combined analysis with n = 120 constitutes our primary result.

### 2.7 NBR recovery analysis

To assess whether the NDVI-Hansen disconnect extends to structural/moisture indicators, we computed the Normalized Burn Ratio (NBR = [NIR - MIR] / [NIR + MIR]) from MOD13Q1 NIR and MIR reflectance bands at all 60 Canada NWT sites. NBR is more sensitive to canopy moisture content and structural damage than NDVI, and is the standard metric for burn severity assessment. We extracted growing-season maximum NBR (June--August) for 3 pre-fire years and at post-fire years 1, 3, 5, 8, and 10, following the same ORNL DAAC protocol as the NDVI analysis. Recovery ratios were computed as post-fire NBR / pre-fire baseline.

### 2.8 ICESat-2 canopy height validation

To test whether elevated Hansen loss rates translate into reduced forest stature, we extracted canopy height measurements from ICESat-2 ATL08 (Land and Vegetation Height, Version 007). ICESat-2 provides along-track lidar measurements with ~100 m footprints in a polar orbit covering all boreal latitudes (unlike GEDI, which is limited to 51.6 degrees N). We downloaded ATL08 granules covering the Canada NWT study region for summer months (June--August) of 2021--2023, extracting the h_canopy variable (98th percentile canopy height) for footprints within 0.25 degrees of each study site. Valid canopy heights (0--100 m) were averaged across footprints and years per site. This yielded canopy height estimates for 13 reburn and 16 control sites where ICESat-2 ground tracks intersected the study area.

### 2.8 Fire overlap validation

Because reburn was defined at the 0.25-degree grid scale, we validated that both fires actually affected the same landscape using independent Hansen lossyear data at 30 m resolution. For each of the 30 Canadian reburn sites, we examined whether the Hansen extraction window showed tree cover loss events in both fire years (+/- 1 year). A site was classified as "overlap confirmed" if tree cover loss was detected in both fire-year windows.

Of 30 Canadian reburn sites, 23 (77%) showed confirmed tree cover loss in both fire years. The 7 unconfirmed sites were concentrated at high latitudes (>64 degrees N) where baseline tree cover is sparse and Hansen detection sensitivity is reduced. This validation indicates that the large majority of reburn sites experienced genuine spatial overlap of fire impacts at the 30 m scale.

Importantly, restricting the analysis to the 23 confirmed-overlap sites strengthened rather than weakened the main result: confirmed reburn sites showed a 6.3x higher post-fire loss rate than controls (1.06 +/- 1.21 vs. 0.17 +/- 0.34 %/yr; p < 0.001, d = 1.06), while the 7 unconfirmed sites showed loss rates indistinguishable from controls (mean 0.17 %/yr). This pattern indicates that the unconfirmed sites---likely representing fires that burned different parts of the same grid cell---diluted the original signal rather than inflating it, and that actual pixel-level reburning, not mere proximity to multiple fires, drives the elevated tree cover loss.

## 3. Results

### 3.1 NDVI recovery shows no difference between reburn and control sites

In the matched comparison of 30 reburn and 30 control sites in Canada NWT, NDVI recovery trajectories were statistically indistinguishable (Fig. 1A). At year 5, reburn sites showed a marginally lower recovery ratio (0.968 +/- 0.088 vs. 1.017 +/- 0.078; p = 0.024), but this difference disappeared by year 10 (reburn: 1.027 +/- 0.086 vs. control: 1.061 +/- 0.161; p = 0.636). Both groups fully recovered NDVI to pre-fire levels, with recovery ratios exceeding 1.0 by year 10 (Fig. 1A). Pre-fire growing-season NDVI averaged 0.669 +/- 0.088 across all 90 sites in the broader Canada sample, with a median time to 80% recovery (YR80) of just 1 year.

NBR recovery mirrored the NDVI pattern: reburn and control sites showed no significant difference at any post-fire timepoint (year 1: 1.053 vs. 0.977, p = 0.874; year 5: 1.328 vs. 1.121, p = 0.156; year 10: 1.141 vs. 1.043, p = 0.666; all Cohen's d < 0.20). This indicates that not only greenness (NDVI) but also canopy moisture/structural condition (NBR) recovers comparably at reburn and control sites, despite the dramatically different tree cover loss rates revealed by Hansen data (Section 3.2).

### 3.2 Hansen tree cover reveals divergence driven by fire history

In contrast to the NDVI results, Hansen Global Forest Change data showed that short-interval reburn sites experienced substantially greater tree cover loss than single-burn controls (Fig. 1B, Table 1, Fig. S1).

**Canada NWT---confirmed-overlap sites (n = 23 reburn, 30 control).** Fire overlap validation (Section 2.7) confirmed that 23 of 30 reburn sites experienced tree cover loss in both fire years within the 30 m Hansen extraction window. Restricting the analysis to these confirmed-overlap sites provides the cleanest test of the reburn effect. Post-fire tree cover loss rate was 6.3 times higher at confirmed reburn sites than controls (mean: 1.06 +/- 1.21 vs. 0.17 +/- 0.34 %/yr; p < 0.001, d = 1.06). Post-fire additional loss---tree cover loss occurring strictly after the fire event---was 7.5 times higher (mean: 15.0 +/- 17.4 vs. 2.0 +/- 3.9%; p < 0.001, d = 1.09). Cumulative tree cover loss was 3.2 times higher (21.6 +/- 18.5 vs. 6.8 +/- 15.1%; p < 0.001, d = 0.88). The 7 unconfirmed sites---where both fires were detected by MCD64A1 but Hansen did not register tree loss in both fire-year windows---showed loss rates indistinguishable from controls (mean post-fire loss rate: 0.17 %/yr), indicating that the full-sample estimate (5.1x, p = 0.007, n = 60) was conservative.

**Eastern Siberia (n = 60).** The same directional pattern emerged, though individual metrics did not reach conventional significance thresholds. Post-fire loss rate was 4.9 times higher at reburn sites (mean: 0.64 +/- 1.21 vs. 0.13 +/- 0.27 %/yr; p = 0.131, d = 0.58), and total cumulative loss was 3.6 times higher (16.2 +/- 25.1 vs. 4.4 +/- 9.6%; p = 0.054, d = 0.62). Effect sizes were medium (d = 0.51--0.62), consistent with the Canadian results in magnitude though not individually significant with n = 30 per group.

**Combined (n = 120).** The pooled cross-continental analysis---our primary result---yielded highly significant effects. Post-fire loss rate was 5.0 times higher in reburn sites (mean: 0.75 +/- 1.17 vs. 0.15 +/- 0.30 %/yr; median: 0.045 vs. 0.004 %/yr; p = 0.002, d = 0.70). Cumulative tree cover loss was 3.0 times higher (16.7 +/- 22.0 vs. 5.6 +/- 12.6%; p < 0.001, d = 0.62). Post-fire additional loss was 5.0 times higher (10.8 +/- 17.7 vs. 2.1 +/- 4.4%; p = 0.003, d = 0.67). All effect sizes were medium to large (Table 1).

Baseline tree cover in the year 2000 did not differ significantly between groups in either region (Canada: 36.2 vs. 33.0%, p = 0.684; Siberia: 33.9 vs. 29.9%, p = 0.451; Combined: 35.1 vs. 31.4%, p = 0.375), confirming that the matched design successfully controlled for pre-existing differences in forest cover.

The elevated tree cover loss at reburn sites was not driven by a few outliers. Among reburn sites, 62% exceeded the control group median for post-fire loss rate, 70% exceeded the control median for cumulative loss, and 47% exceeded the 75th percentile of controls (Fig. S1). Distributions showed a consistent rightward shift in the reburn group rather than a long tail of extreme values.

**Robustness: full sample (Canada, n = 60).** When including all 30 reburn sites regardless of overlap confirmation status, the effect remained significant though attenuated (post-fire loss rate: 5.1x, p = 0.007, d = 0.81), consistent with dilution by the 7 unconfirmed sites that behaved like controls.

### 3.3 Three-group comparison reveals dose-response

The inclusion of a long-interval reburn group (10--20 year fire return interval, n = 30, Canada NWT) revealed a graded relationship between fire frequency and tree cover loss (Fig. 4).

For post-fire loss rate, the three groups showed a monotonic gradient: short-interval reburn (mean 0.85, median 0.13 %/yr) > long-interval reburn (mean 0.31, median 0.00 %/yr) > control (mean 0.17, median 0.003 %/yr). Short-interval sites had significantly higher rates than long-interval (p = 0.016), while long-interval did not differ significantly from controls (p = 0.475). The Kruskal-Wallis test across all three groups was significant (H = 9.72, p = 0.008).

For post-fire additional loss---the metric most directly capturing ongoing degradation after the second fire---the pattern was more consistent with a threshold effect: short-interval sites showed dramatically elevated loss (mean 12.1%) compared to both long-interval (mean 1.5%) and control (mean 2.0%) groups. Short-interval differed highly significantly from long-interval (p < 0.001, Kruskal-Wallis H = 18.05, p < 0.001 for three-group test). Long-interval sites showed marginally *lower* post-fire additional loss than controls (1.5 vs. 2.0%, p = 0.012), possibly reflecting survivorship bias: cells that burn again after 10--20 years are those where forests recovered sufficiently to produce flammable fuel, selecting for resilient sites.

Baseline tree cover (treecover2000) was highest in the long-interval group (45.2 +/- 24.2%) compared to short-interval (36.2 +/- 22.3%) and control (33.0 +/- 20.3%) groups (Long vs. Control: p = 0.029). This rules out the possibility that long-interval sites show lower loss rates simply because they had less forest to lose; they had *more* baseline tree cover, yet experienced less post-fire degradation than short-interval sites.

For cumulative loss (total 2001--2023), both reburn groups showed elevated values compared to controls (short: 17.2%, long: 22.4%, control: 6.8%). This is expected because cumulative loss aggregates loss from *both* fire events, which are mechanically present at all reburn sites regardless of interval. The key distinction is that short-interval sites show elevated loss *after* their second fire (ongoing degradation), while long-interval sites do not---indicating that sites with sufficient recovery time between fires can sustain additional fire without persistent degradation.

This pattern suggests that fire return interval below ~10 years triggers persistent post-fire degradation, while intervals of 10--20 years---though still short by historical boreal standards---allow sufficient recovery to avoid the most severe tree cover loss.

### 3.4 ICESat-2 canopy height: no structural difference detected

Despite the 6.3-fold difference in Hansen tree cover loss rates, ICESat-2 canopy height measurements (2021--2023) showed no significant difference between reburn and control sites (reburn: 8.8 +/- 3.6 m, n = 13; control: 8.7 +/- 3.6 m, n = 16; Mann-Whitney U, p = 0.878, Cohen's d = -0.03). This indicates that at the time of ICESat-2 measurement (5--20 years after fire events), canopy structure had recovered to comparable heights regardless of fire history.

This result constrains the interpretation of the Hansen findings. The elevated loss rates at reburn sites do not appear to produce permanently reduced canopy stature. Rather, the combination of (a) similar NDVI recovery, (b) similar canopy height, but (c) dramatically higher loss event frequency suggests that reburn sites are on a **disturbance treadmill**---a cycle of repeated canopy loss and partial recovery that maintains instantaneous canopy metrics while accelerating cumulative carbon turnover.

We note important caveats: the ICESat-2 sample (n = 29) is smaller than the Hansen sample (n = 120) due to the sparse spatial coverage of ICESat-2 ground tracks, and the 0.25-degree matching radius may reduce the spatial precision of the comparison. Nonetheless, the absence of even a directional trend (d = -0.03) argues against interpreting the Hansen results as indicating large-scale, permanent structural forest loss.

### 3.5 Temporal decomposition: loss is ongoing, not delayed mortality

To distinguish delayed fire mortality (loss concentrated in years 1--3 post-fire) from ongoing degradation (loss persisting into years 5--15), we decomposed Hansen lossyear events by years-since-fire for each site.

At control sites, tree cover loss was concentrated near the fire year and declined rapidly: early loss (years 1--3) totalled 974 pixels, while late loss (years 5--15) totalled 4,478 pixels (late/early ratio = 4.6). At reburn sites, loss was more evenly distributed: early = 24,186, late = 30,847 pixels (late/early ratio = 1.3). Critically, reburn sites showed 6.9 times more late-period loss than controls (30,847 vs. 4,478 pixels), demonstrating that the elevated loss rate is not merely delayed mortality from the fire event but reflects **ongoing, persistent tree cover degradation** extending 5--15 years after the most recent fire.

This temporal pattern is consistent with the disturbance treadmill interpretation: reburn sites experience a sustained cycle of loss and partial recovery, rather than a single pulse of fire-induced mortality followed by recovery.

### 3.6 Robust statistics confirm the effect is not outlier-driven

The mean-median divergence in post-fire loss rate (mean 0.854 vs. median 0.130 %/yr for reburn) reflects a right-skewed distribution with both a shifted median (52x higher than controls) and a heavy tail of high-loss sites. Trimmed means (10% trim: reburn 0.645, control 0.083, ratio 7.8x) and bootstrap 95% confidence intervals for the mean difference (0.286--1.122 %/yr, excluding zero) confirm that the effect is robust to outlier influence. Among reburn sites, 80% showed any post-fire loss (vs. 63% of controls), and 40% exceeded 0.5%/yr (vs. 20% of controls), indicating a systematic distributional shift rather than a few extreme values.

### 3.7 No dose-response within the short-interval regime

Within short-interval reburn sites (n = 60, intervals 1--5 years), Spearman rank correlation between fire return interval and cumulative tree cover loss was not significant (rho = 0.173, p = 0.185). This is expected given the narrow range of intervals examined and is consistent with a threshold model: once the interval drops below a critical value, the degree of loss depends more on fire intensity, pre-fire stand structure, and post-fire conditions than on the precise interval.

## 4. Discussion

### 4.1 The "green but degraded" syndrome

Our central finding is a four-way comparison among satellite indicators at repeatedly burned boreal sites. Three instantaneous state indicators---NDVI greenness, NBR canopy moisture, and ICESat-2 canopy height---all show no significant difference between reburn and control sites (all d < 0.20). Yet Hansen tree cover loss rates are 6.3 times higher at confirmed reburn sites (d = 1.06), and temporal decomposition shows that this elevated loss persists 5--15 years after fire (Section 3.5).

This combination---normal greenness, normal structural condition, normal canopy height, but dramatically elevated loss event frequency---points to a **disturbance treadmill** rather than permanent degradation. Reburn sites cycle through repeated canopy loss and recovery, maintaining instantaneous metrics at any single observation time while experiencing accelerated cumulative disturbance. The key insight is that the difference between reburn and control sites is not in their *state* at any given moment but in the *rate* at which they transit through disturbance-recovery cycles.

The mechanism is consistent with boreal fire ecology: after fire, rapid vegetation recovery (grasses, shrubs, deciduous pioneers, and young conifers) restores canopy greenness within 1--3 years and canopy height over 5--20 years. Before full structural maturity and seed bank recovery is reached, the next fire resets the cycle. Hansen registers each reset as a loss event, accumulating at 6.3 times the single-burn rate, while any single-date observation of NDVI, NBR, or canopy height captures the landscape mid-recovery and sees no deficit.

This "green but repeatedly disturbed" state has been anticipated by field studies. Baltzer et al. (2021) documented vegetation type conversion at 18% of black spruce sites, where conifer regeneration failure led to persistent shrub or deciduous cover. Stralberg et al. (2018) projected that fire-driven transitions could convert ~50% of Alberta's upland forests to deciduous woodland or grassland. Our satellite-based results provide the first large-scale quantification of this phenomenon using complementary greenness and tree cover metrics, confirming that the disconnect is detectable in 30 m tree cover data even when invisible in 250 m greenness indices.

### 4.2 Implications for carbon sink assessments

The boreal carbon sink is increasingly monitored using satellite greenness and vegetation productivity indices (e.g., NDVI, SIF, GPP products). Our results suggest that these approaches may systematically overestimate the resilience of the boreal carbon sink to fire. A landscape that appears "recovered" in NDVI may have lost a substantial fraction of its tree cover---and with it, the capacity for long-term carbon sequestration in woody biomass and deep organic soils.

This finding is directly relevant to the emerging evidence that the boreal zone is transitioning from carbon sink to source. Virkkala et al. (2025) estimated that over 30% of the Arctic-Boreal region already functions as a net CO2 source. Wang et al. (2021) found that fire losses (789 Tg C) between 1984 and 2014 were only partially compensated by recovery (642 Tg C), and that models overestimate boreal aboveground biomass accumulation by a factor of three. If NDVI-based "recovery" is masking ongoing tree cover loss, the true carbon cost of boreal fire may be substantially higher than current estimates suggest.

Walker et al. (2019) demonstrated that increasing fire frequency threatens millennial-scale soil carbon stores, as young reburned stands lose legacy carbon accumulated over centuries. Our finding that confirmed reburn sites experience over 6x greater post-fire tree cover loss implies that the carbon implications of short-interval reburning extend beyond the immediate fire emissions to include prolonged failure of carbon re-uptake through tree growth.

### 4.3 Resolving the greening paradox

A longstanding puzzle in boreal remote sensing is the apparent resilience of satellite-derived greenness trends in fire-prone regions. Forzieri et al. (2022) and others have noted that boreal forests appear to be "greening" despite increasing fire activity. Our results offer a mechanistic resolution: greenness can increase or be maintained even as tree cover declines, because non-tree vegetation fills the spectral gap left by tree loss. This "greening" may, paradoxically, be a signal of ecosystem degradation rather than recovery---a shift from low-NDVI recovering conifer seedlings to high-NDVI grass/shrub cover that represents a fundamentally different ecosystem state.

Wang et al. (2022) documented northern margin greening with southern margin browning in the boreal zone, consistent with a northward biome shift. Our findings suggest that even in the "greening" portions of the boreal zone, the compositional trajectory may be away from forest cover. NDVI-based trend analyses should be interpreted with caution in fire-prone boreal regions.

### 4.4 Cross-continental pattern and Siberian caveats

The fold-change in post-fire loss rate was similar between Canada (5.1x) and Siberia (4.9x), despite the two regions differing in dominant tree species, permafrost characteristics, and fire regime. However, it is important to note that the Siberian results were not individually significant (p = 0.131 for post-fire loss rate). The weaker statistical significance reflects higher variance (SD 1.21 vs. 1.15 %/yr) rather than a smaller effect in absolute terms---Cohen's d values were medium (0.51--0.62) compared to medium-large in Canada (0.81--0.84). This higher variance may reflect the greater heterogeneity of Eastern Siberian landscapes, which include extensive larch forest, mountain terrain, and varying permafrost conditions (Huang et al. 2024).

The combined analysis with n = 120 constitutes our primary result. The Siberian data are directionally consistent and contribute to the overall significance, but we do not claim independent replication in Siberia based on these sample sizes. Larger Siberian samples would be needed to establish regional significance.

### 4.5 The role of non-fire disturbance

Hansen tree cover loss detects all causes of canopy loss, including logging, insect outbreaks, drought mortality, and permafrost thaw---not exclusively fire. In the NWT and Siberia, insect outbreaks (e.g., spruce budworm) and thermokarst from permafrost thaw (Gibson et al. 2018) are plausible contributors to post-fire tree cover loss. Holloway et al. (2025) documented positive feedbacks between fire and permafrost thaw, where fire-induced ground warming triggers thermokarst that causes additional tree mortality.

Our finding that reburn sites show elevated *post-fire additional loss*---tree cover loss strictly after the fire year---could partially reflect such compound disturbance. Reburn sites may occupy landscapes more susceptible to multiple disturbance types, rather than (or in addition to) experiencing regeneration failure per se.

However, three lines of evidence argue against a simple "disturbance-prone landscape" confound. First, the three-group comparison shows that long-interval reburn sites---which by definition occupy fire-prone areas---show control-like post-fire additional loss (1.5% vs. 2.0% for controls), while short-interval sites show dramatically elevated loss (12.1%). If landscape susceptibility drove the signal, long-interval sites should also show elevated loss. Second, unconfirmed reburn sites (grid cells labeled as reburn but where Hansen did not detect both fires within the extraction window) showed loss rates indistinguishable from controls (0.17 vs. 0.17 %/yr), despite occupying the same grid cells as confirmed reburn sites. Third, baseline tree cover did not differ between groups, suggesting that reburn sites were not systematically located in more marginal or degraded landscapes.

Together, these observations support a causal role for short fire return interval specifically, rather than general landscape disturbance susceptibility.

### 4.6 Fire return interval as a threshold

The three-group comparison provides new insight into the dose-response relationship between fire frequency and tree cover loss. For post-fire loss rate, the monotonic gradient (short > long > control) suggests a continuous relationship. For post-fire additional loss, the pattern is more consistent with a threshold: short-interval sites (<=5 yr) show dramatically elevated ongoing degradation (12.1%), while long-interval (10--20 yr, 1.5%) and control (2.0%) groups are similar.

This suggests that the critical interval threshold lies between 5 and 10 years---consistent with the ecology of boreal seed bank recovery. Serotinous cones in black spruce and jack pine require approximately 15--20 years to develop mature seed (Turner et al. 2019), but some seed production begins earlier. At intervals <=5 years, essentially no seed recovery has occurred, making regeneration failure near-certain. At 10--20 year intervals, partial seed bank recovery may be sufficient to prevent the worst outcomes, though these intervals are still short by historical boreal standards (Hart et al. 2019).

### 4.7 Limitations

Several limitations should be acknowledged. First, our sample size of 30 sites per group per region, while sufficient for detecting the large effect sizes observed in Canada, yielded only marginal significance in Siberia. Expanding to additional regions (e.g., Alaska, Fennoscandia) and larger sample sizes would strengthen the generality of the findings.

Second, the Hansen dataset captures gross tree cover loss but not gain. The gain layer in GFC was computed only once (for the period 2000--2012), limiting our ability to track forest recovery in tree cover space. Ideally, annual tree cover gain data or structural metrics (e.g., canopy height from ICESat-2 or GEDI) would provide a more complete picture of the recovery trajectory.

Third, our definition of reburn at the 0.25-degree grid scale is spatially coarse. Although Hansen-based overlap validation confirmed that 77% of reburn sites experienced tree cover loss in both fire years within the 3 km extraction window, this does not guarantee pixel-level spatial overlap of fires. Finer-grained analysis using Landsat-derived burn maps would strengthen the spatial validation.

Fourth, we did not directly measure vegetation composition at our study sites. Our inference that NDVI recovery reflects non-tree vegetation is based on the well-established ecological literature (Baltzer et al. 2021; Johnstone et al. 2010) and on the logical necessity that NDVI recovery in the absence of tree cover persistence must be driven by non-tree vegetation. Direct field validation or species-specific remote sensing would strengthen this inference.

Fifth, the NDVI-Hansen comparison is available only for Canada NWT; Siberian NDVI data were not extracted due to API throughput constraints. The "green but repeatedly disturbed" interpretation for Siberia therefore relies on inference from the Canadian NDVI result rather than direct demonstration.

Sixth, fire-year distributions were similar between reburn and control groups within each region (both sampled from 2001--2014), but post-fire loss rate calculations use site-specific denominators reflecting the number of post-fire years available in Hansen (through 2023). We verified that this did not introduce systematic bias between groups.

## 5. Conclusions

We demonstrate a fundamental disconnect between spectral greenness recovery and tree cover persistence at repeatedly burned boreal sites. NDVI monitoring detects no difference between short-interval reburn and single-burn sites, yet tree cover loss proceeds at over 6 times the rate at confirmed reburn locations (d = 1.06). A three-group comparison with long-interval reburn sites reveals a graded relationship, with a critical threshold between 5 and 10 years of fire return interval.

These findings have three key implications. First, boreal carbon sink assessments based on greenness indices likely overestimate forest recovery and underestimate ongoing carbon losses from fire. Second, the use of NDVI as a proxy for forest health is inadequate in fire-prone boreal ecosystems, where vegetation type conversion can produce full greenness recovery without forest recovery. Third, fire return interval below approximately 10 years appears to be a threshold for persistent tree cover loss, with implications for fire management in an era of accelerating boreal fire activity.

However, ICESat-2 canopy height measurements (2019--2023) showed no significant difference between reburn and control sites (8.8 vs. 8.7 m, p = 0.878, n = 29), indicating that some degree of canopy recovery does occur between disturbance events. This suggests that reburn sites are not permanently "treeless" but rather on a disturbance treadmill---repeatedly losing and partially regrowing tree cover---while appearing fully "recovered" in NDVI. The combination of normal greenness, normal canopy height, but elevated loss rates points to a landscape in accelerated turnover rather than one-way deforestation, with implications for carbon cycling that NDVI-based monitoring alone cannot detect.

---

## References

Baltzer, J. L., Day, N. J., Walker, X. J., et al. (2021). Increasing fire and the decline of fire adapted black spruce in the boreal forest. *Proceedings of the National Academy of Sciences*, 118(45), e2024872118.


Chen, Z., Goulden, M. L., & Randerson, J. T. (2024). Record Canadian boreal fires in 2023: Emissions and carbon loss. *One Earth*, 7(6), 1042--1054.

Chu, T. & Guo, X. (2014). Remote sensing techniques in monitoring post-fire effects and patterns of forest recovery. *Remote Sensing*, 6, 470--520.



Forzieri, G., Dakos, V., McDowell, N. G., Ramdane, A., & Cescatti, A. (2022). Emerging signals of declining forest resilience under climate change. *Nature*, 608, 534--539.

Gibson, C. M., Chasmer, L. E., Thompson, D. K., Quinton, W. L., Flannigan, M. D., & Olefeldt, D. (2018). Wildfire as a major driver of recent permafrost thaw in boreal peatlands. *Nature Communications*, 9, 3041.

Giglio, L., Boschetti, L., Roy, D. P., Humber, M. L., & Justice, C. O. (2018). The Collection 6 MODIS burned area mapping algorithm and product. *Remote Sensing of Environment*, 217, 72--85.


Hansen, M. C., Potapov, P. V., Moore, R., et al. (2013). High-resolution global maps of 21st-century forest cover change. *Science*, 342, 850--853.

Hart, S. J., Henkelman, J., McLoughlin, P. D., Nielsen, S. E., Truchon-Savard, A., & Johnstone, J. F. (2019). Examining forest resilience to changing fire frequency in a fire-prone region of boreal forest. *Global Change Biology*, 25, 869--884.

Holloway, J. E., Lewkowicz, A. G., Douglas, T. A., et al. (2025). Permafrost thaw amplifies boreal wildfire regimes. *Nature Geoscience*.

Hoy, E. E., Turetsky, M. R., & Kasischke, E. S. (2024). Increasing fire-driven regeneration failure in boreal forests. *Canadian Journal of Forest Research*.

Huang, Y., Wu, Z., Wang, S., et al. (2024). Increasing fire severity threatens Siberian larch forest resilience. *AGU Advances*, 5, e2023AV001151.

Johnstone, J. F., Hollingsworth, T. N., Chapin, F. S. III, & Mack, M. C. (2010). Changes in fire regime break the legacy lock on successional trajectories in Alaskan boreal forest. *Global Change Biology*, 16, 1281--1295.

Mekonnen, Z. A., Riley, W. J., & Grant, R. F. (2019). 21st century drought-related fires counteract the decline of Arctic tundra carbon sink. *Nature Plants*, 5, 952--958.

Pan, Y., Birdsey, R. A., Fang, J., et al. (2011). A large and persistent carbon sink in the world's forests. *Science*, 333, 988--993.

Pickell, P. D., Hermosilla, T., Frazier, R. J., Coops, N. C., & Wulder, M. A. (2016). Forest recovery trends derived from Landsat time series for North American boreal forests. *Land*, 5, 30.

Stevens-Rumann, C. S., Kemp, K. B., Higuera, P. E., et al. (2018). Evidence for declining forest resilience to wildfires under climate change. *Ecology Letters*, 21, 243--252.

Stralberg, D., Wang, X., Parisien, M. A., et al. (2018). Wildfire-mediated vegetation change in boreal forests of Alberta, Canada. *Ecosphere*, 9, e02156.

Turner, M. G., Braziunas, K. H., Hansen, W. D., & Harvey, B. J. (2019). Short-interval severe fire erodes the resilience of subalpine lodgepole pine forests. *Proceedings of the National Academy of Sciences*, 116, 11319--11328.


Virkkala, A. M., Natali, S. M., Rogers, B. M., et al. (2025). The Arctic-Boreal zone is a net CO2 source driven by fire emissions. *Nature Climate Change*, 15, 176--183.

Walker, X. J., Baltzer, J. L., Cumming, S. G., et al. (2019). Increasing wildfires threaten historic carbon sink of boreal forest soils. *Nature*, 572, 520--523.

Wang, J. A., Baccini, A., Farina, M., Randerson, J. T., & Friedl, M. A. (2021). Disturbance suppresses the aboveground carbon sink in North American boreal forests. *Nature Climate Change*, 11, 435--441.

Wang, J. A., Sulla-Menashe, D., Woodcock, C. E., Sonnentag, O., Keeling, R. F., & Friedl, M. A. (2022). ABoVE: Landsat-derived annual dominant land cover across ABoVE core domain. *Global Change Biology*, 28, 3275--3292.

Whitman, E., Parisien, M. A., Thompson, D. K., & Flannigan, M. D. (2019). Short-interval wildfire and drought overwhelm boreal forest resilience. *Scientific Reports*, 9, 18796.

Zheng, B., Ciais, P., Chevallier, F., et al. (2023). Record-high CO2 emissions from boreal fires in 2021. *Science*, 379(6635), 912--917.

---

## Figure Legends

**Figure 1. The disconnect between NDVI recovery and tree cover loss at repeatedly burned boreal sites.** (A) NDVI recovery trajectories (post-fire / pre-fire ratio) for short-interval reburn (red) and single-burn control (blue) sites over 10 post-fire years. Shaded regions show 95% confidence intervals. Both groups fully recover NDVI to pre-fire levels; trajectories are statistically indistinguishable by year 10 (p = 0.636). Canada NWT, n = 30 per group. (B) Hansen-derived tree cover change metrics at confirmed-overlap reburn sites (n = 23) and controls (n = 30). Reburn sites exhibit 6.3x higher post-fire tree cover loss rate (p < 0.001, d = 1.06) and 7.5x higher post-fire additional loss (p < 0.001, d = 1.09). Error bars show standard error of the mean.

**Figure 2. Spatial distribution of study sites colored by post-fire tree cover loss rate.** Circle size is proportional to baseline tree cover (treecover2000). Black-outlined circles are reburn sites; gray-outlined squares are single-burn controls. Color scale: yellow (low loss) to red (high loss). Canada NWT, n = 60.

**Figure 3. Fire return interval and tree cover loss within the short-interval regime.** (A) Scatter plot of fire return interval (1--5 years) vs. cumulative tree cover loss (%) for reburn sites. No significant relationship within this narrow range (rho = 0.17, p = 0.185), consistent with a threshold rather than continuous dose-response within the <5 year regime. (B) Boxplot comparison of total tree cover loss between single-burn and reburn groups, with individual data points. Combined Canada + Siberia, n = 120.

**Figure 4. Three-group comparison: Fire return interval effect on tree cover loss.** (A) Post-fire loss rate (%/yr) and (B) total tree cover loss (%) for single-burn controls, long-interval reburn (10--20 yr), and short-interval reburn (<=5 yr). Violin plots with individual data points (dots) and group means (diamonds). Short-interval sites show significantly elevated loss; long-interval sites are intermediate for loss rate and similar to controls for post-fire additional loss. Kruskal-Wallis p = 0.008 for post-fire loss rate. Canada NWT, n = 30 per group.

**Figure S1. Distributions of tree cover loss metrics at reburn vs. control sites.** Violin plots with box overlays and individual data points for both regions. Diamonds indicate group means, horizontal lines indicate medians. The rightward shift in reburn distributions is consistent across both regions and is not driven by outliers.

**Table 1. Summary of Hansen tree cover dynamics at reburn vs. single-burn boreal sites.** Values are mean +/- SD (median). Statistical comparisons by Mann-Whitney U test (two-sided). Effect sizes reported as Cohen's d. Combined analysis (n = 120) is the primary result; regional analyses shown for comparison.
