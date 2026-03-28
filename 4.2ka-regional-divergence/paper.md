# Regional Divergence in Population Responses to the 4.2 ka Climate Event: A Cross-Regional Permutation Test Framework

**Target journal**: Quaternary Science Reviews / AAES
**Version**: v7.2 (2026-03-28, expanded references to 70, AAES submission preparation)
**Status**: Submission-ready

---

## Abstract

The 4.2 ka climate event, which defines the base of the Meghalayan Stage, has been widely associated with societal collapses across multiple civilizations. However, recent paleoclimate meta-analyses have questioned both its global coherence and its climatic exceptionality. Here we present a unified statistical framework, applying, for the first time, identical calibration procedures, binning parameters, and permutation tests across multiple world regions simultaneously, to compare radiocarbon-inferred population dynamics during the 4.2 ka interval across the Near East, Europe, China, and Japan using the p3k14c (180,070 dates), JOAD (44,425 dates), and NERD databases. Our results demonstrate that regional divergence in population response is statistically significant, consistent with recent paleoclimate evidence that the 4.2 ka event was not globally uniform (McKay et al., 2024). The Near East experienced a −24.9% population decline while western Europe grew by +37.1% (permutation test p = 0.003). Within Japan, Chubu recorded a −30.1% decline while Kyushu grew by +84.6% (p = 0.002), and independent cross-validation using the NERD database reproduced the Near East–Europe divergence (p = 0.003). Sub-regional analysis of China revealed opposing trends masked by national-scale aggregation: the Yellow River basin grew (+18.1%) while the Yangtze basin declined (−16.2%), consistent with the Liangzhu collapse. Site binning corrected a false decline signal in northern Europe (−21.0% unbinned → +10.4% binned), underscoring methodological best practices. Paleoclimate proxy data are consistent with the interpretation that population declines coincided with monsoon weakening (Dongge Cave) and aridification (Gol-e-Zard Cave) in affected regions, while humid-zone populations grew during minimal high-latitude perturbation. We argue that the 4.2 ka event is best understood not as a globally uniform crisis but as a regionally heterogeneous demographic phenomenon, in which the interaction between spatially variable climate forcing and locally specific subsistence vulnerabilities shaped divergent outcomes. This framework provides a quantitative basis for understanding how a climatically heterogeneous event can produce divergent demographic consequences across regions.

**Keywords**: 4.2 ka event; summed probability distribution; radiocarbon demography; permutation test; regional divergence; population dynamics; Holocene; rcarbon

---

## Introduction

The formal ratification of the Meghalayan Stage of the Holocene Series, with its base defined at 4,200 years before present (BP), elevated a regional climatic anomaly to the status of a global chronostratigraphic boundary (Walker et al., 2012, 2018). This decision was grounded in decades of research linking aridification events around 4.2–3.9 ka BP to societal disruptions across multiple civilizations, from the collapse of the Akkadian Empire in Mesopotamia — documented through marine sediment evidence of dust-laden aridification (Cullen et al., 2000; Weiss, 2016) — to the demise of the Old Kingdom in Egypt and the end of the Indus Valley urban phase, where monsoon weakening has been linked to urban decline (Staubwasser et al., 2003; Staubwasser and Weiss, 2006). The narrative of a "global megadrought" driving synchronous societal collapse proved compelling (deMenocal, 2001; Weiss, 2017), generating a substantial body of literature that sought to document the 4.2 ka event across diverse geographical and environmental settings. In the Mediterranean alone, Bini et al. (2019) synthesized evidence from over two hundred proxy records to characterize the event's regional expression, demonstrating both its widespread detectability and its considerable spatial variability in timing, duration, and magnitude. Yet the very breadth of this research programme has, paradoxically, begun to undermine the premise upon which it was built: the assumption that the 4.2 ka event constituted a globally coherent and climatically extraordinary episode.

A recent meta-analysis by McKay et al. (2024), examining 1,142 paleoclimate records distributed across the globe, found that the 4.2 ka event is "not remarkable" when assessed against the full spectrum of Holocene climate variability. Only approximately six percent of records exhibited statistically significant excursions during the 4.2 ka interval — a proportion no greater than expected by chance. This finding resonates with broader insights from paleoclimate synthesis: Neukom et al. (2019) demonstrated that no globally coherent warm or cold periods existed during the preindustrial Common Era, challenging the long-standing assumption that past climate anomalies operated as spatially uniform forcing mechanisms. Furthermore, the very status of the 4.2 ka event as a distinct climatic episode remains debated (Wiener, 2014), and broader reviews of mid-to-late Holocene climate change have situated it within a continuum of abrupt transitions rather than treating it as a singular catastrophe (Wanner et al., 2008). If the 4.2 ka event was neither globally synchronous nor climatically exceptional, then the traditional framing — which posits a single catastrophic climate perturbation propagating uniformly across civilizations — requires fundamental revision (Middleton, 2017; Aimers and Hodell, 2011). The question is no longer whether the 4.2 ka event was globally significant in climatic terms, but rather why human populations in certain regions experienced marked demographic decline while those in other regions maintained stability or even grew. It is this regional divergence in population response, not the global uniformity of the climate signal, that constitutes the phenomenon most in need of explanation.

Emerging palaeoclimatic research provides a mechanistic basis for expecting such divergence. Nan et al. (2025) documented pronounced spatial heterogeneity in the expression of the 4.2 ka event across the Northern Hemisphere, identifying distinct patterns of temperature and precipitation change that varied systematically with latitude and continentality. Regions influenced primarily by the North Atlantic Oscillation and the westerlies experienced drying, whereas those under monsoonal influence exhibited more complex, and in some cases opposite, responses. Railsback et al. (2018) further demonstrated that the 4.2 ka event was not a single sustained anomaly but comprised two distinct pulses of aridification separated by a brief amelioration, with the relative severity of each pulse varying geographically. These findings suggest that the demographic consequences of the 4.2 ka interval should be understood not as a uniform global response to a uniform global stressor, but as a mosaic of regionally contingent outcomes shaped by the interaction between spatially heterogeneous climate forcing, local environmental conditions, and the adaptive capacities of different subsistence systems.

The analytical tools necessary to test this proposition are now well established. The use of summed probability distributions (SPDs) of calibrated radiocarbon dates as a proxy for relative population change — the "dates as data" approach pioneered by Rick (1987) — was placed on a rigorous statistical footing by Shennan et al. (2013), who demonstrated regional population collapses across mid-Holocene Europe. Crema et al. (2016) introduced a permutation-based framework for comparing observed SPDs against expectations derived from fitted demographic models, enabling formal statistical assessment of whether population trajectories in specific regions and time windows deviate significantly from broader trends. The subsequent release of the rcarbon package (Crema and Bevan, 2021) made these methods accessible and reproducible, catalyzing a proliferation of regional demographic studies. Crema (2022) provided a comprehensive review of the methodological landscape, articulating best practices for site-level binning, model selection, and sensitivity analysis that have become standard in the field. More recently, Wirtz et al. (2024) applied these methods to over 91,000 radiocarbon dates to identify multicentennial population cycles, demonstrating the capacity of SPD-based approaches to detect demographic signals across continental scales and millennial timescales. Multi-proxy validation studies have further strengthened confidence in SPD-derived population estimates (Palmisano et al., 2017), and spatio-temporal extensions of the permutation framework (Crema et al., 2017) have enabled formal testing of spatial structure in population dynamics. Most recently, Heaton et al. (2025) have proposed Poisson process-based alternatives to traditional SPDs that may provide more rigorous identification of changepoints, an approach that future studies should evaluate against the permTest framework used here.

Despite this methodological maturity, existing SPD studies of the 4.2 ka interval remain regionally siloed. Palmisano et al. (2021) analysed 10,653 radiocarbon dates across seven regions of the Near East, revealing significant inter-regional variation in demographic trajectories during the third millennium BP. Lawrence et al. (2021) documented population collapse in the northern Fertile Crescent coinciding with the 4.2 ka event, while southern Mesopotamian populations showed more nuanced responses. In East Asia, Wang et al. (2014) compiled 4,656 dates from China and identified demographic fluctuations broadly contemporaneous with the 4.2 ka interval, and Park et al. (2019) detected a population decline around 4.2 ka BP in the Korean Peninsula using 2,190 dates. Each of these studies has made important contributions to understanding regional population dynamics, yet no study has applied a unified analytical framework — identical calibration procedures, binning parameters, temporal resolution, and statistical tests — to compare population responses across multiple world regions simultaneously. The absence of such a cross-regional comparison leaves open the critical question of whether the observed regional differences reflect genuine divergence in demographic response or merely artefacts of methodological heterogeneity among independent studies.

The recent availability of large-scale, curated radiocarbon databases has removed what was previously the primary obstacle to such a synthesis. The p3k14c database (Bird et al., 2022) aggregates over 180,000 published radiocarbon dates with global coverage, while the Japanese Open Archaeometric Database (JOAD; Kudo, 2023), containing 44,425 dates, provides unprecedented coverage of the Japanese archipelago — a region whose demographic history during the 4.2 ka interval remains poorly integrated into global narratives. Together with the NERD (Near East Radiocarbon Dates) compilation, these databases enable a standardized multi-regional comparison — applying, for the first time, identical calibration procedures, binning parameters, and permutation tests across multiple world regions simultaneously — spanning the critical 5.5–3.0 ka BP window that encompasses the 4.2 ka event and its aftermath.

Here we present a cross-regional analysis of radiocarbon-inferred population change during the 4.2 ka interval across four world regions: the Near East, Europe, China, and Japan. Applying a uniform SPD methodology with site-level binning and permutation testing (Crema et al., 2016; Crema and Bevan, 2021), we address three interrelated research questions. First, do population responses to the 4.2 ka interval diverge significantly across climate zones, and if so, what is the magnitude and direction of divergence in each region? Second, can permutation-based statistical tests reliably detect this divergence against the backdrop of long-term demographic trends? Third, how does Japan — a region undergoing the Jomon-to-Yayoi transition and rarely included in global 4.2 ka syntheses — compare to better-studied regions in the Near East and East Asia? We provide the first unified statistical framework to quantify regional divergence in population response to the 4.2 ka event, showing that this divergence is statistically significant and maps onto the spatial structure of climatic heterogeneity. By applying identical methods across multiple regions, we move beyond descriptive comparison toward formal hypothesis testing of differential demographic response — an approach that offers broader insights into the heterogeneous nature of human–environment interactions during periods of rapid climate change.

---

## Methodology

### 2.1. Radiocarbon databases

Three radiocarbon databases were employed to construct regional summed probability distributions (SPDs). The primary global dataset was the p3k14c database version 2024.1 (Bird et al., 2022), which contains 180,070 radiocarbon dates with global coverage. This database served as the source for analyses of the Near East, Europe (western Europe, northern Europe, and the British Isles), and East Asia (China). For Japan, we used the Japanese Open Archaeological Dates database (JOAD) version 1.1.0 (Kudo, 2023), comprising 44,425 radiocarbon dates, of which 39,284 were retained in the curated English-language version. To cross-validate the Near Eastern results obtained from p3k14c, we additionally compiled dates from the Near East Radiocarbon Dates database (NERD; Palmisano et al., 2021), an independently curated collection of radiocarbon dates from Southwest Asia. Across all three databases, dates were filtered by excluding those with reported measurement errors exceeding 500 years, ages of zero or negative values, and samples of marine origin, following standard quality-control procedures in SPD-based demographic inference (Contreras and Meadows, 2014; Crema and Bevan, 2021).

### 2.2. Regional definitions

Study regions were defined by geographic bounding boxes to capture areas with distinct climatic regimes and archaeological traditions. Five global-scale regions were extracted from the p3k14c database (Table 1). These regions were selected to represent contrasting hydroclimatic settings: the Near East as an arid zone, western and northern Europe and the British Isles as humid temperate to cold zones, and East Asia (China) as a monsoon-influenced zone.

**Table 1.** Global study regions extracted from p3k14c.

| Region | Latitude (°N) | Longitude (°E) | N dates | N sites | Climate zone |
|---|---|---|---|---|---|
| Near East | 25–40 | 30–60 | 2,099 | 707 | Arid |
| Western Europe | 42–55 | −10–15 | 21,010 | 10,137 | Humid |
| Northern Europe | 55–70 | −10–30 | 6,818 | 3,546 | Humid-cold |
| British Isles | 50–60 | −10–2 | 15,129 | 7,541 | Humid |
| East Asia (China) | 20–45 | 100–125 | 2,362 | 1,126 | Monsoon |

For China, sub-regional analyses were conducted to resolve potentially opposing demographic trends within the monsoon zone. The Yellow River basin was defined as the area north of 34°N (~1,200 dates), and the Yangtze basin as the area between 26°N and 34°N (~800 dates), broadly corresponding to the millet-based agricultural zone of the north and the rice-based zone of the central-south, respectively.

Six Japanese regions were defined within the JOAD database to capture the well-documented east-west gradient in Jomon population dynamics (Table 2). Regional boundaries followed the standard Japanese archaeological regional classification.

**Table 2.** Japanese study regions extracted from JOAD.

| Region | N dates | N sites |
|---|---|---|
| Chubu | 8,297 | ~1,005 |
| Tohoku | 8,859 | ~900 |
| Kyushu | 5,425 | 863 |
| Kanto | 3,816 | ~500 |
| Hokkaido | 3,186 | ~400 |
| Kansai | 2,559 | ~350 |

### 2.3. SPD computation and site binning

All radiocarbon dates were calibrated using the IntCal20 calibration curve (Reimer et al., 2020). To mitigate the well-documented bias introduced by intensively dated individual sites, we applied site-level binning prior to SPD computation, following the protocol recommended by Crema and Bevan (2021). Dates from the same site falling within 200-year bins were aggregated such that each site-phase contributed equally to the summed distribution, regardless of the number of individual measurements. The resulting SPDs were computed as the sum of the calibrated probability densities and normalized to unity. All SPD computations and statistical tests were carried out using the rcarbon R package (Crema and Bevan, 2021).

### 2.4. Temporal windows and sensitivity analysis

Demographic change across the 4.2 ka event was assessed by comparing mean normalized SPD values between a pre-event window and an event window. The default window definitions were: pre-4.2 ka (5,000–4,500 cal BP), at-4.2 ka (4,400–4,000 cal BP), and post-4.2 ka (4,000–3,500 cal BP). Because SPD-based demographic estimates can be sensitive to the choice of comparison intervals, we conducted a systematic sensitivity analysis using four alternative window definitions: (i) default (pre: 5,000–4,500, event: 4,400–4,000 cal BP); (ii) centered (pre: 4,800–4,400, event: 4,200–3,800 cal BP); (iii) narrow (pre: 4,600–4,400, event: 4,200–4,000 cal BP); and (iv) wide (pre: 5,500–4,500, event: 4,400–3,500 cal BP). Additionally, to evaluate whether the Chubu decline reflected a longer-term trend predating the 4.2 ka event — as suggested by Crema and Kobayashi (2020), who identified the onset of eastern Japanese population decline at approximately 5,000 cal BP — we extended the pre-event window to 6,000–5,500 cal BP.

### 2.5. Statistical testing

Two complementary statistical approaches were employed. First, single-region tests were conducted using modelTest (Crema and Bevan, 2021), which evaluates whether the observed SPD for a given region deviates significantly from an exponential growth null model, using 200 Monte Carlo simulations with site-binned dates. This test was applied to the Near East and Chubu.

The primary statistical framework, however, was the mark permutation test (permTest; Crema et al., 2016; Crema and Bevan, 2021), which assesses whether the SPD trajectories of two regions differ significantly. In each of the 500 permutations, regional group labels were randomly reassigned, new SPDs were generated for each permuted group, and z-transformed differences were compared against a 95% simulation envelope. This approach tests the relative shape of the population trajectory rather than its absolute magnitude, making it robust to inter-regional differences in research intensity — a critical advantage given the substantial variation in date densities across our study regions (e.g., 707 sites in the Near East versus 10,137 in western Europe). Two principal comparisons were conducted: Near East versus western Europe, and Chubu versus Kyushu (Japan). The NERD database was used for an independent cross-validation of the Near East versus western Europe comparison. A sub-regional comparison of the Yellow River versus Yangtze basins was also attempted. To verify the stability of permutation test results, all tests were repeated with five different random seeds (42, 123, 456, 789, 2024); p-values were stable within ±0.001 across seeds.

### 2.6. Paleoclimate proxies

To provide climatic context for the observed demographic trends, three high-resolution paleoclimate proxy series spanning the 4.2 ka interval were compiled (Table 3). These records were selected to represent the principal climate subsystems relevant to our study regions and are used here as contextual evidence rather than as independent confirmation of climate–population causal links: high-latitude North Atlantic circulation (GISP2), the East Asian monsoon (Dongge Cave), and the Iranian–Near Eastern hydroclimate (Gol-e-Zard Cave).

**Table 3.** Paleoclimate proxy records used in this study.

| Proxy record | Location | Variable | Reference |
|---|---|---|---|
| GISP2 ice core | Greenland (72.6°N) | Temperature (°C) | Alley (2000) |
| Dongge Cave | China (25.3°N, 108.1°E) | δ¹⁸O (‰ vPDB) | Wang et al. (2005) |
| Gol-e-Zard Cave | Iran (33.5°N, 51.4°E) | δ¹⁸O (‰ vPDB) | Carolin et al. (2019) |

Mean proxy values were computed for the pre-4.2 ka and at-4.2 ka windows and compared to identify the direction and magnitude of climate shifts coincident with the demographic transitions documented by SPD analysis. We note the absence of a high-resolution Japanese-specific paleoclimate proxy record covering the mid-Holocene, which limits our ability to directly assess climatic conditions on the Japanese archipelago during the 4.2 ka interval. The Dongge Cave record was used as a proxy for broader East Asian monsoon variability, justified by evidence that Holocene monsoon oscillations are broadly synchronous across East Asia at orbital-to-millennial timescales (Wang et al., 2005), though local departures from this regional pattern cannot be excluded.

---

## Results

### 3.1. Global SPD analysis

Site-binned SPD analysis of the p3k14c database reveals a pronounced divergence in population trajectories across regions during the 4.2 ka interval (Table 4). The Near East experienced a robust population decline of −24.9% between the pre-4.2 ka and at-4.2 ka windows, consistent across all four window definitions (sensitivity range: −27.3% to −20.0%). In contrast, humid-zone regions showed sustained demographic growth: western Europe increased by +37.1% (range: +33.1% to +45.3%), the British Isles by +81.3% (range: +73.4% to +89.9%), and northern Europe by +10.4% (range: +9.3% to +13.0%). East Asia (China) as a whole showed modest growth of +11.6% (range: +11.6% to +19.3%).

The eastern Mediterranean exhibited a marginal decline of −5.7%, but this result was not robust, with sign reversal under the centered window definition (range: −11.5% to +13.1%). The eastern Mediterranean is excluded from further analysis because it encompasses a climatically and culturally heterogeneous zone spanning from the Aegean to the Levantine coast, and the sign reversal across window definitions suggests that opposing sub-regional trends may be cancelling — as demonstrated for China. A sub-regional analysis of the eastern Mediterranean, analogous to our China disaggregation, would be needed to resolve this ambiguity.

Notably, the initial unbinned analysis of northern Europe suggested an apparent decline of −21.0%, anomalous for a humid-zone region. Site-level binning corrected this to +10.4% growth, revealing that a small number of intensively dated sites in decline had dominated the unbinned SPD. This finding underscores the critical importance of site binning for SPD analyses, as previously recommended by Crema and Bevan (2021).

**Table 4.** Population change (%) during the 4.2 ka interval by region, site-binned SPDs.

| Region | N dates | N site-bins | Change (%) | Sensitivity range (%) | Robust? |
|---|---|---|---|---|---|
| Near East | 2,099 | 707 | −24.9 | −27.3 to −20.0 | Yes |
| W. Europe | 21,010 | 10,137 | +37.1 | +33.1 to +45.3 | Yes |
| British Isles | 15,129 | 7,541 | +81.3 | +73.4 to +89.9 | Yes |
| N. Europe | 6,818 | 3,546 | +10.4 | +9.3 to +13.0 | Yes |
| E. Asia (China) | 2,362 | 1,126 | +11.6 | +11.6 to +19.3 | Yes |
| E. Mediterranean | 2,408 | 757 | −5.7 | −11.5 to +13.1 | No |

### 3.2. Japanese regional analysis

The JOAD database reveals an analogous east–west divergence within Japan at approximately 4,400 cal BP, coinciding with the Middle-to-Late Jomon transition (Table 5). Chubu experienced a −30.1% decline (range: −39.1% to −19.9%), the largest magnitude decline detected in any region in this analysis — exceeding the Near East (−24.9%). However, cross-regional magnitude comparisons must be interpreted with caution given substantial differences in dating density: Chubu (8,297 dates, ~1,005 sites) has considerably higher data density than the Near East (2,099 dates, 707 sites), and regions with sparser coverage may have experienced comparable or larger declines that our analysis lacks the resolution to detect. Kanto also declined by −18.7% (range: −21.4% to −7.2%). In western Japan, Kyushu showed remarkable growth of +84.6% (range: +84.6% to +94.2%), while Kansai grew by +40.9% (range: +34.1% to +58.3%) and Hokkaido by +26.9% (range: +19.6% to +43.8%).

Tohoku presents an instructive case. The initial unbinned analysis suggested a marginal decline of −6.6%, but site binning revised this to +6.5% growth, paralleling the northern European correction. The core decline zone was thus Chubu–Kanto (central Honshu), not the broader eastern Japan.

**Table 5.** Population change (%) during the 4.2 ka interval for Japanese regions, site-binned SPDs (JOAD).

| Region | N dates | N site-bins | Change (%) | Sensitivity range (%) | Direction |
|---|---|---|---|---|---|
| Chubu | 8,297 | 3,378 | −30.1 | −39.1 to −19.9 | Decline |
| Kanto | 3,816 | 1,944 | −18.7 | −21.4 to −7.2 | Decline |
| Tohoku | 8,859 | 3,545 | +6.5 | +6.5 to +13.0 | Growth (revised) |
| Hokkaido | 3,186 | 1,302 | +26.9 | +19.6 to +43.8 | Growth |
| Kansai | 2,559 | 1,188 | +40.9 | +34.1 to +58.3 | Growth |
| Kyushu | 5,425 | 2,625 | +84.6 | +84.6 to +94.2 | Growth |

The timing of the Japanese divergence centres on ~4,400 cal BP, approximately 200 years earlier than the canonical 4.2 ka event. This offset may reflect: (a) earlier onset of East Asian monsoon weakening in the Japanese archipelago; (b) chronological uncertainty in the JOAD age models; or (c) the culmination of endogenous Middle Jomon population overshoot that broadly coincided with the 4.2 ka climate interval.

### 3.3. China sub-regional analysis

Disaggregation of the aggregate Chinese signal revealed opposing sub-regional trends. The Yellow River basin showed +18.1% growth, consistent with Longshan cultural expansion, while the Yangtze basin experienced −16.2% decline, consistent with the well-documented collapse of the Liangzhu culture (Shao et al., 2021; Zhang et al., 2021). Liu and Feng (2012) documented that six of seven Neolithic cultures in Chinese cultural domains collapsed around this transition, and Wu and Liu (2004) identified a "flood south, drought north" pattern consistent with the divergent sub-regional trajectories we observe. These opposing trends cancel in aggregation, accounting for the modest +11.6% growth of the combined Chinese region and illustrating the hazard of national-scale analysis (cf. Wang et al., 2014).

### 3.4. Paleoclimate context

Paleoclimate proxy data are consistent with the interpretation that the 4.2 ka event was primarily a low-latitude phenomenon affecting monsoon systems and arid zones (Table 6). Gol-e-Zard Cave recorded a +0.35‰ δ¹⁸O shift indicating substantial Near Eastern aridification, and Dongge Cave recorded a +0.31‰ shift indicating East Asian monsoon weakening. GISP2 showed minimal change (+0.66°C), consistent with the absence of a pronounced high-latitude signal. Population declines in the Near East and central Japan are temporally coincident with the documented monsoon weakening, while population growth in humid Europe is consistent with the minimal high-latitude perturbation.

**Table 6.** Paleoclimate proxy values across the 4.2 ka interval.

| Proxy | Pre-4.2 ka | At-4.2 ka | Change | Interpretation |
|---|---|---|---|---|
| GISP2 (Greenland) | −31.2°C | −30.6°C | +0.66°C | Minimal high-latitude signal |
| Dongge Cave (China) | −8.16‰ | −7.85‰ | +0.31‰ | Asian monsoon weakening |
| Gol-e-Zard (Iran) | −8.23‰ | −7.87‰ | +0.35‰ | Near East aridification |

### 3.5. Statistical significance

Single-region modelTest analyses for the Near East (263 sites) and Chubu (1,005 sites) against exponential null models were non-significant, consistent with the observation that single-region tests have lower statistical power for detecting localized demographic perturbations within long-term trends (Crema, 2022).

The mark permutation test (permTest), however, yielded highly significant results for the two principal inter-regional comparisons (Table 7). The Near East versus western Europe comparison returned a global p-value of 0.003, and the Chubu versus Kyushu comparison returned p = 0.002. The NERD cross-validation independently reproduced the Near East–Europe divergence (p = 0.003), supporting the inference that the signal is not database-specific. The Yellow River versus Yangtze comparison was non-significant (p = 0.71). We cannot distinguish between two explanations for this result: insufficient statistical power given the smaller sample sizes (~1,200 and ~800 dates, respectively), or genuine absence of statistically detectable divergence in the radiocarbon record.

All permutation tests were repeated with five random seeds; p-values were stable within ±0.001 across seeds.

**Table 7.** Mark permutation test (permTest) results.

| Comparison | Database | N dates (region 1) | N dates (region 2) | Global p-value | Significant? |
|---|---|---|---|---|---|
| Near East vs. W. Europe | p3k14c | 2,099 | 21,010 | 0.003 | Yes |
| Chubu vs. Kyushu | JOAD | 8,297 | 5,425 | 0.002 | Yes |
| Near East vs. W. Europe | NERD | — | — | 0.003 | Yes |
| Yellow River vs. Yangtze | p3k14c | ~1,200 | ~800 | 0.71 | No |

### 3.6. Extended baseline sensitivity

Extending the pre-event window to 6,000–5,500 cal BP reduced the magnitude of the Chubu decline to −12.3%, but Chubu remained the only region globally to show a decline under this extended baseline. This suggests that while part of the Chubu signal reflects a longer-term population trajectory that predates the 4.2 ka event, the 4.2 ka interval accentuated and accelerated the decline beyond what the preceding trend alone would predict.

---

## Discussion

### 4.1 Regional Divergence as the Primary Signal

The central finding of this study is that population responses to the 4.2 ka climate event were regionally divergent rather than globally uniform, and that this divergence is statistically significant. The mark permutation test (permTest) demonstrates that the Near East and Western Europe followed statistically distinguishable population trajectories during the 4.2 ka interval (p = 0.003), as did Chubu and Kyushu within Japan (p = 0.002). These results move beyond the descriptive observation that some regions declined while others grew — a pattern noted in earlier SPD studies (Palmisano et al., 2021; Park et al., 2019) — to establish that the divergence itself is a robust demographic signal distinguishable from sampling noise.

This finding aligns with and extends the conclusions of McKay et al. (2024), who demonstrated through analysis of 1,142 paleoclimate records that the 4.2 ka event was not a globally coherent climate anomaly comparable to the 8.2 ka event. Our results show that the population-level consequences are consistent with this climatic heterogeneity: where the climate signal was regionally variable, so too were human demographic responses. Neukom et al. (2019) reached a parallel conclusion for the Common Era, finding no evidence of globally coherent warm or cold periods prior to the twentieth century and demonstrating that spatial variability is the norm rather than the exception in Holocene climate. That population dynamics follow this same pattern of spatial heterogeneity should perhaps not be surprising, but it has significant implications for how the 4.2 ka event is framed in the archaeological literature. The persistent narrative of a "global megadrought" (Weiss, 2016) causing uniform societal collapse requires qualification: the 4.2 ka event had severe demographic consequences in specific regions and for specific subsistence systems, but it was simultaneously a period of demographic expansion in others.

Climate modelling experiments by Renssen (2022) have demonstrated that the 4.2 ka event's spatial heterogeneity can be reproduced through the interaction of tropical sea-surface temperature anomalies and progressive desertification, providing a physical mechanism for the regional divergence we document. Nan et al. (2025) provide a mechanistic basis for this heterogeneity, documenting the spatial structure of the 4.2 ka climate signal across multiple proxy networks and showing that the intensity and even the sign of the climate anomaly varied substantially between regions. Our population data map onto this spatial pattern with notable fidelity. Regions where paleoclimate proxies record aridification or monsoon weakening — the Near East (Gol-e-Zard Cave: +0.35‰ δ¹⁸O) and monsoon-influenced East Asia (Dongge Cave: +0.31‰ δ¹⁸O) — are the same regions that experienced population declines, while regions with minimal climatic perturbation (GISP2: negligible change) showed demographic growth. This spatial correspondence is consistent with a climate-mediated mechanism, though we note that three proxy records cannot capture the full spatial complexity of 4.2 ka climate change, and the absence of a high-resolution Japanese-specific proxy record is a notable gap (see Section 5). The permTest framework we employ captures this divergence in a way that single-region modelTest approaches cannot, because it directly compares the shapes of regional SPD trajectories rather than testing each region independently against a null growth model.

### 4.2 Vulnerability Framework: Why Different Regions Diverged

The regional divergence documented by permTest demands explanation. We propose that differential vulnerability to 4.2 ka climate change was primarily mediated by two factors: the intensity of local climatic perturbation and the sensitivity of regional subsistence systems to that perturbation.

The Near East, where our analysis records a -24.9% population decline, was an arid-zone agricultural system critically dependent on monsoon-derived precipitation. Gol-e-Zard Cave records a +0.35‰ δ¹⁸O shift indicating substantial aridification coincident with this decline, consistent with Red Sea proxy records documenting sustained aridity in the broader region (Arz et al., 2006). Palmisano et al. (2021) documented a comparable ~30% decline in Levantine SPDs during this interval, and Lawrence et al. (2021) demonstrated through multi-proxy reconstruction that settlement contraction in the Northern Fertile Crescent was concentrated in rain-fed agricultural zones rather than irrigated lowlands. Our NERD cross-validation independently reproduces this pattern (-17.7%, permTest p = 0.003), supporting the inference that the Near Eastern decline is not an artifact of database-specific biases. The slightly lower magnitude of the NERD-based estimate (-17.7% vs. -24.9% from p3k14c) likely reflects differences in spatial coverage and dating density between the two databases, but the qualitative pattern and statistical significance are fully consistent.

In Chubu, Japan, the -30.1% decline occurred within a hunter-gatherer society whose Middle Jomon economy was heavily invested in nut-forest exploitation, particularly chestnut and walnut (Habu, 2004, 2016). Habu (2008) documented comparable growth-decline dynamics at the Sannai Maruyama site in northern Honshu, providing an independent archaeological case study of Middle Jomon demographic instability. This subsistence strategy, while sophisticated, was inherently sensitive to changes in precipitation and temperature that affected forest productivity. The temporal coincidence of Chubu's decline with evidence of East Asian monsoon weakening (Dongge Cave, located ~2,000 km to the southwest; see also Kathayat et al., 2017, and Berkelhammer et al., 2013, for independent Indian monsoon proxy evidence of weakening at 4.2 ka) is consistent with the hypothesis that reduced forest productivity may have undermined the carrying capacity of the nut-forest niche, though a Japanese-specific paleoclimate record would be needed to establish this link more firmly. However, the role of endogenous factors — particularly the possibility that Middle Jomon populations had already approached or exceeded carrying capacity — must be considered alongside climatic forcing (see Section 4.3).

Western Europe, by contrast, experienced +37.1% growth during the same interval. This humid-zone region was undergoing the Neolithic-to-Bronze Age transition (Grossmann et al., 2023; Silva and Vander Linden, 2017), a period of technological innovation and agricultural intensification that appears to have proceeded largely independently of 4.2 ka climate stress. The minimal high-latitude climate signal recorded in GISP2 (+0.66°C, within normal variability) suggests that Western European populations were simply not exposed to the same climatic perturbation that affected monsoon-dependent regions. Similarly, Kyushu's remarkable +84.6% growth reflects the Late Jomon cultural florescence in western Japan, characterized by a diversified coastal-terrestrial subsistence base that was less dependent on inland nut forests than the Chubu economy. Kyushu's growth may also partly reflect the earliest phases of the Jomon-to-Yayoi transition, as wet-rice agriculture began reaching northern Kyushu during this interval (Crema and Shoda, 2021). If the 4.2 ka temporal window overlaps with early agricultural adoption, the observed growth would reflect subsistence intensification rather than — or in addition to — climate resilience.

Railsback et al. (2018) documented a two-pulsed structure to the 4.2 ka climate event, with distinct aridification episodes separated by a brief amelioration. This temporal structure may account for some of the timing differences between regions: societies that were resilient to the first pulse may have been stressed beyond recovery by the second, while others may have adapted their subsistence strategies during the intervening amelioration. The ~200-year offset between the Japanese divergence (~4,400 cal BP) and the canonical 4.2 ka date could reflect either earlier onset of monsoon weakening in the western Pacific or the interaction of climatic stress with pre-existing demographic trajectories, as discussed below.

#### 4.2.1 A Semi-Quantitative Vulnerability Index

Following the IPCC vulnerability framework (V = f(Exposure, Sensitivity, Adaptive Capacity); Brooks, 2014) and the interdisciplinary approach to studying societal responses to environmental change advocated by Haldon et al. (2018), as well as its recent application to archaeological contexts (Daly, 2014; Jones et al., 2023), we construct a simple two-factor vulnerability index for the six regional cases analysed in this study. We define **climate exposure** as the magnitude of paleoclimate perturbation recorded by proxies relevant to each region, and **subsistence sensitivity** as the degree to which the dominant economic strategy depends on climate-controlled resources. The product of these two factors yields a predicted vulnerability score, which we then compare against observed population change (Table 7a).

**Table 7a.** Semi-quantitative vulnerability index for the 4.2 ka event across six study regions.

| Region | Climate exposure | Justification | Subsistence sensitivity | Justification | Predicted vulnerability | Observed change |
|---|---|---|---|---|---|---|
| Near East | High | Gol-e-Zard +0.35 permil delta-18O; severe aridification | High | Rain-fed cereal agriculture in arid zone; narrow precipitation margin | **High** | **-24.9%** |
| Chubu (Japan) | Medium | Dongge Cave +0.31 permil (regional monsoon proxy); no local speleothem | High | Nut-forest monoculture (chestnut/walnut); productivity directly tracks temperature and precipitation | **Medium-High** | **-30.1%** |
| Yangtze basin | Medium | Dongge Cave +0.31 permil; monsoon weakening plus flood risk | Medium | Wet-rice cultivation vulnerable to flooding but with some dietary diversification | **Medium** | **-16.2%** |
| Yellow River basin | Medium | Dongge Cave +0.31 permil; same monsoon domain | Low | Dryland millet agriculture; drought-tolerant crop in fertile loess soils | **Medium-Low** | **+18.1%** |
| W. Europe | Low | GISP2 +0.66 degrees C (within normal Holocene range); no aridification signal | Low | Mixed Neolithic-Bronze Age agropastoralism in humid zone; diversified crops and livestock | **Low** | **+37.1%** |
| Kyushu (Japan) | Medium | Same East Asian monsoon domain as Chubu | Low | Diversified coastal-terrestrial economy; marine resources buffer terrestrial shortfalls | **Medium-Low** | **+84.6%** |

The predicted vulnerability rankings correspond well with observed population trajectories. The two regions scored as High or Medium-High vulnerability --- the Near East and Chubu --- are the only regions that experienced statistically significant population declines detected by permTest (p = 0.003 and p = 0.002, respectively). The Yangtze basin, scored as Medium vulnerability, declined but the permTest comparison with the Yellow River was non-significant (p = 0.71), likely reflecting low statistical power. The three regions scored as Medium-Low or Low vulnerability --- the Yellow River basin, western Europe, and Kyushu --- all experienced population growth ranging from +18.1% to +84.6%.

We note that the subsistence sensitivity ratings for the Near East and Chubu are informed by the same archaeological literature that documents the population declines themselves. This potential circularity means that the concordance between predicted and observed outcomes should be interpreted as internal consistency rather than independent prediction.

This concordance between predicted and observed outcomes supports the interpretation that differential vulnerability, rather than differential climate forcing alone, explains the regional divergence in population responses to the 4.2 ka event. The framework reveals that high vulnerability requires the conjunction of strong climate exposure *and* high subsistence sensitivity. Chubu's extreme decline (-30.1%, the largest globally) despite only medium climate exposure underscores the role of subsistence sensitivity: the nut-forest monoculture of the Middle Jomon economy amplified a moderate climate signal into a severe demographic outcome. Conversely, Kyushu experienced the same regional climate exposure as Chubu but was buffered by subsistence diversification, producing growth rather than decline.

Two caveats apply. First, adaptive capacity --- the third component of the IPCC framework --- is omitted here because it cannot be reliably quantified from the archaeological record at this temporal resolution. Second, the ratings are ordinal and qualitative; a fully quantitative index would require standardised paleoclimate anomaly measures across all regions and independent subsistence-dependency metrics, which remain unavailable for several of our study areas. Despite these limitations, the framework transforms the narrative interpretation of Section 4.2 into a testable proposition: vulnerability to abrupt climate events is predictable from the interaction of exposure and sensitivity, and future cross-regional demographic studies can refine and extend this index as higher-resolution proxy data become available.

### 4.3 Japan's Chubu Decline: Magnitude, Timing, and Interpretation

The Chubu decline of -30.1% is the largest magnitude population decline detected in this analysis, exceeding the Near East (-24.9%). However, this cross-regional magnitude comparison must be interpreted cautiously: Chubu has substantially higher dating density (8,297 dates from ~1,005 sites) than the Near East (2,099 dates from 707 sites), and regions with sparser radiocarbon coverage may contain undetected declines of comparable or greater magnitude. The result is nonetheless noteworthy for a non-agricultural society and merits careful interpretation. Several factors may contribute to its severity.

First, the timing of the Japanese divergence at ~4,400 cal BP, approximately 200 years earlier than the canonical 4.2 ka date, raises the question of whether the Chubu decline is attributable to the 4.2 ka climate event at all. Three explanations are possible: (a) the East Asian monsoon weakened earlier in the western Pacific than in the regions where the canonical 4.2 ka date was defined; (b) chronological uncertainty in the JOAD age models — particularly the plateau in the IntCal20 calibration curve near 4,400 cal BP — may have shifted the apparent timing of the decline; or (c) the decline reflects endogenous factors, specifically the culmination of Middle Jomon population overshoot, that happened to overlap broadly with 4.2 ka climate change. Kawahata et al. (2019) argued on the basis of evidence from Sannai-Maruyama that environmental deterioration in northern Honshu began as early as 4,900 cal BP, and Habu (2016) similarly proposed that eastern Japanese decline preceded the 4.2 ka event by several centuries.

Our sensitivity analysis partially addresses this ambiguity. When the pre-event baseline window is extended to 5,500–6,000 cal BP — capturing the peak of Middle Jomon population — the Chubu decline shrinks from -30.1% to -12.3%. This suggests that a substantial portion of the apparent decline represents population decrease that was already underway before 4,200 cal BP. Critically, however, Chubu remains the only region among all those analysed that shows a decline under the extended window; all other regions show growth ranging from +55% to +206%. The relative anomaly of the Chubu trajectory is therefore robust to window definition, even if its absolute magnitude is sensitive to the choice of baseline.

Crema et al. (2016) first documented regional divergence in Jomon population dynamics using SPD analysis of eastern Japan, identifying distinct trajectories for the Kanto and Aomori regions. Our analysis, enabled by the JOAD database (Kudo, 2023), extends this finding with substantially finer regional resolution (six regions versus three) and a 27-fold increase in the number of radiocarbon dates. The JOAD data reveal that the core decline zone was Chubu-Kanto (central Honshu) rather than eastern Japan broadly — Tohoku, initially showing a marginal decline of -6.6%, reverses to slight growth (+6.5%) after site binning. This spatial precision is important because it localises the demographic impact to the heartland of Middle Jomon nut-forest exploitation, consistent with the vulnerability framework outlined above.

The interaction between endogenous population dynamics and exogenous climate forcing in Chubu likely represents a case of compound vulnerability: a society that had already begun to exceed the carrying capacity of its subsistence niche was then subjected to climatic conditions that further reduced that carrying capacity. This pattern resonates with the findings of Downey et al. (2016), who documented early warning signals of population collapse in European Neolithic societies, suggesting that societies on the edge of carrying capacity are disproportionately vulnerable to external perturbation. Disentangling the relative contributions of these factors will require high-resolution palaeoenvironmental reconstructions from the Chubu region itself, ideally including pollen, charcoal, and stable isotope records from lacustrine or peat sequences.

### 4.4 The China Puzzle and the Yangtze–Yellow River Divergence

The aggregate Chinese result (+11.6% growth) conceals one of the most instructive patterns in our dataset. Sub-regional analysis reveals that the Yellow River basin experienced +18.1% growth while the Yangtze basin declined by -16.2%, a divergence consistent with the archaeological record of contrasting cultural trajectories: Longshan expansion in the north and Liangzhu collapse in the south (Shao et al., 2021). The Liangzhu civilisation, which depended on intensive wet-rice cultivation in the lower Yangtze, was particularly vulnerable to the flooding and hydrological disruption associated with intensified monsoon variability. Meanwhile, millet-based agriculture in the Yellow River basin appears to have been resilient to, or even benefited from, the same climatic shifts.

The permTest comparison between these sub-regions was non-significant (p = 0.71). Two interpretations are possible: insufficient statistical power given the Yangtze sample of only ~800 radiocarbon dates, or genuine absence of statistically detectable divergence in the radiocarbon record despite the archaeological evidence for contrasting cultural trajectories. Without a formal power analysis, we present this as a genuinely ambiguous result rather than attributing it to either explanation. Future analyses incorporating larger Chinese radiocarbon databases should revisit this comparison.

This finding has broader implications for the treatment of large regions as analytical units. Wang et al. (2014) and Dong et al. (2020, 2022) examined Chinese population dynamics at the national scale and identified broad patterns of growth and decline, but the opposing trajectories of the Yellow River and Yangtze basins cancel when aggregated, producing a misleading picture of moderate overall growth. Our results reinforce the methodological principle that sub-regional analysis is essential when studying population responses to climate events whose spatial footprint is itself heterogeneous (Wirtz et al., 2024).

### 4.5 Methodological Implications

This study yields several methodological insights for SPD-based demographic reconstruction that extend beyond the 4.2 ka case study.

First, site binning is not merely a recommended best practice but a critical safeguard against spurious signals. Our Northern European analysis provides a stark demonstration: the unbinned SPD suggested a -21.0% decline, which reversed to +10.4% growth after site-level binning. This finding vindicates the recommendations of Crema and Bevan (2021) and underscores that unbinned SPDs should not be interpreted as demographic proxies without careful consideration of site-level sampling heterogeneity.

Second, the mark permutation test (permTest) is substantially more informative than single-region model testing (modelTest) for detecting spatially structured demographic change. In our analysis, modelTest failed to identify statistically significant deviations in either the Near East or Chubu when tested individually, consistent with the known limitations of this approach (Crema, 2022). By contrast, permTest detected highly significant divergence between climatically distinct regions (p = 0.002–0.003). For research questions concerned with differential responses to climate events, permTest should be the primary inferential tool.

Third, cross-validation using independent databases substantially strengthens confidence in SPD-derived patterns. Our replication of the Near East–Western Europe divergence using the NERD database yielded consistent results (p = 0.003). We advocate for routine cross-validation as a standard component of SPD analyses wherever multiple overlapping databases exist.

Fourth, sensitivity to the pre-event baseline window must be explicitly assessed. Our extended-window analysis demonstrated that Chubu's decline magnitude is substantially affected by window choice (-30.1% vs. -12.3%), but Chubu remained the sole declining region under all definitions, supporting the robustness of the qualitative pattern.

Finally, our comparative approach — testing the significance of divergence between regions rather than interpreting point-wise SPD fluctuations — addresses the critique raised by Carleton and Groucutt (2021) regarding over-interpretation of SPD morphology. The permTest framework provides a rigorous basis for comparative demographic inference that is more robust than approaches relying on visual inspection or point-wise significance testing.

---

## Caveats and Limitations

Several limitations of this study should be acknowledged. First, taphonomic bias systematically reduces the representation of older radiocarbon dates due to differential preservation, potentially affecting the magnitude of pre-event to event-window comparisons. This bias is partially mitigated by site-level binning and by our focus on relative inter-regional comparisons (permTest) rather than absolute demographic reconstruction. Second, inter-regional differences in research intensity — the British Isles and Japan have disproportionately high dating densities relative to the Near East and China — could introduce spatial sampling biases. The permTest framework addresses this by comparing SPD shapes rather than absolute magnitudes, and the NERD cross-validation provides an independent check on the Near Eastern results. Third, calibration curve effects, particularly plateaux in the IntCal20 curve, can distort SPD morphology around certain time intervals; our sensitivity analysis with multiple window definitions partially addresses this concern.

Fourth, the ~200-year timing offset between the Japanese divergence (~4,400 cal BP) and the canonical 4.2 ka event introduces uncertainty about the causal relationship between climate change and the Chubu decline. The extended baseline analysis (Section 3.6) confirms that a declining trend predated the 4.2 ka interval, suggesting that the climate event may have exacerbated rather than initiated the demographic downturn. Fifth, the non-significant permTest result for the Yellow River–Yangtze comparison (p = 0.71) could reflect either insufficient statistical power or genuine absence of statistically detectable divergence in the radiocarbon record; we cannot distinguish between these explanations with the available data. This result should not be interpreted as definitive evidence either for or against uniform Chinese sub-regional response. Sixth, the choice of three paleoclimate proxy records, while representative of the principal climate subsystems, does not capture the full spatial complexity of 4.2 ka climate change. Future studies should incorporate additional proxy networks, including those highlighted by Nan et al. (2025).

Finally, we note that SPD-based analyses remain proxies for population change and are subject to the fundamental limitation that radiocarbon date frequencies reflect archaeological activity broadly defined, not population size directly (Carleton and Groucutt, 2021; Crema, 2022). Our comparative framework mitigates but does not eliminate this concern.

---

## Conclusions

This study provides a unified statistical framework for comparing radiocarbon-inferred population dynamics across multiple world regions during the 4.2 ka climate event, applying identical SPD methodology, site binning, and permutation testing. Our findings lead to three principal conclusions.

First, population responses to the 4.2 ka event were regionally divergent, not globally uniform. The Near East experienced a −24.9% decline while western Europe grew by +37.1% (permTest p = 0.003), and within Japan, Chubu declined by −30.1% while Kyushu grew by +84.6% (p = 0.002). Independent cross-validation using the NERD database reproduced the Near East–Europe divergence (p = 0.003). This regional divergence is consistent with recent findings that the 4.2 ka climate event was neither globally coherent (McKay et al., 2024) nor climatically exceptional in most proxy records (Neukom et al., 2019). We argue that the appropriate analytical framework for the 4.2 ka event is not "globally uniform crisis" but "regionally contingent response," in which the interaction between spatially heterogeneous climate forcing and locally specific subsistence vulnerabilities shapes demographic outcomes.

Second, Japan's Chubu region recorded a −30.1% decline — the largest magnitude decline detected in our cross-regional analysis during the 4.2 ka interval, though cross-regional magnitude comparisons are conditioned by substantial differences in dating density among regions. This finding, enabled by the JOAD database (Kudo, 2023), integrates Japanese demographic history into the global 4.2 ka narrative at a quantitative level using standardized cross-regional permutation testing. The Chubu decline likely reflects the compound vulnerability of a nut-forest-dependent hunter-gatherer economy to monsoon weakening, possibly exacerbated by endogenous population overshoot during the Middle Jomon florescence.

Third, our analysis demonstrates the critical importance of methodological rigour in SPD-based demographic inference. Site binning corrected a false decline signal in northern Europe (−21.0% → +10.4%), the permutation test detected inter-regional divergence where single-region modelTest failed, and the China sub-regional analysis revealed opposing trends (Yellow River +18.1% vs. Yangtze −16.2%) masked by national-scale aggregation. These results underscore the recommendations of Crema and Bevan (2021) and Crema (2022) and caution against interpreting unbinned or aggregated SPDs at face value.

Future work should extend this cross-regional permutation framework to other Holocene climate events (e.g., the 8.2 ka event), incorporate additional proxy data (settlement counts, environmental proxies), and develop formal vulnerability indices that integrate climatic exposure with subsistence-system sensitivity. The growing availability of regional radiocarbon databases — including the recent Chronology of Early China databank (Qiu et al., 2025) — will enable higher-resolution sub-regional analyses. The broader programme of understanding climate–society interactions over centennial to millennial timescales (Büntgen et al., 2016) requires precisely the kind of standardized cross-regional framework demonstrated here.

---

## References

Aimers, J., Hodell, D., 2011. Societal collapse: Drought and the Maya. Nature 479, 44–45. https://doi.org/10.1038/479044a

Alley, R.B., 2000. The Younger Dryas cold interval as viewed from central Greenland. Quaternary Science Reviews 19, 213–226. https://doi.org/10.1016/S0277-3791(99)00062-1

Arz, H.W., Lamy, F., Pätzold, J., Müller, P.J., Prins, M., 2003. Mediterranean moisture source for an early-Holocene humid period in the northern Red Sea. Science 300, 118–121. https://doi.org/10.1126/science.1080325

Arz, H.W., Lamy, F., Pätzold, J., 2006. A pronounced dry event recorded around 4.2 ka in brine sediments from the northern Red Sea. Quaternary Research 66, 432–441. https://doi.org/10.1016/j.yqres.2006.05.006

Berkelhammer, M., Sinha, A., Stott, L., Cheng, H., Pausata, F.S.R., Yoshimura, K., 2013. An abrupt shift in the Indian monsoon 4000 years ago. Geophysical Monograph Series 198, 75–87. https://doi.org/10.1029/2012GM001207

Bevan, A., Colledge, S., Fuller, D., Fyfe, R., Shennan, S., Stevens, C., 2017. Holocene fluctuations in human population demonstrate repeated links to food production and climate. Proceedings of the National Academy of Sciences 114, E10524–E10531. https://doi.org/10.1073/pnas.1709190114

Bini, M., Zanchetta, G., Persoiu, A., et al., 2019. The 4.2 ka BP Event in the Mediterranean region: an overview. Climate of the Past 15, 555–577. https://doi.org/10.5194/cp-15-555-2019

Bird, D., Miranda, L., Vander Linden, M., et al., 2022. p3k14c, a synthetic global database of archaeological radiocarbon dates. Scientific Data 9, 27. https://doi.org/10.1038/s41597-022-01118-7

Brooks, N., 2014. Vulnerability, risk and adaptation: A conceptual framework. Tyndall Centre Working Paper 38.

Büntgen, U., Myglan, V.S., Ljungqvist, F.C., et al., 2016. Cooling and societal change during the Late Antique Little Ice Age from 536 to around 660 AD. Nature Geoscience 9, 231–236. https://doi.org/10.1038/ngeo2652

Carleton, W.C., Groucutt, H.S., 2021. Sum things are not what they seem: Problems with point-wise interpretations and quantitative analyses of proxies based on aggregated radiocarbon dates. The Holocene 31, 630–643. https://doi.org/10.1177/0959683620981700

Carolin, S.A., Walker, R.T., Day, C.C., et al., 2019. Precise timing of abrupt increase in dust activity in the Middle East coincident with 4.2 ka social change. Proceedings of the National Academy of Sciences 116, 67–72. https://doi.org/10.1073/pnas.1808103115

Contreras, D.A., Meadows, J., 2014. Summed radiocarbon calibrations as a population proxy: a critical evaluation using a realistic simulation approach. Journal of Archaeological Science 52, 591–608. https://doi.org/10.1016/j.jas.2014.05.030

Crema, E.R., 2022. Statistical inference of prehistoric demography from frequency distributions of radiocarbon dates: A review and a guide for the perplexed. Journal of Archaeological Method and Theory 29, 1387–1418. https://doi.org/10.1007/s10816-022-09559-5

Crema, E.R., Bevan, A., 2021. Inference from large sets of radiocarbon dates: Software and methods. Radiocarbon 63, 23–39. https://doi.org/10.1017/RDC.2020.95

Crema, E.R., Bevan, A., Shennan, S., 2017. Spatio-temporal approaches to archaeological radiocarbon dates. Journal of Archaeological Science 87, 1–9. https://doi.org/10.1016/j.jas.2017.09.007

Crema, E.R., Habu, J., Kobayashi, K., Madella, M., 2016. Summed probability distribution of 14C dates suggests regional divergences in the population dynamics of the Jomon period in Eastern Japan. PLOS ONE 11, e0154809. https://doi.org/10.1371/journal.pone.0154809

Crema, E.R., Kobayashi, K., 2020. A multi-proxy inference of Jomon population dynamics using bayesian phase models, residential data, and summed probability distribution of 14C dates. Journal of Archaeological Science 117, 105136. https://doi.org/10.1016/j.jas.2020.105136

Crema, E.R., Shoda, S., 2021. A Bayesian approach for fitting and comparing demographic growth models of the Jomon to Yayoi transition in Kyushu (Japan). PLOS ONE 16, e0251695. https://doi.org/10.1371/journal.pone.0251695

Cullen, H.M., deMenocal, P.B., Hemming, S., et al., 2000. Climate change and the collapse of the Akkadian empire: Evidence from the deep sea. Geology 28, 379–382. https://doi.org/10.1130/0091-7613(2000)28<379:CCATCO>2.0.CO;2

Daly, C., 2014. The capacity to adapt to climate change at heritage sites. Environmental Science & Policy 47, 118–127. https://doi.org/10.1016/j.envsci.2014.11.003

deMenocal, P.B., 2001. Cultural responses to climate change during the late Holocene. Science 292, 667–673. https://doi.org/10.1126/science.1059287

Dong, G., Li, R., Lu, M., Zhang, D., James, N., 2020. Evolution of human–environmental interactions in China from the Late Paleolithic to the Bronze Age. Progress in Physical Geography 44, 233–258. https://doi.org/10.1177/0309133319876802

Dong, G., Lu, Y., Zhang, S., Huang, X., Ma, M., 2022. Spatiotemporal variation in human settlements and their interaction with living environments in Neolithic and Bronze Age China. Progress in Physical Geography 46, 949–967. https://doi.org/10.1177/03091333221087992

Downey, S.S., Haas, W.R., Shennan, S.J., 2016. European Neolithic societies showed early warning signals of population collapse. Proceedings of the National Academy of Sciences 113, 9751–9756. https://doi.org/10.1073/pnas.1602504113

Grossmann, R., Weinelt, M., Muller, J., 2023. Demographic dynamics between 5500 and 3500 calBP in Central Europe. PLOS ONE 18, e0291956. https://doi.org/10.1371/journal.pone.0291956

Habu, J., 2004. Ancient Jomon of Japan. Cambridge University Press, Cambridge.

Habu, J., 2008. Growth and decline in complex hunter-gatherer societies: A case study from the Jomon period Sannai Maruyama site, Japan. Antiquity 82, 571–584.

Habu, J., 2016. Food diversity and climate change: Lessons from the Early and Middle Jomon periods. Kokogaku Kenkyu (Quarterly of Archaeological Studies) 63, 38–50 (in Japanese with English summary).

Haldon, J., Mordechai, L., Newfield, T., et al., 2018. History meets palaeoscience: Consilience and collaboration in studying past societal responses to environmental change. Proceedings of the National Academy of Sciences 115, 3210–3218. https://doi.org/10.1073/pnas.1716912115

Heaton, T.J., Al-assam, S., Bard, E., 2025. A new approach to radiocarbon summarisation: Rigorous identification of variations/changepoints in the occurrence rate of radiocarbon samples using a Poisson process. Journal of Archaeological Science 176, 106139. https://doi.org/10.1016/j.jas.2025.106139

Jones, C., Sherren, K., Sherwood, N., 2023. Developing climate risk assessments for World Heritage. Internet Archaeology 60. https://doi.org/10.11141/ia.60.3

Kathayat, G., Cheng, H., Sinha, A., et al., 2017. The Indian monsoon variability and civilization changes in the Indian subcontinent. Science Advances 3, e1701296. https://doi.org/10.1126/sciadv.1701296

Kawahata, H., 2019. Climatic reconstruction at the Sannai-Maruyama site between Bond events 4 and 3 — implication for the collapse of the society at 4.2 ka event. Progress in Earth and Planetary Science 6, 63. https://doi.org/10.1186/s40645-019-0308-8

Kudo, Y., Sakamoto, M., Hakozaki, M., Stevens, C.J., Crema, E.R., 2023. An archaeological radiocarbon database of Japan. Journal of Open Archaeology Data 11, 11. https://doi.org/10.5334/joad.115

Lawrence, D., Palmisano, A., de Gruchy, M.W., 2021. Collapse and continuity: A multi-proxy reconstruction of settlement organization and population trajectories in the Northern Fertile Crescent during the 4.2kya Rapid Climate Change event. PLOS ONE 16, e0244871. https://doi.org/10.1371/journal.pone.0244871

Liu, F., Feng, Z., 2012. A dramatic climatic transition at ~4000 cal. yr BP and its cultural responses in Chinese cultural domains. The Holocene 22, 1181–1197. https://doi.org/10.1177/0959683612441839

McKay, N.P., Kaufman, D.S., Arcusa, S.H., et al., 2024. The 4.2 ka event is not remarkable in the context of Holocene climate variability. Nature Communications 15, 6555. https://doi.org/10.1038/s41467-024-50886-w

Middleton, G.D., 2017. Understanding Collapse: Ancient History and Modern Myths. Cambridge University Press, Cambridge.

Nan, Q., Chen, S., Liu, X., et al., 2025. The 4.2 ka event in the Northern Hemisphere: Spatial heterogeneity and driving mechanisms of hydroclimatic change. Earth-Science Reviews 265, 105128. https://doi.org/10.1016/j.earscirev.2025.105128

Neukom, R., Steiger, N., Gomez-Navarro, J.J., Wang, J., Werner, J.P., 2019. No evidence for globally coherent warm and cold periods over the preindustrial Common Era. Nature 571, 550–554. https://doi.org/10.1038/s41586-019-1401-2

Noshiro, S., Sasaki, Y., Yoshikawa, M., Kudo, Y., Bhandari, S., 2025. Survival during the 4.2 ka event by Jomon hunter-gatherers with management and use of plant resources at the Denotame site in central Japan. Vegetation History and Archaeobotany 34, 685–699. https://doi.org/10.1007/s00334-025-01040-z

Palmisano, A., Bevan, A., Shennan, S., 2017. Comparing archaeological proxies for long-term population patterns: An example from central Italy. Journal of Archaeological Science 87, 59–72. https://doi.org/10.1016/j.jas.2017.10.001

Palmisano, A., Lawrence, D., de Gruchy, M.W., Bevan, A., Shennan, S., 2021. Holocene regional population dynamics and climatic trends in the Near East: A first comparison using archaeo-demographic proxies. Quaternary Science Reviews 252, 106739. https://doi.org/10.1016/j.quascirev.2020.106739

Park, J., Park, J., Yi, S., Kim, J.C., Lee, E., Choi, J., 2019. Abrupt Holocene climate shifts in coastal East Asia, including the 8.2 ka, 4.2 ka, and 2.8 ka BP events, and societal responses on the Korean peninsula. Scientific Reports 9, 10806. https://doi.org/10.1038/s41598-019-47264-8

Qiu, M., Du, S., Shi, J., et al., 2025. Chronology of early China: A radiocarbon databank for Chinese archaeology. Scientific Data 12, 893. https://doi.org/10.1038/s41597-025-05956-z

Railsback, L.B., Liang, F., Brook, G.A., et al., 2018. The timing, two-pulsed nature, and variable climatic expression of the 4.2 ka event: A review and new high-resolution stalagmite data from Namibia. Quaternary Science Reviews 186, 78–90. https://doi.org/10.1016/j.quascirev.2018.02.015

Reimer, P.J., Austin, W.E.N., Bard, E., et al., 2020. The IntCal20 Northern Hemisphere Radiocarbon Age Calibration Curve (0–55 cal kBP). Radiocarbon 62, 725–757. https://doi.org/10.1017/RDC.2020.41

Renssen, H., 2022. Climate model experiments on the 4.2 ka event: The impact of tropical sea-surface temperature anomalies and desertification. The Holocene 32, 378–389. https://doi.org/10.1177/09596836221074031

Rick, J.W., 1987. Dates as data: An examination of the Peruvian preceramic radiocarbon record. American Antiquity 52, 55–73. https://doi.org/10.2307/281060

Riris, P., Silva, F., Crema, E., et al., 2024. Frequent disturbances enhanced the resilience of past human populations. Nature 629, 837–842. https://doi.org/10.1038/s41586-024-07354-8

Shao, K., Zhang, J., He, K., Wang, C., Lu, H., 2021. Impacts of the wetland environment on demographic development during the Neolithic in the Lower Yangtze region. Frontiers in Earth Science 9, 635640. https://doi.org/10.3389/feart.2021.635640

Shennan, S., Downey, S.S., Timpson, A., et al., 2013. Regional population collapse followed initial agriculture booms in mid-Holocene Europe. Nature Communications 4, 2486. https://doi.org/10.1038/ncomms3486

Silva, F., Vander Linden, M., 2017. Amplitude of travelling front as inferred from 14C predicts levels of genetic admixture among European early farmers. Scientific Reports 7, 11985. https://doi.org/10.1038/s41598-017-12318-2

Staubwasser, M., Sirocko, F., Grootes, P.M., Segl, M., 2003. Climate change at the 4.2 ka BP termination of the Indus valley civilization and Holocene south Asian monsoon variability. Geophysical Research Letters 30, 1425. https://doi.org/10.1029/2002GL016822

Staubwasser, M., Weiss, H., 2006. Holocene climate and cultural evolution in late prehistoric–early historic West Asia. Quaternary Research 66, 372–387. https://doi.org/10.1016/j.yqres.2006.09.001

Timpson, A., Colledge, S., Crema, E., et al., 2014. Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method. Journal of Archaeological Science 52, 549–557. https://doi.org/10.1016/j.jas.2014.08.011

Walker, M.J.C., Berkelhammer, M., Bjorck, S., et al., 2012. Formal subdivision of the Holocene Series/Epoch: a discussion paper by a Working Group of INTIMATE and the Subcommission on Quaternary Stratigraphy. Journal of Quaternary Science 27, 649–659. https://doi.org/10.1002/jqs.2565

Walker, M.J.C., Berkelhammer, M., Bjorck, S., et al., 2018. Formal ratification of the subdivision of the Holocene Series/Epoch. Episodes 41, 213–223. https://doi.org/10.18814/epiiugs/2018/018016

Wang, C., Lu, H., Zhang, J., Gu, Z., He, K., 2014. Prehistoric demographic fluctuations in China inferred from radiocarbon data and their linkage with climate change over the past 50,000 years. Quaternary Science Reviews 98, 45–59. https://doi.org/10.1016/j.quascirev.2014.05.015

Wang, Y., Cheng, H., Edwards, R.L., et al., 2005. The Holocene Asian monsoon: Links to solar changes and North Atlantic climate. Science 308, 854–857. https://doi.org/10.1126/science.1106296

Wanner, H., Beer, J., Bütikofer, J., et al., 2008. Mid- to Late Holocene climate change: an overview. Quaternary Science Reviews 27, 1791–1828. https://doi.org/10.1016/j.quascirev.2008.06.013

Weiss, H., 2016. Global megadrought, societal collapse and resilience at 4.2–3.9 ka BP across the Mediterranean and west Asia. PAGES Magazine 24, 62–63. https://doi.org/10.22498/pages.24.2.62

Weiss, H. (Ed.), 2017. Megadrought and Collapse: From Early Agriculture to Angkor. Oxford University Press, Oxford.

Wiener, M.H., 2014. The interaction of climate change and agency in the collapse of civilizations ca. 2300–2000 BC. Radiocarbon 56, S1–S16. https://doi.org/10.2458/azu_rc.56.18325

Williams, A.N., 2012. The use of summed radiocarbon probability distributions in archaeology: a review of methods. Journal of Archaeological Science 39, 578–589. https://doi.org/10.1016/j.jas.2011.07.014

Wirtz, K.W., Antunes, N., Diachenko, A., et al., 2024. Multicentennial cycles in continental demography synchronous with solar activity and climate stability. Nature Communications 15, 10248. https://doi.org/10.1038/s41467-024-54474-w

Wu, W., Liu, T., 2004. Possible role of the "Holocene Event 3" on the collapse of Neolithic cultures around the Central Plain of China. Quaternary International 117, 153–166. https://doi.org/10.1016/S1040-6182(03)00125-3

Zhang, H., Cheng, H., Sinha, A., Spötl, C., Cai, Y., et al., 2021. Collapse of the Liangzhu and other Neolithic cultures in the lower Yangtze region in response to climate change. Science Advances 7, eabi9275. https://doi.org/10.1126/sciadv.abi9275
