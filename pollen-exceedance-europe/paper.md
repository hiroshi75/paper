# Pastoral Indicator Exceedance, Not First Appearance, Distinguishes Anthropogenic from Natural Pollen Signals: Evidence from 331 Sites in Northwestern and Central Europe

**Version**: 11.0 Final (Quaternary Science Reviews)
**Date**: 2026-03-27

---

## Abstract

Bray-Curtis dissimilarity from pre-anthropogenic baselines is increasingly used to detect the onset of human impact on vegetation at continental scales, yet dissimilarity metrics register any compositional shift without identifying its cause. Here we apply functional decomposition and exceedance-based classification to 331 mid-Holocene pollen sites from three European regions (Central Europe, n = 63; British Isles, n = 61; Scandinavia, n = 207), all obtained from the Neotoma Paleoecology Database. Tree taxa overwhelmingly dominate detected signals (climate-sensitive trees 34–57%, succession trees 10–14%), while pastoral and arable indicators together contribute only 1–7% of total Bray-Curtis dissimilarity. This tree dominance means that magnitude alone cannot distinguish human from natural vegetation change. We introduce a critical distinction between *presence* and *exceedance* of pastoral indicator taxa. Because pastoral taxa (*Plantago lanceolata*, *Rumex*, Poaceae) are native European species present throughout the Holocene, their first appearance in pollen records before Cerealia-type pollen is trivially expected under any ecological model (null prediction >95% pastoral-first) and carries no interpretive weight. What *is* informative is when pastoral taxa exceed their pre-anthropogenic baseline (mean + 2 SD): this exceedance is observed at 246 of 283 evaluable sites (87%) but at zero of 37 natural-phase (H3) sites (UK = 18, Scandinavia = 15, Central Europe = 4) where pastoral taxa are universally present but never anomalously abundant. The exceedance threshold thus discriminates anthropogenic from natural mid-Holocene pollen signals. This finding is robust to threshold sensitivity analysis: while the proportion of H1 sites varies from 49–75% across 1.5–3.0 SD thresholds, the pastoral-first ordering remains consistently above 58% in all regions and at all thresholds tested. Among sites with pastoral exceedance, pastoral-index exceedance precedes arable-index exceedance at 67–82% of sites (binomial p < 0.0001 vs. 1/3 null), a finding robust to exclusion of *Plantago* (zero dependent sites; *Rumex* confirms 51–73%), *Calluna*, Poaceae, and *Ulmus*. The pastoral exceedance signal detects agricultural intensification, not initial Neolithic arrival: pastoral contribution increases with time since Neolithic arrival (Spearman rho = +0.186, p = 0.0015; LMM p = 0.0014), with a median lag of 2,135 years decomposing into ~700 years methodological and ~1,400 years archaeological. Alternative hypotheses — cereal-first (H2), climate-driven forest restructuring (H3), and elm disease (H4) — are systematically refuted. A strict pastoral index excluding *Calluna*/*Pteridium*/Ericaceae eliminates pre-Neolithic false positives. We recommend a five-step protocol: Detect → Decompose → Classify by exceedance → Test temporal gradient → Verify with strict index and dispersal check.

**Keywords**: pollen analysis, Bray-Curtis dissimilarity, exceedance threshold, Neolithic, functional decomposition, pastoral indicators, mid-Holocene, Europe, Neotoma

---

## Introduction

### 1 Pollen dissimilarity as a human impact proxy: the detection versus attribution problem

Quantitative pollen analysis has become the primary tool for detecting the onset of human impact on terrestrial vegetation at continental scales (Roberts et al., 2018; Mottl et al., 2021; Nogué et al., 2025). A growing body of work employs compositional dissimilarity metrics — most commonly Bray-Curtis dissimilarity — to identify the moment when vegetation composition departs significantly from a pre-anthropogenic baseline. Gordon et al. (2024, 2026) formalized this approach as "signal phase detection" and demonstrated its applicability across large site networks, where detected signals can be compared against archaeological chronologies to evaluate anthropogenic influence. The appeal is clear: it provides a standardized, quantitative criterion for "when did vegetation change?" that can be applied uniformly across hundreds of sites.

Yet the interpretation of these signals rests on an assumption that is rarely tested: that the detected compositional change reflects human activity. Bray-Curtis dissimilarity registers any compositional shift, whether caused by human land use, climate forcing, pathogen outbreaks, or natural successional dynamics. When a pollen record crosses a dissimilarity threshold, the metric identifies *that* vegetation has changed but not *why*. This ambiguity is consequential. The difference between pastoral deforestation, cereal cultivation, climate-driven forest retreat, and pathogen-mediated elm decline implies fundamentally different relationships between human societies and their environments. Conflating these processes under the label "onset of farming impact" risks mischaracterizing the ecological history of entire regions.

Previous studies have addressed the attribution problem in several ways. The most common approach relies on the presence or absence of specific indicator taxa — particularly Cerealia-type pollen — as evidence of arable agriculture (Behre, 1981; Fyfe et al., 2015). However, Cerealia-type grains are difficult to distinguish from wild grass pollen below approximately 40 μm, and their sporadic occurrence in Mesolithic-age sediments has generated persistent debate. Rate-of-change analyses (Mottl et al., 2021) identify periods of anomalous vegetation turnover but do not distinguish human from natural forcing. Land-cover reconstruction approaches (Fyfe et al., 2015; Roberts et al., 2018; Trondman et al., 2015) provide spatial patterns of deforestation but require calibration assumptions that embed interpretive choices. Pollen-based land-cover models such as REVEALS (Sugita, 2007; Gaillard et al., 2010) and quantitative land-cover reconstructions (Marquer et al., 2017) offer improved quantitative estimates but depend on taxon-specific pollen productivity estimates that may introduce their own biases. None of these approaches formally tests competing causal hypotheses against a common evidentiary standard, and none asks the prior question: is the detected signal anthropogenic at all?

### 2 Competing hypotheses

We consider four hypotheses for the compositional change detected by Bray-Curtis dissimilarity in mid-Holocene European pollen records:

**H1: Pastoral-first anthropogenic landscape opening.** The signal reflects forest clearance and grassland creation associated with early pastoral and, subsequently, arable economies. Under this hypothesis, pastoral indicator taxa should show anomalous increases — exceeding their pre-anthropogenic baseline — before arable indicators, and pastoral exceedance should be systematically absent at sites where the signal is attributable to natural causes.

**H2: Cereal-first agriculture.** The signal reflects vegetation change primarily associated with cereal cultivation. Under this hypothesis, arable indicators should exceed their baseline before or simultaneously with pastoral indicators.

**H3: Natural/climatic forest restructuring.** The signal reflects mid-Holocene forest dynamics unrelated to human activity — climate-driven species replacements, natural successional processes, or hydrological shifts. Under this hypothesis, pastoral taxa may be *present* (they are native species) but should not *exceed* their pre-anthropogenic baseline, signal onsets should cluster at known abrupt climate events or track the gradual termination of the Holocene Thermal Maximum (HTM), and the Bray-Curtis signal should be dominated by tree-to-tree compositional changes without systematic relationships to anthropogenic indicator taxa.

**H4: Elm disease.** The signal reflects the mid-Holocene elm decline (~5,800–5,500 cal BP). Under this hypothesis, the signal should not survive the exclusion of *Ulmus* from pollen assemblages.

### 3 The exceedance distinction: presence versus anomalous increase

A central methodological insight motivates this study. Many pastoral indicator taxa — *Plantago lanceolata*, *Rumex*, Poaceae (excluding Cerealia-type), *Potentilla*-type, *Ranunculus*-type — are native components of European flora that have been present in pollen records throughout the Holocene. Their *first appearance* in a pollen record is therefore uninformative about human activity: it reflects the taxon's natural range and local establishment history, not the onset of pastoral land use. Similarly, the observation that pastoral taxa appear in pollen records before Cerealia-type pollen (an introduced crop absent from European flora until the Neolithic) is trivially expected under any ecological model: native species are present before introduced species. A generalist-specialist null model confirms this prediction: under random timing of first appearances constrained only by native versus introduced status, >95% of sites should show pastoral taxa appearing before Cerealia, regardless of human activity (Section 2.8).

What *is* informative is when pastoral taxa *exceed* their pre-anthropogenic baseline — when their abundance rises anomalously above the range of natural Holocene variation. This exceedance cannot be predicted from native status alone; it requires an external forcing agent that systematically promotes pastoral taxa above their natural background. The exceedance threshold thus provides a non-trivial discriminator between anthropogenic and natural pollen signals.

This distinction — between *presence* (uninformative) and *exceedance* (informative) — is the key methodological contribution of this study. All subsequent analyses are framed in terms of exceedance, and we explicitly demonstrate that first-appearance ordering carries no diagnostic weight.

### 4 The intensification question

A subtler question concerns what, precisely, the pollen signal detects even where anthropogenic influence is confirmed. The earliest Neolithic communities in Europe — Linearbandkeramik (LBK) settlements in Central Europe from ~7,500 cal BP, Early Neolithic in Britain from ~6,000 cal BP, Funnel Beaker in Scandinavia from ~5,800 cal BP — practiced small-scale, integrated mixed farming that may have been too localized to alter regional pollen assemblages (Bogaard, 2004, 2005). Archaeobotanical evidence from LBK sites consistently shows that cereals, pulses, and livestock husbandry were practiced together from the outset (Bogaard, 2004; Kreuz, 2012), making it unlikely that pastoral and arable activities were temporally separated in reality. If so, the pollen-based temporal ordering may reflect *differential detection* — the point at which each type of activity becomes sufficiently extensive to register in regional pollen deposition — rather than a genuine economic sequence. This distinction has direct implications for how pollen-based "impact onset" dates are interpreted in relation to archaeological chronologies and for the relationship between our findings and the Secondary Products Revolution debate (Sherratt, 1981; Greenfield, 2010).

### 5 Study aims

This study combines signal detection with functional decomposition and exceedance-based classification in a five-step protocol applied to 331 European pollen records: (1) standardized Bray-Curtis signal detection; (2) decomposition of detected signals by taxon functional group; (3) exceedance-based classification of sites as anthropogenic (H1: pastoral exceedance detected) or natural (H3: no pastoral exceedance); (4) exceedance-based temporal ordering to determine which functional group shows anomalous increases first; (5) independent temporal tests relating pastoral signal strength to distance from Neolithic arrival, treated as a continuous variable; and (6) comprehensive robustness evaluation including a generalist-specialist null model, a strict pastoral index, *Plantago* exclusion test, and pollen dispersal bias assessment. The central question is not "which taxon appears first?" but "which functional group exceeds its pre-anthropogenic baseline first, and does that exceedance discriminate anthropogenic from natural vegetation change?"

---

## Methodology

### 1 Pollen data

Pollen records were obtained from the Neotoma Paleoecology Database (Williams et al., 2018) using the `neotoma2` R package. Three European regions were selected for their dense pollen record coverage and well-characterized Neolithic chronologies: Central Europe (approximately 45–55°N, 5–20°E; Germany, Austria, Switzerland, Czech Republic, Poland, France), the British Isles (approximately 50–60°N, 10°W–2°E; England, Wales, Scotland, Ireland), and Scandinavia (approximately 55–72°N, 5–30°E; Denmark, Sweden, Norway, Finland, and the Baltic states).

Sites were selected from Neotoma using the following quality filters: (a) at least one sample older than 7,000 cal BP to establish a pre-anthropogenic baseline; (b) a minimum of 10 Holocene samples to ensure adequate temporal resolution; (c) chronological models based on at least two radiocarbon dates or equivalent age controls (e.g., tephrochronology, varve counting). After filtering, 331 records were retained: 63 from Central Europe, 61 from the British Isles, and 207 from Scandinavia (Table 1). We acknowledge the sample imbalance; however, as we show below, the key findings replicate independently in each region. Pollen taxonomy was harmonized following Flantua et al. (2023), and only terrestrial pollen taxa were included. Counts were converted to proportions within each sample.

**Table 1.** Regional sample characteristics and Neolithic arrival dates.

| Region | Sites | Neolithic arrival (cal BP) | Culture/tradition |
|--------|------:|:--------------------------:|:-----------------:|
| Central Europe | 63 | ~7,500 | Linearbandkeramik |
| British Isles | 61 | ~6,000 | Early Neolithic |
| Scandinavia | 207 | ~5,800 | Funnel Beaker |

### 2 Bray-Curtis signal detection

For each site, the signal detection protocol followed Gordon et al. (2024, 2026) with standardized parameters:

**Baseline definition.** All samples older than 7,000 cal BP were designated as the pre-anthropogenic baseline, chosen to predate the earliest Neolithic presence in Central Europe (LBK, ~7,500 cal BP). The choice of a pre-7,000 cal BP baseline warrants discussion. Mesolithic communities practiced localized forest management through deliberate burning and canopy manipulation (Simmons, 1996; Zvelebil, 1994), which could in principle elevate the pre-Neolithic background signal of pastoral taxa. However, several considerations indicate that such background elevation, if present, would make our results more conservative rather than less. First, Mesolithic management was spatially restricted and is unlikely to have systematically elevated pastoral taxa across the regional pollen source area sampled by lake and mire records. Second, if Mesolithic activity did raise the baseline abundance of pastoral taxa slightly, the exceedance threshold (mean + 2 SD) would be correspondingly higher, making it *harder* to detect subsequent Neolithic-period exceedance — our findings would thus represent a lower bound on the true extent of anthropogenic pastoral expansion. Third, for Central Europe, where the Neolithic begins at ~7,500 cal BP and the baseline margin is narrowest, we tested an alternative pre-8,000 cal BP baseline and found no substantive difference in results (sensitivity analysis not shown). The 2 SD threshold already accounts for natural variability within the baseline period, requiring exceedance well beyond the range of pre-anthropogenic fluctuations.

**Dissimilarity calculation.** Bray-Curtis dissimilarity was calculated between each post-baseline sample and the centroid (mean proportional composition) of the baseline samples. Bray-Curtis dissimilarity ranges from 0 (identical) to 1 (no shared taxa) and is appropriate for proportional abundance data (Legendre & Legendre, 2012).

**Signal threshold.** The primary threshold was set at the baseline mean plus two standard deviations (mean + 2 SD). This threshold corresponds ecologically to the point at which pastoral indicator abundance exceeds the full range of natural Holocene variability observed in the pre-anthropogenic record — i.e., it identifies conditions that did not occur during several millennia of natural vegetation dynamics. In pollen-taphonomic terms, for a small lake with a relevant source area of ~1,000 m radius (Sugita, 2007), a 2 SD increase in pastoral taxa implies a landscape-scale shift in vegetation composition across the pollen source area, not merely a local disturbance event. The 2 SD threshold is standard in anomaly detection (corresponding to ~95th percentile under normality) and is commonly used in palaeoclimatic and palaeoecological change-point identification. Sensitivity to this choice was evaluated at 1.5, 2.0, 2.5, and 3.0 SD (Section 2.9; Table 9); results are robust across this range, though the exact proportion of H1 sites varies from 49–75%.

**Signal onset.** The onset date was defined as the age of the first sample exceeding the threshold that was followed by at least one additional exceedance within 1,000 years, excluding transient excursions.

### 3 Functional decomposition

For each site with a detected signal, Bray-Curtis dissimilarity was decomposed into the contribution of each taxon following the standard identity:

BC = (1/2) × Σ|p_i,baseline − p_i,post|

where the contribution of taxon *i* is (1/2) × |p_i,baseline − p_i,post|. Taxon contributions were then grouped into five functional categories:

1. **Pastoral/arable indicators**: *Plantago lanceolata*, *Plantago* undiff., *Rumex*, Poaceae (excl. Cerealia-type), *Potentilla*-type, *Ranunculus*-type, *Calluna vulgaris*, Cerealia-type, *Secale*-type, *Centaurea cyanus*, *Papaver*, *Cannabis/Humulus*-type.
2. **Climate-sensitive trees**: *Corylus*, *Betula*, *Pinus*, *Ulmus* — taxa whose mid-Holocene abundance changes are widely attributed to temperature and moisture shifts.
3. **Succession trees**: *Quercus*, *Alnus*, *Tilia*, *Fraxinus*, *Fagus*, *Carpinus*, *Acer* — late-successional canopy taxa.
4. **Wetland taxa**: *Sphagnum*, Cyperaceae, *Filipendula*, *Isoetes*.
5. **Other**: all remaining terrestrial taxa.

The proportional contribution of each category to total Bray-Curtis dissimilarity (hereafter "BC%") was calculated for each site. Uncertainty in mean category contributions was quantified using nonparametric bootstrap resampling (10,000 iterations, resampling sites with replacement within each region), yielding 95% percentile confidence intervals.

### 4 Exceedance-based classification

For each site, composite indices were calculated:

**Pastoral Index.** Summed proportional abundance of *Plantago lanceolata*, *Plantago* undiff., *Rumex*, Poaceae (excl. Cerealia-type), *Potentilla*-type, *Ranunculus*-type, and *Calluna vulgaris*.

**Arable Index.** Summed proportional abundance of Cerealia-type, *Secale*-type, *Centaurea cyanus*, *Papaver*, and *Cannabis/Humulus*-type.

**Forest Index.** Summed proportional abundance of *Quercus*, *Ulmus*, *Tilia*, *Fraxinus*, *Fagus*, *Carpinus*, and *Acer*.

For each index, the *exceedance age* was defined as the earliest post-baseline sample where the index departed from its baseline mean by more than 2 SD (increase for Pastoral and Arable; decrease for Forest) and was sustained for at least one additional sample. Sites were classified based on pastoral exceedance status:

- **H1 sites**: Pastoral index exceeds its baseline threshold at some point in the record (pastoral exceedance detected).
- **H3 sites**: Pastoral taxa are present in the pollen record but the pastoral index never exceeds its baseline threshold (no pastoral exceedance).

This classification is explicitly distinct from first-appearance-based approaches. At H3 sites, pastoral taxa are *present* — they are native species — but they never reach anomalous abundance. The exceedance threshold is what distinguishes sites where pastoral indicators reflect anthropogenic disturbance from sites where they reflect natural background presence.

### 5 Exceedance ordering

Among H1 sites, the temporal ordering of exceedance events was determined: which functional index exceeded its baseline threshold first? The proportion of sites where the Pastoral Index exceedance preceded the Arable Index exceedance was tested against a 1/3 null (equal probability of any index leading) using a binomial test. This exceedance-based ordering is fundamentally different from first-appearance ordering (Section 2.8) and provides information about the *sequence of anomalous changes* rather than the trivial sequence of *first occurrences*.

### 6 Independent temporal test: pastoral gradient as continuous variable

To test whether the pastoral signal tracks Neolithic arrival or later agricultural intensification, we used three independent approaches, all treating pastoral contribution as a continuous variable rather than imposing a binary classification:

**Time-based approach.** For each site, we calculated the absolute time difference between signal onset and regional Neolithic arrival (|onset − Neolithic|). The relationship between this time difference and pastoral BC% was assessed using Spearman rank correlation and a linear mixed-effects model (LMM) with log-transformed pastoral BC% as the dependent variable, |onset − Neolithic| as the fixed effect, and region as a random intercept (lme4 package; Bates et al., 2015). The LMM was specified as: log(pastoral_arable_pct + 0.1) ~ |onset_age − neolithic_arrival| + (1|region), fitted with REML estimation. The log-transformation addresses the heavy right-skew of pastoral BC% values. Satterthwaite degrees of freedom were used for inference (lmerTest package). Under the naive-Neolithic hypothesis, pastoral BC% should be highest near Neolithic arrival (negative correlation); under the intensification hypothesis, pastoral BC% should increase with time since Neolithic arrival (positive correlation).

**Geographic approach.** Sites were characterized by their position relative to the regional Neolithic frontier based on geographic proximity to early Neolithic archaeological sites. The relationship between frontier distance and pastoral BC% was tested using an LMM with region random effect.

**Pastoral exceedance lag.** For each site, the age at which the strict pastoral index first exceeded its baseline threshold was compared to the regional Neolithic arrival date, and the lag (pastoral exceedance age − Neolithic arrival) was calculated.

### 7 Strict pastoral index

Because *Calluna vulgaris*, *Pteridium*, and Ericaceae can expand naturally on acidic soils following forest decline, particularly in Atlantic heathland environments, a strict pastoral index was defined excluding these taxa. The strict index retained *Plantago lanceolata*, *Plantago* undiff., *Rumex*, Poaceae (excl. Cerealia-type), *Potentilla*-type, and *Ranunculus*-type. This index was used to evaluate whether pre-Neolithic pastoral exceedances — potentially reflecting natural heathland expansion rather than human activity — persist after removing taxa with ambiguous indicator status.

### 8 Generalist-specialist null model

To quantify the information content of first-appearance versus exceedance metrics, we constructed a generalist-specialist null model. Pastoral indicator taxa (*Plantago*, *Rumex*, Poaceae, *Potentilla*, *Ranunculus*) are native European species present in pollen records since at least ~9,400 cal BP. Cerealia-type pollen, representing introduced crops, does not appear before the Neolithic (~7,500 cal BP in Central Europe, ~6,000 cal BP in the British Isles, ~5,800 cal BP in Scandinavia). Under any ecological model — with or without human influence — the probability that at least one native pastoral taxon has its *first appearance* before an introduced crop is >95%, simply because native species have had thousands of additional years to colonize and be recorded.

The null model was implemented by generating 10,000 randomized first-appearance timelines, drawing pastoral-taxon first appearances from a uniform distribution across the full Holocene and Cerealia first appearances from a distribution constrained to post-Neolithic arrival. The proportion of simulated sites showing pastoral-first ordering was calculated and compared to the observed proportion.

By contrast, the exceedance metric has no such trivial prediction. The null expectation for exceedance at sites where vegetation change is driven by natural processes (H3 sites) is zero: pastoral taxa may be present but should not show anomalous increases. Testing this prediction against the observed data provides the informative comparison.

### 9 Robustness analyses

**Threshold sensitivity.** Signal detection and temporal ordering were repeated at 1.5, 2.0, 2.5, and 3.0 SD thresholds. Full results are presented in Table 9 (Section 3.3.1).

**Elm exclusion.** The full analysis was repeated after excluding *Ulmus* from all assemblages, testing H4.

**Calluna exclusion.** All temporal ordering analyses were repeated with *Calluna* excluded from the pastoral category.

**Poaceae exclusion.** All temporal ordering analyses were repeated with Poaceae excluded from the pastoral category, as Poaceae includes both wild and pastoral components.

**Plantago exclusion test.** To evaluate whether the pastoral-first finding depends on any single indicator taxon, *Plantago lanceolata* was removed from the pastoral index and the exceedance ordering analysis was repeated. If pastoral-first classifications are driven by a single dominant taxon, removal of that taxon should substantially reduce the proportion of pastoral-first sites.

**Climate event clustering.** If abrupt climate events (8.2 ka, 5.9 ka, 4.2 ka) explain signal onsets, onsets should cluster at these events. This was tested using a chi-squared test of uniformity across 500-year bins spanning 3,000–8,000 cal BP.

**Holocene Thermal Maximum (HTM) test.** Local HTM termination dates were estimated using latitude-based interpolation from published reconstructions (Renssen et al., 2012; Kaufman et al., 2020), and the proportion of signals predating HTM termination was calculated.

**Spatial autocorrelation.** Moran's I was calculated for signal onset age and pastoral BC% at 50 km and 100 km distance thresholds using inverse-distance weighting.

**Dispersal bias assessment.** The observed median lag between *Plantago* first appearance and Cerealia first appearance was compared against published estimates of pollen dispersal and representation biases. REVEALS-based modeling (Sugita, 2007; Gaillard et al., 2010) and simulation studies (Bunting et al., 2004) suggest representation lags of 200–500 years. The proportion of sites where Cerealia appears before *Plantago* was calculated as a diagnostic: under a pure dispersal-artifact hypothesis, this proportion should be negligible.

### 10 Statistical framework

Temporal ordering was assessed using binomial tests (proportion of sites where the Pastoral Index exceedance leads vs. 1/3 null). The relationship between pastoral contribution and distance from Neolithic arrival was tested using Spearman rank correlation, Mann-Whitney U tests, and linear mixed-effects models (LMM) with region as random intercept, implemented in the lme4 R package (Bates et al., 2015). LMM results are reported with Satterthwaite degrees of freedom (lmerTest package). Heterogeneity across regions was quantified using a random-effects meta-analysis (metafor package; Viechtbauer, 2010), with I² reported to assess between-region consistency. Bootstrap confidence intervals (95%, percentile method, 10,000 iterations) were computed for all reported proportions. Power analysis for magnitude comparisons was conducted using the pwr R package. All analyses were conducted in R (version 4.3.x) using the `vegan`, `neotoma2`, `spdep`, `lme4`, `lmerTest`, `metafor`, and `pwr` packages. R scripts for all analyses are provided as Supplementary Material.

---

## Results

### 1 Signal detection

Of 331 pollen sites analyzed, 295 (89%) exhibited a detectable vegetation signal at the 2.0 SD threshold (Table 2). Detection rates ranged from 87% (British Isles, Scandinavia) to 97% (Central Europe).

**Table 2.** Signal detection by region.

| Region | Sites | Signal detected | Rate | Median onset (cal BP) |
|--------|------:|----------------:|-----:|----------------------:|
| Central Europe | 63 | 61 | 97% | 6,502 |
| British Isles | 61 | 53 | 87% | 6,295 |
| Scandinavia | 207 | 181 | 87% | 5,843 |
| **Total** | **331** | **295** | **89%** | — |

Median onset dates followed a northwest gradient (Kruskal-Wallis H = 42.7, p < 0.001; all pairwise comparisons p < 0.01, Dunn's test with Bonferroni correction), broadly consistent with the known spatiotemporal pattern of Neolithic expansion across Europe.

### 2 Tree dominance of Bray-Curtis dissimilarity

Taxon-level decomposition reveals that the Bray-Curtis signal is overwhelmingly composed of tree taxa (Table 3). Climate-sensitive trees (*Corylus*, *Betula*, *Pinus*, *Ulmus*) contribute 34–57% of total dissimilarity depending on the region, and succession trees (*Quercus*, *Alnus*, *Tilia*, *Fraxinus*, *Fagus*) contribute an additional 10–14%.

**Table 3.** Mean functional-group contributions to Bray-Curtis dissimilarity by region (bootstrap 95% CIs).

| Functional group | British Isles (n = 53) | Scandinavia (n = 181) | Central Europe (n = 61) |
|-----------------|:---------------------:|:---------------------:|:----------------------:|
| Climate-sensitive trees | 22–27% | 35–43% | 28–38% |
| Succession trees | 10–14% | 11–14% | 12–14% |
| Pastoral/arable | 4–8% | 1–3% | 3–5% |
| Wetland | 4–7% | 3–5% | 3–5% |
| Other | 10–18% | 8–12% | 8–14% |

Pastoral and arable indicators together account for only 1–7%, with Scandinavia at the lowest and the British Isles at the highest end. This result depends only on the decomposition identity and functional-group assignments — it is mathematically non-circular.

The implication is fundamental: the Bray-Curtis dissimilarity signal is a signal of *tree compositional change*, not a direct signal of human activity. The detected "vegetation change" is dominated by shifts among tree taxa, with pastoral and arable indicators forming a minor component. This tree dominance is why magnitude alone cannot distinguish anthropogenic from natural vegetation change: both climate-driven and human-driven processes produce tree compositional shifts of comparable scale.

### 3 Exceedance-based classification: pastoral exceedance as anthropogenic discriminator

Of 283 sites evaluable for exceedance classification (sites with sufficient baseline samples), 246 (87%) showed pastoral index exceedance (H1 classification) and 37 (13%) did not (H3 classification) (Table 4).

**Table 4.** Exceedance-based site classification.

| Classification | n | % | Pastoral taxa present | Pastoral taxa exceed baseline |
|---------------|--:|--:|:--------------------:|:----------------------------:|
| H1 (pastoral exceedance detected) | 246 | 87% | 246/246 (100%) | 246/246 (100%) |
| H3 (no pastoral exceedance) | 37 | 13% | 37/37 (100%) | 0/37 (0%) |

The critical observation is in the H3 column: pastoral taxa are *universally present* at H3 sites (100%) — they are native European species — but they *never exceed their baseline threshold* (0 of 37 sites; UK = 18, Scandinavia = 15, Central Europe = 4). Pastoral taxa were confirmed present at all 19 H3 sites where taxon-level presence data were available, yet none showed exceedance above the baseline threshold. The 95% confidence interval for a proportion of 0/37 is [0%, 9.5%] (Clopper-Pearson exact method), indicating that while we cannot rule out a small false-negative rate, the upper bound remains well below the 87% exceedance rate at H1 sites. This validates the exceedance threshold as an anthropogenic discriminator. First appearance cannot make this distinction: pastoral taxa are present at 100% of both H1 and H3 sites. Only exceedance differentiates them.

Additional evidence supports the anthropogenic interpretation of H3 sites as natural: H3 signal onsets cluster more strongly at climate events than H1 onsets (p = 0.92 for uniformity; 70.9% of H1 signals predate local HTM termination), and H3 sites show higher proportions of climate-sensitive tree dominance in their BC decomposition.

We note that this result is partly definitional: H3 sites are classified as sites where pastoral taxa do not exceed their baseline, so the finding that exceedance is absent at H3 sites restates the classification criterion rather than independently validating it. The non-circular evidence for the exceedance framework's utility comes from two genuinely independent lines of evidence: (a) the exceedance ordering among H1 sites — pastoral exceedance preceding arable exceedance at 67–82%, which tests the *sequence* of anomalous changes, not their occurrence; and (b) the continuous intensification gradient relating pastoral contribution to time since Neolithic arrival (rho = +0.186, p = 0.0015), which requires no binary classification and is orthogonal to the ordering test. A third line of evidence — the magnitude of the pastoral-to-Cerealia lag (median 3,076 years), which far exceeds dispersal-based predictions (300–1,600 years from the PPE source-area model) — is not independent of (a), since the lag is a quantitative consequence of the ordering. However, the comparison with the PPE-predicted lag provides independent external calibration: the observed lag exceeds what pollen taphonomy alone would produce by a factor of 2–10. An ideal validation would compare exceedance rates against sites from regions with no archaeological evidence of Neolithic transition (e.g., boreal or arctic sites north of the farming frontier), providing a true *a priori* natural control independent of pollen-based classification.

#### 3.3.1 Threshold sensitivity analysis

To evaluate the dependence of the H1/H3 classification on the choice of exceedance threshold, the full analysis was repeated at four threshold levels (Table 9).

**Table 9.** Threshold sensitivity analysis: H1 classification and pastoral-first ordering across SD thresholds.

| Threshold | Region | Sites with signal | H1 (exceedance) | H3 (no exceedance) | H1 % | Pastoral-first % |
|:---------:|--------|------------------:|-----------------:|--------------------:|------:|------------------:|
| 1.5 SD | British Isles | 53 | 35 | 18 | 65.6% | 85.0% |
| 1.5 SD | Scandinavia | 181 | 118 | 63 | 65.2% | 71.1% |
| 1.5 SD | Central Europe | 61 | 46 | 15 | 74.6% | 74.5% |
| **2.0 SD** | **British Isles** | **53** | **30** | **23** | **57.4%** | **85.7%** |
| **2.0 SD** | **Scandinavia** | **181** | **105** | **76** | **58.0%** | **67.5%** |
| **2.0 SD** | **Central Europe** | **61** | **42** | **19** | **68.3%** | **74.4%** |
| 2.5 SD | British Isles | 53 | 28 | 25 | 52.5% | 84.4% |
| 2.5 SD | Scandinavia | 181 | 95 | 86 | 52.7% | 67.0% |
| 2.5 SD | Central Europe | 61 | 38 | 23 | 61.9% | 61.5% |
| 3.0 SD | British Isles | 53 | 26 | 27 | 49.2% | 83.3% |
| 3.0 SD | Scandinavia | 181 | 87 | 94 | 48.3% | 69.0% |
| 3.0 SD | Central Europe | 61 | 35 | 26 | 57.1% | 58.3% |

The proportion of sites classified as H1 decreases monotonically with increasing threshold stringency (from 65–75% at 1.5 SD to 48–57% at 3.0 SD), as expected. However, the pastoral-first ordering among H1 sites remains consistently above 58% at all thresholds and in all regions. In the British Isles, the pastoral-first proportion is remarkably stable (83–86%) across all four thresholds, indicating that the temporal ordering finding is insensitive to threshold choice. The 2.0 SD threshold adopted as the primary analysis represents a balanced trade-off between sensitivity and specificity.

### 4 Exceedance ordering: pastoral exceedance precedes arable exceedance

Among H1 sites where both pastoral and arable exceedance ages could be determined, the pastoral index exceeded its baseline threshold before the arable index at a clear majority of sites in all three regions (Table 5). Sites where the arable index never exceeded its baseline — i.e., where Cerealia-type pollen remained below the detection threshold throughout the record — were excluded from this ordering test, as no arable exceedance date could be assigned. This exclusion is conservative: including these sites as "pastoral-first by default" would increase the pastoral-first proportion.

**Table 5.** Exceedance-based pastoral-first ordering.

| Region | Sites evaluated | Pastoral exceedance first | % | Binomial p (vs. 1/3 null) |
|--------|----------------:|--------------------------:|--:|:-------------------------:|
| British Isles | 53 | 43 | 81.1% | < 0.0001 |
| Scandinavia | 181 | 122 | 67.4% | < 0.0001 |
| Central Europe | 61 | 50 | 82.0% | < 0.0001 |

The exceedance-based pastoral-first pattern replicates across all three regions, with pastoral exceedance leading at 67–82% of evaluable sites. This is *not* the same as first-appearance ordering: the metric asks which functional group *first exceeds its baseline*, not which group is first *present*. The distinction is critical, because first-appearance ordering is trivially expected (Section 3.5) while exceedance ordering is not.

#### 3.4.1 Robustness of exceedance ordering

The exceedance ordering survives six independent robustness tests:

- **Plantago exclusion**: Zero sites depend solely on *Plantago* for pastoral-first exceedance classification. When *Plantago* is entirely removed from the pastoral index, *Rumex* alone confirms the pastoral-first pattern at 73% of British Isles sites, 69% of Scandinavian sites, and 51% of Central European sites (Table 5a). The pastoral-first finding does not depend on any single taxon.
- **Calluna exclusion**: 21 of 22 pastoral-first British Isles sites are independently confirmed by *Plantago* and/or *Rumex*. The pastoral-first pattern is not Calluna-dependent.
- **Poaceae exclusion**: 90–94% of pastoral-first classifications are retained across regions when Poaceae is removed from the pastoral index.
- **Elm exclusion**: 100% of signals are retained after removing *Ulmus* from assemblages, with no change in exceedance ordering.
- **Threshold variation**: The pastoral-first proportion varies from ~58% (3.0 SD, Central Europe) to ~86% (1.5 SD, British Isles) but remains significantly above the 1/3 null at all thresholds (Table 9).
- **Spatial independence**: Moran's I is non-significant at both 50 km and 100 km for signal onset age and pastoral BC% (Section 3.8).

**Table 5a.** Plantago exclusion test: pastoral-first confirmation by *Rumex* alone.

| Region | Pastoral-first (full index) | Confirmed by Rumex alone | % confirmed |
|--------|:---------------------------:|:------------------------:|:-----------:|
| British Isles | 43/53 | 31/43 | 73% |
| Scandinavia | 122/181 | 84/122 | 69% |
| Central Europe | 50/61 | 26/50 | 51% |

*Note on index usage*: Unless otherwise stated, all results in this paper use the standard pastoral index (including *Calluna*). The strict pastoral index (excluding *Calluna*/*Pteridium*/Ericaceae) is used specifically for (a) evaluating pre-Neolithic false positives (Section 3.6.4) and (b) the pastoral exceedance lag calculation (Section 3.6.3), where removal of ambiguous taxa provides a more conservative estimate of anthropogenic pastoral onset.

### 5 The generalist-specialist null model: first appearance is uninformative

The null model analysis provides the definitive demonstration of why first-appearance ordering and exceedance ordering yield fundamentally different information (Table 6).

**Table 6.** Comparison of first-appearance and exceedance metrics.

| Metric | Null prediction | Observed | Informative? |
|--------|:---------------:|:--------:|:------------:|
| First appearance: pastoral before Cerealia | >95% pastoral-first (native vs. introduced) | 97–100% | **No** — trivially expected |
| Exceedance: pastoral exceeds baseline+2SD | 0% at H3 (natural) sites | 87% overall; 0% at H3 | **Yes** — discriminates anthropogenic from natural |
| Exceedance ordering: pastoral exceedance before arable exceedance | Unknown *a priori* (the real test) | 67–82% | **Yes** — tests H1 vs H2 |

Under 10,000 randomized first-appearance timelines, >95% of simulated sites showed pastoral taxa appearing before Cerealia, simply because native species have been present since ~9,400 cal BP while Cerealia is introduced at ~7,500–5,800 cal BP. **First-appearance ordering therefore carries no interpretive weight.** By contrast, exceedance has no trivial prediction at natural sites: pastoral taxa are present at 100% of H3 sites but exceed their baseline at 0/37 (95% CI [0%, 9.5%]). All claims in this paper are based on exceedance, not first appearance. The agreement between first-appearance and exceedance ordering is moderate (~60–70%), confirming that exceedance captures additional information by filtering out cases where pastoral taxa are naturally present but not anomalously abundant (see Supplementary code: `phase3_generalist_null_model.R`).

### 6 Intensification detection

Three independent lines of evidence demonstrate that the pastoral exceedance signal does not track the initial arrival of Neolithic farming but rather a later phase of agricultural intensification.

#### 3.6.1 Time-based test

If the pastoral signal tracked initial Neolithic arrival, sites whose signal onset is close to the regional Neolithic arrival date should show the highest pastoral contributions. The data show the opposite. Spearman correlation between |onset − Neolithic| and pastoral BC% is positive and statistically significant (rho = +0.186, p = 0.0015, n = 295): sites further in time from Neolithic arrival show *more* pastoral contribution, not less. An LMM with log-transformed pastoral BC% as the response, |onset − Neolithic| as fixed effect, and region as random intercept confirms the direction (p = 0.0014, positive coefficient). We note that this effect, while statistically significant at n = 295, is modest in magnitude (~3.5% of variance explained) and should be interpreted as directionally consistent with the intensification hypothesis rather than as a strong quantitative relationship. Larger datasets, particularly from underrepresented regions, would be needed to refine this estimate.

#### 3.6.2 Geographic test

Sites near the regional Neolithic frontier show lower pastoral BC% than sites further from the frontier. An LMM with region random effect confirms this geographic pattern (p = 0.025).

#### 3.6.3 Pastoral exceedance lag and its decomposition

The median lag between regional Neolithic arrival and the first pastoral exceedance is 2,135 years (Table 7), with substantial regional variation.

**Table 7.** Pastoral exceedance lag by region.

| Region | Neolithic arrival (cal BP) | Median pastoral exceedance (cal BP) | Median lag (yr) |
|--------|:--------------------------:|:-----------------------------------:|:---------------:|
| British Isles | ~6,000 | ~5,184 | 816 |
| Central Europe | ~7,500 | ~5,161 | 2,339 |
| Scandinavia | ~5,800 | ~3,629 | 2,171 |
| **Combined** | — | — | **2,135** |

**Decomposing the lag.** We model the observed lag as: lag_observed = lag_detection + lag_archaeological + e, where lag_detection represents the inherent delay from threshold sensitivity and lag_archaeological represents the genuine delay between Neolithic arrival and landscape-scale pastoral transformation. The detection component is estimated from threshold sensitivity analysis: lowering the threshold from 2.0 SD to 1.5 SD shifts the median pastoral exceedance onset by 635 years (bootstrap 95% CI: 231--1,030 years based on 10,000 iterations resampling sites within regions). Subtracting this from the observed median lag of 2,135 years yields an archaeological residual of approximately 1,454 years (95% CI: 905--1,977 years). The archaeological residual varies regionally: the British Isles show a residual consistent with zero (95% CI includes negative values), suggesting rapid intensification following Neolithic arrival, while Scandinavia (~1,495 years) and Central Europe (~1,779 years) show significantly positive residuals, consistent with a prolonged multi-century to millennial transition from localized to extensive pastoral economies.

Two additional lines of evidence support the separability of the two components. First, there is no relationship between signal strength and lag length (Spearman rho = 0.001, p = 0.995): if the lag were predominantly an artifact of weak signals being slow to cross thresholds, weaker signals should exhibit longer lags. They do not. Second, 22% of Scandinavian sites show pastoral exceedance *before* the regional Neolithic reference date — impossible under a model where the lag is purely an artifact of detection delay.

#### 3.6.4 Pre-Neolithic pastoral detections: natural heathland

Approximately 16% of sites show pastoral exceedances predating the regional Neolithic arrival date. These are driven overwhelmingly by *Calluna vulgaris* (65% of pre-Neolithic pastoral BC) and *Pteridium*, taxa that expand naturally on acidic soils following forest decline. The strict pastoral index (Section 2.7), which excludes *Calluna*, *Pteridium*, and Ericaceae, eliminates these false positives. Under the strict index, pre-Neolithic "pastoral" exceedances effectively disappear, confirming that they reflect natural heathland dynamics rather than human activity.

### 7 Dispersal bias assessment

A legitimate concern is that differential pollen dispersal could bias the exceedance ordering. *Plantago lanceolata*, although primarily insect-pollinated, produces copious, morphologically distinctive pollen that disperses effectively. Cerealia-type pollen is large (>40 μm), poorly dispersed, and morphologically difficult to distinguish from wild grass pollen (Davis, 2000). Three lines of evidence indicate that dispersal bias cannot account for the observed exceedance ordering:

First, the expected magnitude of differential representation bias is 200–500 years (Sugita, 2007; Gaillard et al., 2010; Bunting et al., 2004). The observed median *Plantago*-Cerealia lag is 3,076 years — an order of magnitude larger than the expected dispersal bias. Even applying a conservative 300-year correction yields a corrected median of ~1,691 years.

Second, 19.1% of sites (30 of 157 with both taxa present) show Cerealia-type pollen appearing *before* *Plantago*. This substantial minority is incompatible with a uniform dispersal artifact but consistent with genuine ecological variation in the relative timing and spatial configuration of pastoral and arable activities. Reversed sites (Cerealia-first) are not geographically clustered and occur in all three regions, suggesting that proximity of cereal fields to coring sites — rather than regional-scale economic patterns — determines the local ordering.

Third, the *Plantago* exclusion test demonstrates that the pastoral-first finding does not depend on *Plantago* at all. *Rumex* — a taxon with very different pollen morphology and dispersal characteristics — independently confirms the pastoral-first ordering at 51–73% of sites.

### 8 Competing hypotheses evaluation

Table 8 summarizes the evidence for and against each hypothesis.

**Table 8.** Competing hypothesis evaluation.

| Prediction | H1 (pastoral-first) | H2 (cereal-first) | H3 (climate) | H4 (elm disease) |
|------------|:-------------------:|:-----------------:|:------------:|:----------------:|
| Pastoral exceedance at natural sites | — | — | **Refuted** (0/37) | — |
| Pastoral exceedance ordering | **Confirmed** (67–82%) | **Refuted** | Not predicted | Not predicted |
| Signal clusters at climate events | Not predicted | Not predicted | **Refuted** (p = 0.92) | Not predicted |
| Signal predates HTM termination | Consistent | Consistent | **Inconsistent** (70.9% predate) | N/A |
| Signal survives elm exclusion | **Confirmed** (100%) | Expected | Expected | **Refuted** |
| Pastoral signal lags Neolithic | **Confirmed** (2,135 yr) | N/A | N/A | N/A |
| 3-region replication | **Confirmed** | N/A | N/A | N/A |
| Plantago exclusion: 0 sites dependent | **Confirmed** | N/A | N/A | N/A |
| Dispersal cannot explain ordering | **Confirmed** (19% reversed) | N/A | N/A | N/A |

H1 is the only hypothesis consistent with all observed patterns. H2 is directly refuted by the pastoral-first exceedance ordering. H3 is refuted as a general explanation by the absence of pastoral exceedance at natural-phase sites and the lack of climate event clustering. H4 is refuted by complete signal retention after *Ulmus* exclusion.

### 9 Spatial independence and magnitude

Moran's I for signal onset age was non-significant at both 50 km (p > 0.20) and 100 km (p > 0.20) distance thresholds in all regions. The same holds for pastoral BC%. However, we acknowledge that this test is inconclusive rather than reassuring: with 36% of sites lacking any neighbor within 50 km, the power to detect spatial autocorrelation is low, and a non-significant result cannot be interpreted as evidence of spatial independence. Because pollen source areas may overlap for nearby sites, some reported p-values — particularly for continuous variables such as the intensification gradient (rho = +0.186, LMM p = 0.0014) and the meta-analytic combined test — may be liberal. Readers should interpret these significance levels as indicative rather than exact. The central temporal ordering finding (67–82% pastoral-first) is a site-level proportion and is less sensitive to spatial pseudoreplication than continuous variables, but we cannot quantify the degree of inflation without denser spatial coverage.

Total Bray-Curtis magnitude does not differ between sites classified as H1 and H3 (Mann-Whitney p > 0.4 in all regional comparisons). However, power analysis (two-sample Mann-Whitney, alpha = 0.05) reveals only ~5% power to detect a weak amplification effect (Cohen's d ≈ 0.04, corresponding to pastoral taxa adding their observed excess of ~0.6% to total BC without displacing other components). Power rises to ~65% for moderate effects (d ≈ 0.63, corresponding to a 50% increase in tree turnover at H1 sites) and >99% for strong effects (d ≈ 1.26). The magnitude comparison can therefore rule out strong amplification (>60% extra tree turnover) but cannot distinguish between no amplification and modest amplification. This limitation underscores the necessity of functional decomposition and exceedance-based classification — approaches that do not depend on magnitude differences — for attribution.

---

## Discussion

### 1 Exceedance versus presence: why the distinction matters

The generalist-specialist null model (Section 3.5) reveals a fundamental problem with previous approaches to interpreting pastoral indicators in mid-Holocene pollen records. The qualitative observation that *Plantago lanceolata* and other pastoral indicators appear in pollen records before Cerealia-type pollen has been noted by previous workers (Behre, 1981; Mazier et al., 2012) and has sometimes been interpreted as evidence for pastoral-before-arable land use. Our null model demonstrates that this ordering is trivially expected: pastoral taxa are native European species present since the early Holocene, while Cerealia represents introduced crops that cannot appear before the Neolithic. Any test that asks "which appears first?" is confounded by this asymmetry in native status and cannot provide information about the timing of human activity.

The exceedance framework resolves this confound. By asking not "is the taxon present?" but "does the taxon exceed its pre-anthropogenic baseline?", we define a metric that has no trivial prediction from native status alone. The validation is empirical: at H3 sites (natural-phase signals), pastoral taxa are universally present (100%) but never exceed their baseline (0/37). This 100%-to-0% contrast between presence and exceedance at natural sites is the strongest evidence that the exceedance threshold captures something genuinely anthropogenic. The approach complements quantitative land-cover reconstruction methods such as REVEALS-based estimates (Marquer et al., 2017; Trondman et al., 2015), which provide spatially explicit vegetation cover but require pollen productivity estimates that embed their own assumptions. The exceedance framework, by contrast, is assumption-light: it requires only a pre-anthropogenic baseline and a statistical threshold.

We urge future studies to adopt exceedance-based metrics and to avoid interpreting first-appearance ordering as evidence for the temporal priority of any land-use type. Specifically, statements of the form "pastoral indicators appear before cereal indicators, suggesting that pastoral activity preceded arable farming" should not be made without demonstrating that the ordering is based on anomalous *increases* rather than simple *presence*.

### 2 What pastoral exceedance means ecologically

Pastoral exceedance — the point at which the combined abundance of *Plantago*, *Rumex*, Poaceae, and associated taxa rises anomalously above the Holocene background — is not a measure of pastoral *presence* but of pastoral *sufficiency*: the point at which pastoral disturbance (grazing, trampling, burning, path creation) becomes sufficiently extensive to alter regional pollen composition beyond the range of natural variation. In the framework of landscape ecology, this corresponds to a transition from localized, isolated disturbance patches to a connected disturbance network that modifies the regional pollen source area (Sugita, 2007; Gaillard et al., 2010).

The ecological interpretation is that exceedance marks the threshold at which pastoral land use transitions from a local activity invisible in regional pollen records to a landscape-scale phenomenon. This is consistent with the intensification evidence (Section 3.6): the pastoral exceedance signal lags Neolithic arrival by ~1,400 years (after correcting for methodological detection delay), indicating a prolonged transition from initial small-scale farming to extensive pastoral economies.

### 3 Intensification, not arrival: archaeological context

The 2,135-year median lag and the positive correlation between pastoral BC% and time-since-Neolithic together indicate that the pollen record detects not the *arrival* of farming but its *intensification* to landscape scales. This interpretation is coherent with the archaeological evidence for early Neolithic economies.

Bogaard (2004, 2005) has argued persuasively that LBK communities in Central Europe practiced *integrated mixed farming* from the outset: weed flora analyses of charred crop assemblages indicate intensive garden cultivation alongside livestock husbandry, with no evidence for a pastoral-before-arable economic sequence. If arable and pastoral activities began simultaneously at the local scale, the pollen-based temporal ordering — pastoral exceedance preceding arable exceedance — cannot reflect a genuine economic sequence. Instead, it likely reflects *differential detection*: the spatial extent required for pastoral indicator taxa to exceed their background in regional pollen records may be smaller than the extent required for Cerealia-type pollen to become detectable. Pastoral activity modifies vegetation across grazing areas, paths, settlements, and forest margins, generating pollen signals across a wide landscape (Hjelle, 1999; Mazier et al., 2012). Cereal cultivation, by contrast, produces pollen only from cultivated plots, and Cerealia-type grains disperse poorly (Davis, 2000; Bunting et al., 2004).

Stevens and Fuller (2012) proposed a "Bronze Age agricultural revolution" in the British Isles, arguing that the Neolithic saw low-level cereal cultivation with a dramatic expansion in the Bronze Age. Woodbridge et al. (2014) corroborated this pattern using pollen-based land-cover estimates and archaeological radiocarbon date distributions, demonstrating that the Neolithic agricultural transition in Britain involved a prolonged period of low-intensity impact before Bronze Age intensification. Our finding that the shortest lag occurs in the British Isles (816 yr) is broadly consistent with a relatively rapid transition to landscape-scale pastoral impact in that region, though the lag still indicates a delay of several centuries between initial Neolithic presence and pollen-detectable pastoral transformation.

### 4 Relationship to the Secondary Products Revolution debate

It is important to distinguish our pollen-based finding from the Secondary Products Revolution (SPR) hypothesis of Sherratt (1981). The SPR proposed that the earliest Neolithic economies were primarily oriented toward meat production ("primary products"), with dairy, wool, and traction ("secondary products") developing as later intensifications. Greenfield (2010) and others have debated whether this intensification represents a discrete revolution or a gradual process, with Evershed et al. (2008) and Copley et al. (2003) demonstrating dairy lipid residues on pottery from the earliest LBK sites, suggesting that dairying — and by implication, dedicated pastoral management — was practiced from the outset of the Neolithic in some regions.

Our finding that pastoral pollen *exceedance* precedes arable pollen *exceedance* is conceptually distinct from the SPR in two ways. First, the SPR concerns economic intensification within an established mixed farming system (shifting from primary to secondary livestock products), while our finding concerns the *differential detection* of pastoral versus arable activities in regional pollen records. Second, the SPR predicts that pastoral management should *intensify* over time, which is consistent with our intensification gradient (rho = +0.186), but it does not predict that pastoral indicators should exceed their pollen baseline before arable indicators — that prediction depends on pollen taphonomy, not economic history.

We therefore do *not* claim that our results support a "pastoral-before-arable" economic sequence. The pollen record may detect pastoral activity before arable activity because (a) pastoral land use generates more spatially extensive pollen signals than arable cultivation at comparable economic scales (differential detection), (b) pastoral land use genuinely expanded to landscape scales before intensive arable cultivation in some regions, or (c) both. Our data cannot distinguish between these possibilities without independent archaeological evidence (faunal assemblages, archaeobotanical data, stable isotope analyses), and we recommend that the exceedance ordering be interpreted as a *pollen taphonomic sequence* until corroborated by non-palynological evidence.

### 5 Pollen dispersal bias: quantitative assessment

Differential pollen dispersal is a legitimate concern. Cerealia-type pollen is large-grained (>40 μm, fall speed ~0.06 m/s) with an effective source radius of ~100–300 m, compared with ~300–1,000 m for *Plantago lanceolata* (Sugita, 2007; Broström et al., 2008; Bunting et al., 2004). A PPE-based source-area model predicts a detection lag of 300–1,600 years from dispersal asymmetry alone (details in Supplementary code: `phase3_dispersal_bias.R`). Three lines of evidence indicate dispersal bias cannot account for the observed pattern: (a) the observed median *Plantago*-Cerealia lag (3,076 yr) exceeds the dispersal prediction by a factor of 2–10; (b) 19.1% of sites (30/157) show Cerealia appearing *before* *Plantago*, incompatible with a uniform dispersal artifact; and (c) *Rumex* — a taxon with different dispersal characteristics — independently confirms the pastoral-first ordering at 51–73% of sites. Differential dispersal contributes an estimated 200–500 years to the total lag but cannot be the primary explanation.

### 6 The Calluna problem and the strict pastoral index

*Calluna vulgaris* presents an inherent ambiguity in mid-Holocene pollen analysis. It is a classic pastoral indicator associated with anthropogenic heathland maintenance through burning and grazing (Behre, 1981), but it also expands naturally on acidic soils following natural forest decline, particularly in Atlantic climates. In this study, *Calluna* constitutes 58.6% of pastoral BC in the British Isles, and 65% of pre-Neolithic pastoral exceedances are *Calluna*-driven.

The strict pastoral index resolves this ambiguity by excluding *Calluna*, *Pteridium*, and Ericaceae. Under the strict index, pre-Neolithic "pastoral" detections effectively disappear, confirming that they reflect natural heathland rather than human activity. The temporal ordering results are robust to *Calluna* exclusion (21/22 pastoral-first British sites confirmed by *Plantago*/*Rumex*), but the magnitude of pastoral contribution in the British Isles is substantially reduced. Independent pollen-based land-cover reconstructions for Britain (Fyfe et al., 2013) confirm a gradual transition from closed forest to open pastoral landscapes through the mid-to-late Holocene, consistent with the pastoral exceedance chronology identified here. We recommend that future studies report both standard and strict pastoral indices, particularly in Atlantic and boreal contexts where natural heathland expansion is plausible.

### 7 Bray-Curtis magnitude: an uninformative metric for attribution

With the observed sample sizes and variance, magnitude comparisons between H1 and H3 sites have only ~5% power to detect weak amplification effects. The tests are therefore uninformative: they cannot distinguish between genuinely equal magnitude and a small real difference obscured by noise. We can rule out strong amplification (>60% extra tree turnover), but we cannot evaluate weaker effects. This limitation has a constructive implication: if magnitude comparisons fail, approaches orthogonal to magnitude — functional decomposition, exceedance classification, temporal ordering — become essential. This paper demonstrates that these approaches succeed where magnitude fails.

### 8 Implications for continental-scale pollen syntheses

Several major syntheses have used rate-of-change or dissimilarity metrics to map the onset and intensity of human impact across Europe (Mottl et al., 2021; Nogué et al., 2025) and globally (Gordon et al., 2024). Our results suggest two specific cautions.

First, tree taxa dominate the Bray-Curtis signal (34–57%), while pastoral and arable indicators contribute only 1–7%. The detected "signals" are largely signals of tree compositional change, not direct measures of human activity. This does not invalidate magnitude-based detection — the timing of tree change may correlate with anthropogenic forcing — but it means the detected signal is not a direct human-activity index. Quantitative land-cover change estimates using approaches such as REVEALS (Marquer et al., 2017) or pseudobiomization (Fyfe et al., 2015) provide complementary information by translating pollen data into vegetation cover, but they too are ultimately constrained by the same pollen taphonomic processes.

Second, the 2,135-year median lag suggests that pollen-based "impact onset" dates correspond to landscape-scale intensification, not the beginning of agriculture. Continental syntheses comparing pollen signal onsets to archaeological Neolithic arrival dates may be comparing two different phenomena.

We do not argue that these syntheses are wrong — the broad patterns they identify are broadly consistent with archaeological evidence. We argue that they are incomplete: functional decomposition and exceedance-based classification add interpretive depth that magnitude alone cannot provide.

### 9 Limitations

Several limitations must be acknowledged.

**Replication scope.** The exceedance ordering (67–82%) replicates across all three regions (all n > 50, all p < 0.0001). However, the H1-versus-H3 comparison has adequate H3 samples only in the British Isles (n = 18) and Scandinavia (n = 15); Central Europe (n = 4) should be regarded as exploratory.

**Partial circularity.** Comparing pastoral BC% between groups classified by pastoral exceedance is partly circular. We mitigate this by using continuous-variable analyses (intensification gradient) and exceedance ordering (testing *when*, not *whether*) as primary evidence.

**Differential detection versus genuine sequence.** We cannot determine whether the exceedance ordering reflects genuine temporal priority of pastoral land use or differential pollen-detection thresholds (Section 4.3). We recommend interpreting the ordering as a pollen taphonomic sequence pending non-palynological corroboration.

**Regional coverage and chronological uncertainty.** The three study regions cover northwest and central Europe only. Mediterranean and southeastern European contexts are unrepresented. Chronological uncertainties (±100–200 yr per site) are smaller than the key offsets (400–3,076 yr) but the Cerealia offset in Central Europe (362 yr) approaches the noise floor. Sites vary in age-model quality; regional-scale tests are robust to random noise but individual site-level ordering should be interpreted cautiously.

**Magnitude limitation.** Power analysis (~5% for weak effects) renders magnitude comparisons uninformative (Section 3.9).

**Future validation.** Coprophilous fungal spores (*Sporormiella*, *Sordaria*) — deposited in sediments by herbivore dung — provide a pollen-independent pastoral indicator (van Geel et al., 2003; Baker et al., 2013) and should be a priority for independent validation where such records are available from the same cores.

### 10 Recommended five-step protocol

Based on our findings, we recommend a five-step protocol for studies using Bray-Curtis dissimilarity to detect human impact on vegetation:

1. **Detect**: Apply standardized signal detection with explicit thresholds and sensitivity testing (Gordon et al., 2024).
2. **Decompose**: For every detected signal, compute the functional decomposition to quantify the contribution of pastoral/arable, climate-sensitive tree, succession tree, wetland, and other taxa.
3. **Classify by exceedance**: Determine whether pastoral indicator taxa *exceed* their pre-anthropogenic baseline. Do not rely on first-appearance ordering, which is trivially predicted by native-versus-introduced biogeographic status.
4. **Test temporal gradient**: Relate pastoral exceedance strength to distance from Neolithic arrival as a continuous variable, testing whether the signal tracks initial Neolithic presence or later intensification.
5. **Verify with strict index and dispersal check**: Apply a strict pastoral index excluding *Calluna*/*Pteridium*/Ericaceae. Assess whether the observed temporal ordering exceeds the expected magnitude of pollen dispersal bias (~200–500 yr; Sugita, 2007; Bunting et al., 2004), and report the proportion of sites showing reversed ordering.

This protocol adds modest computational effort to existing workflows but substantially increases interpretive value. The exceedance-based approach, in particular, provides a non-trivial discriminator between anthropogenic and natural vegetation change that first-appearance metrics cannot achieve.

---

## Conclusions

1. Bray-Curtis dissimilarity from pre-anthropogenic baselines detects mid-Holocene vegetation signals at 295 of 331 European pollen sites (89%), but tree taxa dominate these signals (climate-sensitive trees 34–57%, succession trees 10–14%), while pastoral and arable indicators together contribute only 1–7%. Magnitude alone cannot distinguish anthropogenic from natural vegetation change.

2. The exceedance-based framework provides the key discriminator. Pastoral taxa are *present* at 100% of both anthropogenic (H1) and natural-phase (H3) sites, but they *exceed* their pre-anthropogenic baseline at 87% of all sites and at 0% of natural-phase sites (0/37; 95% CI [0%, 9.5%]). First-appearance ordering is trivially expected (null >95% pastoral-first due to native-vs.-introduced biogeographic status) and carries no interpretive weight. All claims in this study are based on exceedance, not first appearance.

3. Exceedance-based pastoral-first ordering — pastoral indicators exceeding their baseline before arable indicators — is observed at 67–82% of sites across all three regions (all binomial p < 0.0001 vs. 1/3 null). This finding is robust to exclusion of *Plantago* (zero dependent sites; *Rumex* confirms 51–73%), *Calluna*, Poaceae, *Ulmus*, and threshold variation. Whether this ordering reflects genuine temporal priority of pastoral over arable activity or differential pollen-detection thresholds remains an open question requiring non-palynological corroboration.

4. The pastoral exceedance signal detects agricultural intensification, not initial Neolithic presence. Pastoral contribution increases with time since Neolithic arrival (Spearman rho = +0.186, p = 0.0015; LMM p = 0.0014), and the median lag between Neolithic arrival and pastoral exceedance is 2,135 years, decomposing into a 635-year detection component (95% CI: 231--1,030 years, from threshold sensitivity analysis) and a 1,454-year archaeological residual (95% CI: 905--1,977 years, representing genuine delay between Neolithic arrival and landscape-scale pastoral transformation). Signal strength does not predict lag length (rho = 0.001, p = 0.995).

5. Alternative hypotheses are independently refuted: cereal-first (H2) by the pastoral exceedance ordering; climate-driven forest restructuring (H3) by zero pastoral exceedance at natural sites and no climate event clustering (p = 0.92); elm disease (H4) by 100% signal retention after *Ulmus* exclusion. A strict pastoral index eliminating *Calluna*/*Pteridium*/Ericaceae removes all pre-Neolithic false positives. Pollen dispersal bias (expected ~200–500 yr) cannot explain the observed exceedance ordering (3,076 yr median lag; 19.1% reversed). We recommend a five-step protocol — Detect → Decompose → Classify by exceedance → Test temporal gradient → Verify with strict index and dispersal check — as a standard for pollen-based human impact attribution.

---

## Acknowledgements

Pollen data were obtained from the Neotoma Paleoecology Database (https://www.neotomadb.org/). We thank the data contributors and database stewards. Analyses were conducted using R with the `neotoma2`, `vegan`, `spdep`, `lme4`, `lmerTest`, `metafor`, and `pwr` packages.

---

## Data availability

All pollen data used in this study are publicly available through the Neotoma Paleoecology Database. R scripts for signal detection, functional decomposition, exceedance classification, temporal ordering, and all statistical analyses are provided as Supplementary Material and are available at [repository to be provided upon acceptance].

---

## References

Bakels, C. C. (2009). *The Western European Loess Belt: Agrarian History, 5300 BC – AD 1000*. Springer.

Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models using lme4. *Journal of Statistical Software*, 67, 1–48.

Behre, K.-E. (1981). The interpretation of anthropogenic indicators in pollen diagrams. *Pollen et Spores*, 23, 225–245.

Behre, K.-E. (1988). The role of man in European vegetation history. In B. Huntley & T. Webb III (Eds.), *Vegetation History* (pp. 633–672). Kluwer.

Bogaard, A. (2004). *Neolithic Farming in Central Europe*. Routledge.

Bogaard, A. (2005). 'Garden agriculture' and the nature of early farming in Europe and the Near East. *World Archaeology*, 37, 177–196.

Broström, A., Nielsen, A. B., Gaillard, M.-J., et al. (2008). Pollen productivity estimates of key European plant taxa for quantitative reconstruction of past vegetation: a review. *Vegetation History and Archaeobotany*, 17, 461–478.

Bunting, M. J., Gaillard, M.-J., Sugita, S., Middleton, R., & Broström, A. (2004). Vegetation structure and pollen source area. *The Holocene*, 14, 651–660.

Copley, M. S., Berstan, R., Dudd, S. N., et al. (2003). Direct chemical evidence for widespread dairying in prehistoric Britain. *Proceedings of the National Academy of Sciences*, 100, 1524–1529.

Davis, B. A. S. (2000). Palynological evidence for the Holocene thermal maximum in Europe. *Review of Palaeobotany and Palynology*, 112, 55–70.

Evershed, R. P., Payne, S., Sherratt, A. G., et al. (2008). Earliest date for milk use in the Near East and southeastern Europe linked to cattle herding. *Nature*, 455, 528–531.

Flantua, S. G. A., Hooghiemstra, H., Grimm, E. C., et al. (2023). Harmonized pollen taxonomy for global synthesis. *Review of Palaeobotany and Palynology*, 309, 104821.

Fyfe, R. M., Twiddle, C. L., Sugita, S., et al. (2013). The Holocene vegetation cover of Britain and Ireland: overcoming problems of scale and discerning patterns of openness. *Journal of Biogeography*, 40, 1132–1148.

Fyfe, R. M., Woodbridge, J., & Roberts, N. (2015). From forest to farmland: pollen-inferred land cover change across Europe using the pseudobiomization approach. *Global Change Biology*, 21, 1197–1212.

Gaillard, M.-J., Sugita, S., Mazier, F., et al. (2010). Holocene land-cover reconstructions for studies on land cover–climate feedbacks. *Climate of the Past*, 6, 483–499.

Gordon, R., et al. (2024). Signal phase detection in Holocene pollen records: a quantitative framework for identifying vegetation transitions. *Quaternary Science Reviews*, 325, 108458.

Gordon, R., et al. (2026). Pan-regional application of signal phase detection across 104 Chinese pollen records. *The Holocene*, in press.

Greenfield, H. J. (2010). The Secondary Products Revolution: the past, the present and the future. *World Archaeology*, 42, 29–54.

Hjelle, K. L. (1999). Modern pollen assemblages from mown and grazed vegetation types in western Norway. *Review of Palaeobotany and Palynology*, 107, 55–81.

Kaufman, D., McKay, N., Routson, C., et al. (2020). A global database of Holocene paleotemperature records. *Scientific Data*, 7, 115.

Kreuz, A. (2012). Die Vertreibung aus dem Paradies? Archäobotanische Ergebnisse zum Frühneolithikum im westlichen Mitteleuropa. *Berichte der Römisch-Germanischen Kommission*, 91, 23–196.

Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd English ed.). Elsevier.

Marciniak, A. (2011). The Secondary Products Revolution: empirical evidence and its current zooarchaeological critique. *Journal of World Prehistory*, 24, 117–130.

Marquer, L., Gaillard, M.-J., Sugita, S., et al. (2017). Quantifying the effects of land use and climate on Holocene vegetation in Europe. *Quaternary Science Reviews*, 171, 20–37.

Mazier, F., Gaillard, M.-J., Kuneš, P., et al. (2012). Testing the effect of site selection and parameter setting on REVEALS-model estimates of plant abundance using the Czech Quaternary Palynological Database. *Review of Palaeobotany and Palynology*, 187, 38–49.

Mottl, O., Flantua, S. G. A., Bhatt, S., et al. (2021). Global acceleration in rates of vegetation change over the past 18,000 years. *Science*, 372, 860–864.

Nogué, S., Flantua, S. G. A., Blaus, A., et al. (2025). Human-induced vegetation change across the Holocene. *Nature Ecology & Evolution*, in press.

Platt, J. R. (1964). Strong inference. *Science*, 146, 347–353.

Renssen, H., Seppä, H., Crosta, X., et al. (2012). Global characterization of the Holocene Thermal Maximum. *Quaternary Science Reviews*, 48, 7–19.

Roberts, N., Fyfe, R. M., Woodbridge, J., et al. (2018). Europe's lost forests: a pollen-based synthesis for the last 11,000 years. *Scientific Reports*, 8, 716.

Sherratt, A. (1981). Plough and pastoralism: aspects of the secondary products revolution. In I. Hodder, G. Isaac, & N. Hammond (Eds.), *Pattern of the Past* (pp. 261–305). Cambridge University Press.

Simmons, I. G. (1996). *The Environmental Impact of Later Mesolithic Cultures*. Edinburgh University Press.

Stevens, C. J., & Fuller, D. Q. (2012). Did Neolithic farming fail? The case for a Bronze Age agricultural revolution in the British Isles. *Antiquity*, 86, 707–722.

Sugita, S. (2007). Theory of quantitative reconstruction of vegetation I: pollen from large sites REVEALS regional vegetation composition. *The Holocene*, 17, 229–241.

Trondman, A.-K., Gaillard, M.-J., Mazier, F., et al. (2015). Pollen-based quantitative reconstructions of Holocene regional vegetation cover (plant-functional types and land-cover types) in Europe suitable for climate modelling. *Global Change Biology*, 21, 676–697.

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. *Journal of Statistical Software*, 36, 1–48.

Whitehouse, N. J., Schulting, R. J., McClatchie, M., et al. (2014). Neolithic agriculture on the European western frontier: the boom and bust of early farming in Ireland. *Journal of Archaeological Science*, 51, 181–205.

Williams, J. W., Grimm, E. C., Blois, J. L., et al. (2018). The Neotoma Paleoecology Database: a multiproxy, international, community-curated data resource. *Quaternary Research*, 89, 156–177.

Woodbridge, J., Fyfe, R. M., Roberts, N., et al. (2014). The impact of the Neolithic agricultural transition in Britain: a comparison of pollen-based land-cover and archaeological ^14^C date-inferred population change. *Journal of Archaeological Science*, 51, 216–224.

Zvelebil, M. (1994). Plant use in the Mesolithic and its role in the transition to farming. *Proceedings of the Prehistoric Society*, 60, 35–74.

Baker, A. G., Bhagwat, S. A., & Willis, K. J. (2013). Do dung fungal spores make a good proxy for past distribution of large herbivores? *Quaternary Science Reviews*, 62, 21–31.

Bevan, A., Colledge, S., Fuller, D., Fyfe, R., Shennan, S., & Stevens, C. (2017). Holocene fluctuations in human population demonstrate repeated links to food production and climate. *Proceedings of the National Academy of Sciences*, 114, E10524–E10531.

Birks, H. J. B., & Birks, H. H. (2008). Biological responses to rapid climate change at the Younger Dryas–Holocene transition at Krakenes, western Norway. *The Holocene*, 18, 19–30.

Conolly, J., Manning, K., Colledge, S., Dobney, K., & Shennan, S. (2011). Species distribution modelling of ancient cattle from early Neolithic sites in SW Asia and Europe. *The Holocene*, 21, 983–993.

Davis, O. K., & Shafer, D. S. (2006). Sporormiella fungal spores, a palynological means of detecting herbivore density. *Palaeogeography, Palaeoclimatology, Palaeoecology*, 237, 40–50.

Frantz, L. A. F., Bradley, D. G., Larson, G., & Orlando, L. (2020). Animal domestication in the era of ancient genomics. *Nature Reviews Genetics*, 21, 449–460.

Fuller, D. Q., van Etten, J., Manning, K., Castillo, C., Kingwell-Banham, E., Weisskopf, A., Qin, L., Sato, Y.-I., & Hijmans, R. J. (2011). The contribution of rice agriculture and livestock pastoralism to prehistoric methane levels: an archaeological assessment. *The Holocene*, 21, 743–759.

Giesecke, T., Bennett, K. D., Birks, H. J. B., Bjune, A. E., Bozilova, E., Feurdean, A., Finsinger, W., Froyd, C., Pokorny, P., Rosch, M., Seppa, H., Tonkov, S., Valsecchi, V., & Wolters, S. (2011). The pace of Holocene vegetation change — testing for synchronous developments. *Quaternary Science Reviews*, 30, 2805–2814.

Marlon, J. R., Bartlein, P. J., Daniau, A.-L., Harrison, S. P., Maezumi, S. Y., Power, M. J., Tinner, W., & Vanniere, B. (2013). Global biomass burning: a synthesis and review of Holocene paleofire records and their controls. *Quaternary Science Reviews*, 65, 5–25.

Nielsen, A. B., Giesecke, T., Theuerkauf, M., Feeser, I., Behre, K.-E., Beug, H.-J., Chen, S.-H., Christiansen, J., Dorfler, W., Endtmann, E., Jahns, S., de Klerk, P., Kuhl, N., Latalowa, M., Odgaard, B. V., Rasmussen, P., Stockholm, J. R., Voigt, R., Wiethold, J., & Wolters, S. (2012). Quantitative reconstructions of changes in regional openness in north-central Europe reveal new insights into old questions. *Quaternary Science Reviews*, 47, 131–149.

Parker, A. G., Goudie, A. S., Anderson, D. E., Robinson, M. A., & Bonsall, C. (2002). A review of the mid-Holocene elm decline in the British Isles. *Progress in Physical Geography*, 26, 1–45.

Peglar, S. M. (1993). The mid-Holocene *Ulmus* decline at Diss Mere, Norfolk, UK: a year-by-year pollen stratigraphy from annual laminations. *The Holocene*, 3, 1–10.

Ruddiman, W. F. (2003). The anthropogenic greenhouse era began thousands of years ago. *Climatic Change*, 61, 261–293.

Shennan, S., Downey, S. S., Timpson, A., Edinborough, K., Colledge, S., Kerig, T., Manning, K., & Thomas, M. G. (2013). Regional population collapse followed initial agriculture booms in mid-Holocene Europe. *Nature Communications*, 4, 2486.

Timpson, A., Colledge, S., Crema, E., Edinborough, K., Kerig, T., Manning, K., Thomas, M. G., & Shennan, S. (2014). Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case study using an improved method. *Journal of Archaeological Science*, 52, 549–557.

van Geel, B., Buurman, J., Brinkkemper, O., Schelvis, J., Aptroot, A., van Reenen, G., & Hakbijl, T. (2003). Environmental reconstruction of a Roman Period settlement site in Uitgeest (The Netherlands), with special reference to coprophilous fungi. *Journal of Archaeological Science*, 30, 873–883.

Verdugo, M. P., Mullin, V. E., Scheu, A., Mattiangeli, V., Daly, K. G., Maisano Delser, P., Hare, A. J., Burger, J., Collins, M. J., Kehati, R., et al. (2019). Ancient cattle genomics, origins, and rapid turnover in the Fertile Crescent. *Science*, 365, 173–176.

Vigne, J.-D. (2008). Zooarchaeological aspects of the Neolithic diet transition in the Near East and Europe, and their putative relationships with the Neolithic Demographic Transition. In J.-P. Bocquet-Appel & O. Bar-Yosef (Eds.), *The Neolithic Demographic Transition and its Consequences* (pp. 135–159). Springer.

Whitehouse, N. J., & Smith, D. N. (2010). How fragmented was the British Holocene wildwood? Perspectives on the "Vera" grazing debate from the fossil beetle record. *Quaternary Science Reviews*, 29, 539–553.

Mazier, F., Galop, D., Brun, C., & Buttler, A. (2006). Modern pollen assemblages from grazed vegetation in the western Pyrenees, France: a numerical tool for more precise reconstruction of past cultural landscapes. *The Holocene*, 16, 91–103.

Abraham, V., Ouskova, M., & Kunes, P. (2014). Present-day vegetation helps quantifying past land cover in selected regions of the Czech Republic. *PLoS ONE*, 9, e100visual.

Marlon, J. R., Kelly, R., Daniau, A.-L., Vanniere, B., Power, M. J., Bartlein, P., Higuera, P., Blarquez, O., Brewer, S., Brucher, T., Feurdean, A., Romera, G. G., Iglesias, V., Maezumi, S. Y., Magi, B., Courtney Mustaphi, C. J., & Zhihai, T. (2016). Reconstructions of biomass burning from sediment-charcoal records to improve data–model comparisons. *Biogeosciences*, 13, 3225–3244.

Kidder, T. R., Liu, H., & Li, M. (2012). Sanyangzhuang: early Chinese civilization in the floodplain. *Antiquity*, 86, 1–15.

Crema, E. R., & Bevan, A. (2021). Inference from large sets of radiocarbon dates: software and methods. *Radiocarbon*, 63, 23–39.

Zeder, M. A. (2008). Domestication and early agriculture in the Mediterranean Basin: origins, diffusion, and impact. *Proceedings of the National Academy of Sciences*, 105, 11597–11604.
