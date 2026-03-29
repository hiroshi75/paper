# Three Diagnostic Principles for Evaluating Paleoecological Proxy Performance

**Version**: 2.0
**Date**: 2026-03-29
**Type**: Perspective / Theory
**Target journal**: Trends in Ecology & Evolution / Quaternary Science Reviews

---

## Abstract

Paleoecological proxies are measurement instruments, yet the discipline lacks an explicit diagnostic framework specifying what proxies detect, what they miss, and how to tell the difference. Across eleven studies analysing more than 1,500 pollen, charcoal, and radiocarbon records from five continents, a recurring finding has emerged: detection outcomes depend as much on the instrument as on the signal. Here we synthesise these findings into a diagnostic framework organised around three principles and three tools. The principles are: (1) detection requires domain overlap between the indicator's measurement space and the transformation's effect space (identifiability); (2) the detection threshold for anthropogenic impact is exceedance above baseline variability, not mere presence (exceedance); and (3) different proxies detect different domains, and their union always exceeds any single proxy (complementarity). The tools are: first-difference detrending as a signal/noise separator, spectral decomposition as a mechanism identifier, and the indicator dependency gap as a measurement quality metric. Together, these elements constitute a practical diagnostic framework for evaluating paleoecological proxy performance --- one that makes the instrument's limitations an explicit, quantifiable part of every inference.

**Keywords**: diagnostic framework, paleoecological proxies, detection function, identifiability, exceedance, multi-proxy complementarity, indicator dependency, detrending, spectral analysis

---

## 1. Introduction: Why Paleoecology Needs Measurement Theory

Every empirical science eventually confronts the distinction between the thing measured and the act of measurement. Physics formalised this through quantum diagnostic framework; clinical medicine through sensitivity, specificity, and receiver operating characteristics; remote sensing through the modulation transfer function. Paleoecology has not.

The consequences are not abstract. Continental-scale syntheses now routinely use pollen-derived metrics to determine when and how intensely human societies transformed terrestrial vegetation (Mottl et al., 2021; Nogué et al., 2025). Charcoal records reconstruct fire regimes. Summed probability distributions (SPDs) of radiocarbon dates serve as demographic proxies. These measurements inform debates about Anthropocene onset, conservation baselines, and the deep-time relationship between human societies and biodiversity (Ellis et al., 2021; Stephens et al., 2019). Yet in each case, investigators report what the proxy shows without formally specifying what the proxy can and cannot detect --- and the gap between these two is often the gap between correct and incorrect inference.

Consider a concrete example. European pastoral indicator taxa (*Plantago lanceolata*, *Rumex*, Cerealia-type) detect agricultural impact at 87% of European pollen sites. Applied to the same analytical framework in eastern North America --- an independent centre of plant domestication with a well-documented agricultural tradition (Smith, 2006) --- the detection rate drops to 0% (Gordon et al., 2026a). The correct interpretation is not that eastern North America lacked agriculture but that the measurement instrument was calibrated against the wrong signal. This is a measurement failure, not an ecological finding, and it can only be recognised as such within a framework that treats the proxy as an instrument with specifiable detection properties.

This paper synthesises results from eleven studies conducted within the ArchEco research programme into a unified diagnostic framework for paleoecological proxies. The theory is built around three elements: a detection function that formalises what proxies measure (Section 2), three principles that govern when detection succeeds and fails (Section 3), and three diagnostic tools for evaluating proxy reliability (Section 4). Our aim is not mathematical elegance but practical utility: to give paleoecologists a language for making their instruments' limitations as explicit as their results.

---

## 2. The Detection Function *D*(*I*, *T*, *s*)

### 2.1 Definition

We define the detection function as the probability that an indicator system *I* detects a transformation *T* at signal strength *s*:

> *D*(*I*, *T*, *s*) = Pr(indicator system *I* classifies transformation *T* as present | signal strength = *s*)

where:

- **I** is the complete indicator specification: which proxy (pollen, charcoal, radiocarbon), which taxa or variables, and which detection criterion (threshold, statistical test, dissimilarity metric).
- **T** is the transformation to be detected: a specific change in the ecological or demographic system (e.g., forest clearance, crop cultivation, fire regime change, population decline).
- **s** is the signal strength: the magnitude of *T* in its native domain (e.g., percentage of landscape cleared, area under cultivation, fire frequency change).

This notation makes three things explicit that are usually implicit. First, detection is a probability, not a certainty. Second, it depends jointly on the indicator system and the transformation --- neither alone determines the outcome. Third, detection probability varies with signal strength, but the *form* of that variation depends on the relationship between *I* and *T*.

### 2.2 Three structural constraints

Across the eleven studies, three structural constraints on the detection function recur:

**Orthogonality.** When the ecological domain perturbed by *T* has no overlap with the domain monitored by *I*, the detection probability is invariant to signal strength: dD/ds approx 0. Increasing the magnitude of a transformation that operates in a dimension the indicator does not observe adds signal in an unmonitored axis. This is a geometric property of the measurement, not an empirical accident. European pastoral indicators monitoring forest-clearance taxa cannot detect Eastern Agricultural Complex cultivation that operates through Amaranthaceae increases, regardless of how extensive that cultivation becomes (Gordon et al., 2026a). The indicator dependency gap on every non-European continent exceeds 60 percentage points (Gordon et al., 2026b), confirming that orthogonality is not a special case but a pervasive structural feature.

**Masking.** When the domain of *T* overlaps with *I* but is confounded by a larger signal from a different process, the transformation may be detectable in principle but invisible in practice. Tree taxa show exceedance at 87--94% of temperate-forest pollen sites on every continent tested, driven by climate-induced species migrations and natural successional dynamics (Gordon et al., 2026a, 2026b). This near-universal exceedance renders tree taxa non-diagnostic for any specific agricultural system: the signal of anthropogenic compositional change is masked by the signal of natural Holocene forest dynamics.

**Redundancy deficit.** When a transformation produces only a single measurable change in the indicator's domain, the detection function depends critically on the taxonomic resolution and taphonomic fidelity of that single channel. EAC cultivation produces a signal primarily through family-level Amaranthaceae increase; because Amaranthaceae includes both cultivated (*Chenopodium*) and wild species, the detection rate (39%) is lower than for European pastoral indicators (87%), which benefit from multiple independently diagnostic taxa (*Plantago*, *Rumex*, Cerealia-type, Poaceae) that provide redundant confirmation (Gordon et al., 2026a). Forward simulation confirms this: detection sensitivity at 80% requires 2--3 times greater signal strength for single-taxon indicators than for multi-taxon syndromes.

### 2.3 The detection function as a unifying lens

We emphasise that the detection function as presented here is a **conceptual framework, not a parametric model**. Unlike Item Response Theory in psychometrics, where detection functions have estimable parametric forms, D(I,T,s) in its current form is a notation for organising empirical observations rather than a theory that generates novel predictions. The forward simulation results (Gordon et al., 2026a) suggest that logistic or step-function forms may be appropriate for specific proxy-transformation pairs, but a fully estimable general form has not yet been developed. We regard formalising D into an estimable model as the most important future direction for this framework.

Despite this limitation, the detection function notation is not merely decorative. It unifies otherwise disparate observations under a common logic. The presence/exceedance distinction (Section 3.2), the indicator dependency gap (Section 4.3), multi-proxy complementarity (Section 3.3), and even the detrending diagnostic (Section 4.1) are all instances of asking: what is *D*(*I*, *T*, *s*) for this specific combination of instrument and signal? When *D* is high, we have a reliable measurement. When *D* is low or zero, we have a systematic blind spot. When *D* varies across indicator choices more than across ecological conditions, we have an instrument problem masquerading as an ecological finding.

---

## 3. Three Principles of Proxy Detection

### 3.1 Principle 1: Detection requires domain overlap (the identifiability law)

**Statement.** A proxy framework detects only those transformations that alter its measurement domain at its resolution. Formally, if indicator set *I* spans a subspace of the compositional or measurement space, and transformation *T* perturbs a different subspace with no overlap, then *T* is operationally unidentifiable under *I*: dD/ds approx 0, regardless of *s*.

**Empirical basis.** The identifiability principle was first demonstrated for within-proxy variation: different pollen indicator sets applied to the same 442 sites produced detection rates ranging from 0% to 94% (Gordon et al., 2026a). The three-domain decomposition --- structural (Type A: forest clearance), crop-specific (Type B: cultivated taxa), and compositional (Type C: tree reorganisation) --- showed that each domain is independently measurable and that indicators designed for one domain are structurally blind to others.

Global extension across 776 sites on five continents confirmed that the indicator dependency gap exceeds 60 percentage points on every non-European continent (Gordon et al., 2026b). European pastoral indicators detect 87% of European sites but only 0--41% elsewhere. Region-specific crop indicators recover additional signal at rates of 13--58%. Tree taxa show near-universal exceedance (87--94%) everywhere. These are not marginal differences in sensitivity; they are categorical differences in what the instrument can see.

The forward detection simulation (Gordon et al., 2026a) provided the formal demonstration. When synthetic pollen assemblages are perturbed by a transformation that operates through taxa outside the indicator set, detection probability remains at baseline (approx 5% false-positive rate) no matter how large the perturbation. The detection function is flat: dD/ds = 0. Conversely, when the perturbation operates through taxa within the indicator set, detection rises monotonically with signal strength, reaching 80% at moderate perturbation levels.

**Practical implication.** Any cross-cultural comparison of agricultural impact that uses a single indicator set without testing for domain overlap is vulnerable to systematic false negatives in regions where the agricultural tradition operates in a different domain. The identifiability principle requires that investigators specify which transformation domain their indicators monitor and acknowledge which domains they do not.

### 3.2 Principle 2: Detection threshold is exceedance, not presence (the exceedance law)

**Statement.** For native taxa, mere presence in the paleoecological record is trivially expected under any ecological model and carries no diagnostic weight for anthropogenic impact. Only exceedance --- anomalous increase above baseline variability (operationally: baseline mean + 2 SD) --- constitutes evidence that an external forcing agent has shifted the taxon's abundance beyond the range of natural Holocene variation.

**Empirical basis.** Pastoral indicator taxa (*Plantago lanceolata*, *Rumex*, Poaceae) are native European species present in pollen records throughout the Holocene. Their first appearance before Cerealia-type pollen is predicted by a simple generalist-specialist null model at >95% of sites (Gordon et al., 2026c). First-appearance ordering thus has no diagnostic power: it reflects native status, not human activity.

Exceedance, by contrast, discriminates sharply. At 246 of 283 evaluable European sites (87%), pastoral taxa exceed their pre-anthropogenic baseline during the agricultural period. At all 37 natural-phase control sites, pastoral taxa are present but never exceed baseline --- a 0% false-positive rate (Gordon et al., 2026c). The exceedance threshold is robust to sensitivity analysis: from 1.5 SD to 3.0 SD, the proportion of H1 sites varies (49--75%) but the discrimination between anthropogenic and natural signals is maintained.

The exceedance principle applies beyond pastoral indicators. In eastern North America, Amaranthaceae are present throughout the Holocene but show zero exceedance before 3,000 cal BP (0% false-positive rate across 93 sites over >2,000 years of pre-agricultural baseline). The first exceedances (3,000--2,000 BP) cluster in the archaeologically documented EAC heartland (Gordon et al., 2026a).

**Formal interpretation.** The distinction between presence and exceedance maps onto the detection function as a threshold property. Below the exceedance threshold, *D* approx 0 even when the taxon is present, because presence is within the range of natural variability. Above the threshold, *D* rises sharply. The threshold defines the boundary between the noise floor (natural Holocene variability) and the signal domain (anthropogenic forcing). The critical insight is that this threshold is not arbitrary: it is the empirically determined boundary of the system's natural state space.

### 3.3 Principle 3: Multi-proxy union exceeds any single proxy (the complementarity law)

**Statement.** Different proxies detect different ecological domains. The union of multiple proxies always detects a larger fraction of total human impact than any single proxy, and the discordant fraction (signal in one proxy but not the other) represents genuine domain-specific impacts, not noise.

**Empirical basis.** At 111 sites where both pollen and charcoal records exist in the same sediment core, the 2 x 2 contingency table reveals: 21% dual exceedance, 20% charcoal-only exceedance, 14% pollen-only exceedance, 46% neither (Gordon et al., 2026d). Fisher's exact test confirms non-independence (p = 0.002, OR = 3.51), but the 33% discordance rate demonstrates that each proxy has a distinct detection domain. The 22 charcoal-only sites represent fire regime changes invisible to pollen --- potentially anthropogenic burning that did not produce lasting deforestation. The 15 pollen-only sites represent vegetation transformations without detectable fire intensification.

This result extends the identifiability framework from intra-proxy to inter-proxy comparison. Within pollen, different indicator sets detect different agricultural domains (Types A, B, C). Between proxies, pollen and charcoal detect different impact mechanisms altogether. The addition of a Type D domain (fire regime, detected by charcoal) to the existing A/B/C framework means that complete impact assessment requires at minimum a two-proxy approach, and the degree of complementarity is empirically measurable through the discordance rate.

**Formal interpretation.** Let *I*_1 and *I*_2 be two indicator systems monitoring domains *M*_1 and *M*_2. If *M*_1 and *M*_2 are not identical --- if each proxy responds to at least some transformations the other does not --- then the detection function of the union is strictly greater than either component:

> *D*(*I*_1 union *I*_2, *T*, *s*) >= max(*D*(*I*_1, *T*, *s*), *D*(*I*_2, *T*, *s*))

with equality only when one domain is a subset of the other. In practice, pollen and charcoal have substantially non-overlapping domains, making the union considerably more powerful. The 33% discordance rate at same-core sites is a direct measure of the complementarity gain.

---

## 4. Diagnostic Tools

The three principles describe the structure of the detection function. But investigators need practical tools to evaluate whether their specific analysis falls into a detection regime where inference is reliable or one where it is compromised. Three diagnostic tools emerge from the research programme, each addressing a different aspect of measurement quality.

### 4.1 Detrending as a signal/noise separator

**The problem.** Holocene paleoecological and archaeological time series are strongly trended: post-glacial vegetation succession, progressive human landscape modification, and taphonomic loss of older material all drive broadly parallel trends. Correlating two trended time series produces spurious associations that can reach |rho| > 0.8 with p < 0.0001 (Gordon et al., 2026e).

**The diagnostic.** First-difference detrending --- correlating changes between successive time steps rather than levels --- removes shared monotonic trends and tests whether two variables genuinely covary on a step-by-step basis. Across five region-proxy combinations (three regions for SPD x pollen richness, two regions for SPD x charcoal), all raw correlations were strong (|rho| = 0.26--0.85), but all first-difference correlations collapsed to near-zero (|rho| = 0.003--0.086), with only one marginal exception (Scandinavia: rho = 0.454, p = 0.052). GAM-residual analyses confirmed the pattern: SPD contributed 0% additional deviance explained beyond the smooth effect of time in the British Isles (Gordon et al., 2026e).

The same detrending applied to a different question --- demographic synchrony across seven regions --- produced the opposite result. After first-differencing, 20 of 21 inter-regional SPD pairs showed significant positive correlation (mean detrended r = 0.458), surviving six independent robustness tests (Gordon et al., 2026f). The contrast is diagnostic: SPD x pollen correlations vanish after detrending (spurious); SPD x SPD correlations survive (genuine). Detrending does not indiscriminately destroy signal; it selectively eliminates associations that were never real.

**Formal role.** In detection function terms, detrending asks: does the relationship between *I* and *T* persist when shared confounds are removed? If *D* drops to baseline after detrending, the original detection was an artifact of shared trends, not a genuine measurement. If *D* survives, the measurement captures real covariation. First-differencing is thus the simplest available test of whether a proxy-based inference reflects the target process or merely the shared background dynamics of Holocene change.

### 4.2 Frequency-domain decomposition as a mechanism identifier

**The problem.** Even when a correlation survives detrending, the mechanism remains ambiguous. Climate forcing, endogenous population dynamics, and cultural diffusion all operate at different characteristic timescales. A correlation that is statistically genuine may still be attributed to the wrong cause if the timescale of the signal is not examined.

**The diagnostic.** Spectral (bandpass) decomposition separates a time series into frequency bands, allowing the investigator to ask: at which timescale is the signal concentrated? For the demographic synchrony result, bandpass filtering revealed a striking frequency structure (Gordon et al., 2026f):

- **High frequency (100--500 years):** Mean r = 0.420, all 21 pairs positive, 20 significant. This is the dominant synchrony signal.
- **Medium frequency (500--1,500 years):** Mean r = 0.057, 11/21 positive, 11/21 significant. Near zero.
- **Low frequency (1,500--5,000 years):** Mean r = 0.095, pair-specific, ranging from r = -0.96 to +0.96.

The medium-frequency band (500--1,500 years) is precisely where climate forcing --- Bond cycles, century-scale droughts --- should produce synchrony if climate were the primary driver. Its near-absence constitutes a frequency-domain falsification of simple climate models. The signal is concentrated at high frequencies, suggesting endogenous demographic processes or rapid cultural transmission rather than slow climate forcing.

**Formal role.** Frequency decomposition constrains the mechanism by asking: at which timescale does *D*(*I*, *T*, *s*) produce non-zero signal? If the detection function is non-zero only at frequencies inconsistent with the hypothesised forcing, the hypothesis is falsified regardless of the aggregate correlation strength. The frequency that carries a signal constrains the process that generates it.

### 4.3 The indicator dependency gap as a measurement quality metric

**The problem.** When multiple indicator sets are available for the same sites, how should we quantify the degree to which detection outcomes depend on the instrument rather than the phenomenon?

**The diagnostic.** The indicator dependency gap is defined as the difference between the maximum and minimum detection rates across indicator sets within a single region. It directly measures how much the conclusion "this region shows agricultural impact" depends on the choice of indicator framework.

Across five continents (Gordon et al., 2026a, 2026b):

| Continent | Min detection (%) | Max detection (%) | Gap (pp) |
|-----------|------------------:|------------------:|---------:|
| Europe | 87 (pastoral) | 93 (trees) | 6 |
| Eastern N. America | 0 (pastoral) | 94 (trees) | 94 |
| China | 16 (pastoral) | 87 (trees) | 71 |
| Latin America | 28 (pastoral) | 89 (trees) | 61 |
| Sub-Saharan Africa | 28 (pastoral) | 91 (trees) | 63 |

Europe's narrow gap (6 pp) reflects the fact that European indicators were designed for European agricultural systems --- the instrument is calibrated for the signal. Every non-European continent shows a gap exceeding 60 percentage points, meaning that the conclusion about agricultural impact is more dependent on indicator choice than on the ecology of the sites.

**Formal role.** The indicator dependency gap is a direct empirical estimate of how much the detection function *D* varies with *I* for a given set of sites and plausible transformations. A large gap signals that the measurement is dominated by instrument properties rather than signal properties --- the paleoecological equivalent of an instrument with poor calibration. A small gap signals that different measurement approaches converge on the same conclusion, which is the minimum requirement for treating a detection outcome as a property of the system rather than of the instrument.

We propose that the indicator dependency gap be reported as a standard quality metric in all cross-cultural comparisons, analogous to reporting inter-rater reliability in observational studies or inter-method agreement in analytical chemistry.

---

## 5. Implications for Practice

### 5.1 A protocol for proxy-based inference

The three principles and three diagnostics converge on a practical protocol for any proxy-based paleoecological inference:

1. **Specify the detection function.** Before reporting results, state explicitly which indicator system *I* is being used, which transformation domain it monitors, and which domains it cannot detect. This is the identifiability audit.

2. **Apply exceedance, not presence.** For any indicator taxon native to the study region, base detection on exceedance above baseline variability, not on first appearance or mere presence. Report the exceedance threshold and its sensitivity.

3. **Test for domain overlap.** If the study involves cross-cultural comparison, verify that the indicator set has been validated against the agricultural tradition of each region. Report the indicator dependency gap.

4. **Detrend before correlating.** Any correlation between a proxy and another time series must survive first-difference or residual-based detrending. Report both raw and detrended statistics.

5. **Decompose by frequency.** When a detrended correlation is significant, examine its frequency structure to constrain the mechanism. Report which frequency band carries the signal.

6. **Combine proxies.** Where possible, use multiple proxies with non-overlapping detection domains and report the concordance/discordance matrix. The discordant fraction is as informative as the concordant fraction.

### 5.2 What this framework does not do

The diagnostic framework proposed here is deliberately minimal. It does not provide a general statistical model for paleoecological inference; existing approaches (Bayesian age-depth modelling, REVEALS land-cover reconstruction, rcarbon SPD analysis) serve that function. It does not resolve specific archaeological or ecological questions. What it does is provide a meta-level language for evaluating whether the instruments being used to address those questions are capable of doing so. It is a theory about the relationship between the instrument and the phenomenon, not about the phenomenon itself.

The framework also has known limitations. The detection function as presented is qualitative --- we have not specified a parametric form for *D*(*I*, *T*, *s*), though the forward simulation results (Gordon et al., 2026a) suggest that logistic or step-function forms are reasonable approximations. The three-domain classification (A/B/C/D) is likely incomplete; as new proxy types and new agricultural systems are examined, additional domains will emerge. The exceedance threshold (mean + 2 SD) is a practical default, not a universal optimum; different ecological contexts may warrant different thresholds.

### 5.3 The deeper point

The recurring lesson across all eleven studies is that paleoecological proxies are not transparent windows onto the past. They are instruments with specifiable detection properties, systematic blind spots, and measurable limitations. A pollen record does not tell us what the past landscape looked like; it tells us what happened within the specific compositional subspace that pollen can detect, at the taxonomic resolution pollen can achieve, above the noise floor of natural Holocene variability. A charcoal record tells us what happened in the fire domain. A radiocarbon SPD tells us something about the intensity of dateable human activity, filtered through taphonomic loss and calibration curve effects.

This is not a counsel of despair. Physics did not abandon measurement when it recognised the observer effect; it developed diagnostic framework. Medicine did not abandon diagnostic tests when it discovered that sensitivity and specificity are properties of the test, not the disease; it developed test evaluation frameworks. Paleoecology can do the same. The path forward is not to pretend that proxies are unbiased but to make their biases explicit, quantifiable, and correctable.

---

## 6. Conclusions

We have proposed a minimal diagnostic framework for paleoecological proxies, built on three elements: the detection function *D*(*I*, *T*, *s*), three principles of proxy detection (identifiability, exceedance, complementarity), and three diagnostic tools (detrending, frequency decomposition, indicator dependency gap). The framework synthesises results from eleven empirical studies spanning 1,500+ sites across five continents, but its claims are structural rather than site-specific: they follow from the geometry of measurement in high-dimensional compositional space.

The three principles can be summarised as a single sentence: **a proxy detects only those transformations that overlap with its measurement domain, exceed its noise floor, and no single proxy spans the full domain of human impact.** The three diagnostics provide the operational tests: detrending separates signal from shared trend, frequency decomposition constrains the mechanism, and the indicator dependency gap quantifies how much the conclusion depends on the instrument.

The practical consequence is a shift in how paleoecological results should be reported. Instead of "this region shows (or does not show) agricultural impact," the measurement-theoretic formulation requires: "indicator system *I* detects (or does not detect) transformation domain *T* at this signal level, with this indicator dependency gap, and the detection survives (or does not survive) detrending at these frequencies." This is more cumbersome. It is also more honest --- and honesty about what our instruments can and cannot see is the first step toward seeing more clearly.

---

## References

*(Selected references --- full bibliography in preparation)*

- Behre, K.E. (1981). The interpretation of anthropogenic indicators in pollen diagrams. *Pollen et Spores*, 23, 225--245.
- Crema, E.R., & Bevan, A. (2021). Inference from large sets of radiocarbon dates: software and methods. *Radiocarbon*, 63, 23--39.
- Ellis, E.C., et al. (2021). People have shaped most of terrestrial nature for at least 12,000 years. *PNAS*, 118, e2023483118.
- Fritz, G.J. (2019). *Feeding Cahokia: Early Agriculture in the North American Heartland*. University of Alabama Press.
- Gordon, A.Y., et al. (2026a). Indicator-dependent detection of agricultural transformation. *[Paper 8, this series]*.
- Gordon, A.Y., et al. (2026b). Global indicator dependency in pollen-based agriculture detection. *[Paper 9, this series]*.
- Gordon, A.Y., et al. (2026c). Pastoral indicator exceedance distinguishes anthropogenic from natural pollen signals. *[Paper 6, this series]*.
- Gordon, A.Y., et al. (2026d). Different proxies, different domains: pollen and charcoal exceedance. *[Paper 10, this series]*.
- Gordon, A.Y., et al. (2026e). Shared Holocene trends inflate correlations between archaeological and paleoecological variables. *[Paper 5, this series]*.
- Gordon, A.Y., et al. (2026f). Pervasive demographic synchrony across Holocene Eurasia. *[Paper 11, this series]*.
- Granger, C.W.J., & Newbold, P. (1974). Spurious regressions in econometrics. *Journal of Econometrics*, 2, 111--120.
- Mottl, O., et al. (2021). Global acceleration in rates of vegetation change over the past 18,000 years. *Science*, 372, 860--864.
- Mueller, N.G. (2017). *Documenting domestication in a lost crop*. In *Enduring Records*, ed. M.T. Scarry, pp. 45--62.
- Nogué, S., et al. (2025). Global patterns of human impact on past vegetation. *Nature Ecology & Evolution*.
- Smith, B.D. (2006). Eastern North America as an independent center of plant domestication. *PNAS*, 103, 12223--12228.
- Stephens, L., et al. (2019). Archaeological assessment reveals Earth's early transformation through land use. *Science*, 365, 897--902.
- Sugita, S. (2007). Theory of quantitative reconstruction of vegetation I. *The Holocene*, 17, 229--241.

---

**Acknowledgements.** This synthesis draws on data from the Neotoma Paleoecology Database, the p3k14c radiocarbon database, and the Global Charcoal Database. We thank the data contributors and database maintainers whose work made this research programme possible.

**Word count.** ~3,800 words (main text, excluding abstract, references, and tables).
