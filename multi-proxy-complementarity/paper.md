# Different Proxies, Different Domains: Pollen and Charcoal Exceedance Reveal Complementary Human Impacts at 111 Same-Core Sites

**Version**: 1.0
**Date**: 2026-03-28
**Target journal**: Journal of Quaternary Science

---

## Abstract

Multi-proxy paleoecological studies routinely combine pollen and charcoal records from the same sediment core, yet the degree to which these proxies detect different dimensions of human impact remains unquantified at a global scale. Here we exploit the Neotoma Paleoecology Database to identify 216 sites where both pollen and charcoal records exist in the same core, of which 111 yield sufficient temporal depth for dual exceedance analysis. We apply a uniform exceedance framework---comparing post-baseline (< 5,000 cal BP) variability against pre-baseline (> 5,000 cal BP) variability with a > 10% threshold---to classify each proxy independently as showing exceedance or no signal. The resulting 2 x 2 contingency table reveals significant but incomplete concordance: 23 sites (20.7%) show dual exceedance (fire and vegetation change), 22 (19.8%) show charcoal-only exceedance, 15 (13.5%) show pollen-only exceedance, and 51 (45.9%) show neither. Fisher's exact test confirms non-independence (p = 0.002, OR = 3.51), but the 33.3% discordance rate demonstrates that each proxy has a distinct detection domain. The 22 charcoal-only sites represent fire regime changes invisible to pollen analysis---potentially anthropogenic burning that did not produce lasting deforestation. The 15 pollen-only sites represent vegetation transformations without detectable fire intensification. These results formalize the concept of proxy-specific detection domains: pollen monitors the vegetation domain while charcoal monitors the fire domain, and human impacts that operate through only one mechanism are systematically missed by single-proxy studies. We propose extending the identifiability framework of indicator-dependent detection (where different pollen indicators detect different agricultural types) to inter-proxy detection (where different proxies detect different impact mechanisms), adding a Type D fire-regime domain to the existing three-domain classification.

**Keywords**: multi-proxy paleoecology, pollen analysis, charcoal analysis, detection domain, exceedance threshold, Neotoma, human impact, fire regime, identifiability

---

## 1. Introduction

### 1.1 The single-proxy problem in human impact detection

Paleoecological reconstructions of human impact on terrestrial ecosystems rely overwhelmingly on single proxies. Pollen-based studies detect changes in vegetation composition and structure; charcoal-based studies detect changes in fire regime. Both are widely used to infer the timing, extent, and intensity of human transformation of landscapes (Marlon et al., 2008; Roberts et al., 2018; Mottl et al., 2021). Yet the assumption underlying most single-proxy studies---that the chosen proxy captures the full spectrum of human impact---is rarely tested empirically.

The problem is not merely one of resolution or sensitivity. It is structural. Different proxies monitor different ecological domains. Pollen records the taxonomic composition and relative abundance of plant communities. Charcoal records the occurrence and intensity of biomass burning. Human impacts that operate through vegetation clearance (deforestation for agriculture, pastoralism) produce strong pollen signals but may produce weak or absent charcoal signals if clearance proceeds through non-fire mechanisms such as ring-barking, selective felling, or grazing pressure. Conversely, human impacts that operate through fire management---deliberate burning to maintain grasslands, drive game, or clear undergrowth---may produce strong charcoal signals while leaving the regional pollen signature largely unchanged if burning does not lead to permanent forest loss (Bowman et al., 2011; Whitlock & Larsen, 2001).

This asymmetry creates an identifiability problem. A site that shows no pollen exceedance may nonetheless have experienced substantial anthropogenic fire management. A site that shows no charcoal exceedance may nonetheless have experienced deforestation through non-fire mechanisms. In both cases, a single-proxy study would classify the site as showing "no detectable human impact," when in fact impact occurred through a mechanism outside the proxy's detection domain.

### 1.2 From indicator dependency to proxy dependency

In a companion study (Paper 8; this series), we demonstrated that agricultural impact detection within pollen records is structurally dependent on indicator selection. Using 442 pollen records from Europe and eastern North America, we showed that detection outcomes vary from 0% to 93.8% at the same sites depending solely on which pollen taxa are included in the indicator framework. We formalized this as a three-domain classification: Type A (structural transformation: forest clearance), Type B (crop-specific transformation: cultivated taxa), and Type C (compositional transformation: tree species reorganisation). European pastoral indicators detect Type A transformations but miss Type B (indigenous crop cultivation) and Type C (fire-mediated compositional change) entirely.

The present study extends this identifiability framework from intra-proxy to inter-proxy comparison. If different indicator sets within a single proxy (pollen) detect different agricultural transformations, then different proxies (pollen vs. charcoal) should detect different impact mechanisms altogether. Where Paper 8 asked "which agricultural transformation does this indicator set detect?", the present study asks "which impact mechanism does this proxy detect?"

This extension is conceptually important because it adds a fourth domain to the framework:

- **Type A** (structural): forest clearance, detected by pollen AP decline
- **Type B** (crop-specific): cultivated taxa, detected by pollen crop indicators
- **Type C** (compositional): tree reorganisation, detected by pollen compositional change
- **Type D** (fire regime): burning intensity/frequency, detected by charcoal accumulation

Types A--C operate within the vegetation domain monitored by pollen. Type D operates in the fire domain monitored by charcoal. A complete assessment of human impact requires coverage of all four domains---which necessarily requires multiple proxies.

### 1.3 The same-core advantage

The strongest test of proxy complementarity uses records from the same sediment core. Same-core pollen and charcoal data share identical chronology, sedimentation history, and catchment area, eliminating the confounds that plague comparisons across sites. If pollen and charcoal from the same core show different exceedance patterns, the difference must reflect genuine differences in what each proxy detects, not differences in site characteristics or chronological control.

The Neotoma Paleoecology Database (Williams et al., 2018) archives both pollen and charcoal records with site-level metadata, enabling systematic identification of same-core pairs at a global scale. Previous multi-proxy studies have been conducted at individual sites or small regional networks (e.g., Whitlock et al., 2007; Colombaroli et al., 2009; Vanniere et al., 2011). No study has systematically quantified pollen-charcoal concordance across a large, globally distributed sample using a uniform detection framework.

### 1.4 Study design and hypotheses

We query the Neotoma database for all sites containing both pollen and charcoal datasets, apply a uniform exceedance detection framework to both proxies independently, and construct a 2 x 2 contingency table classifying each site by the presence or absence of exceedance in each proxy.

We test three hypotheses:

**H1 (Non-independence):** Pollen and charcoal exceedance are statistically associated, because many human impacts produce both vegetation change and fire regime change simultaneously.

**H2 (Incomplete concordance):** The association is significantly less than perfect, because some human impacts operate through only one mechanism (fire without deforestation, or deforestation without fire).

**H3 (Domain specificity):** The discordant cells of the 2 x 2 table (charcoal-only and pollen-only exceedance) represent distinct and interpretable impact mechanisms, not random noise.

If confirmed, these hypotheses would demonstrate that multi-proxy studies are not merely desirable for robustness but structurally necessary for complete impact detection---different proxies recover different fractions of the total human impact signal.

---

## 2. Methods

### 2.1 Data acquisition

We queried the Neotoma Paleoecology Database (accessed March 2026) for all charcoal datasets (dataset type: charcoal), yielding 396 records from 351 unique sites. We then queried all pollen datasets and identified sites where both pollen and charcoal records co-occur. This produced 216 same-site pairs distributed across five broad regions: Europe (n = 49), North America (n = 35), South/Central America (n = 33), Asia/Oceania (n = 5), and unclassified (n = 94).

### 2.2 Temporal filtering

For each same-site pair, we assessed temporal coverage of both the pollen and charcoal records. Sites were retained for exceedance analysis only if they met two criteria: (1) temporal overlap between the pollen and charcoal records, and (2) a pre-baseline period extending beyond 5,000 cal BP, providing sufficient pre-anthropogenic reference data in most regions. Of 216 pairs, 171 showed temporal overlap, and 129 had baselines exceeding 5,000 cal BP. After excluding sites with insufficient sample counts in either record (fewer than 5 samples in the baseline or test period), 111 sites remained for the dual exceedance analysis.

### 2.3 Exceedance framework

We applied a uniform exceedance framework to both pollen and charcoal records independently. For each proxy at each site:

1. **Baseline period**: All samples older than 5,000 cal BP.
2. **Test period**: All samples younger than 5,000 cal BP.
3. **Exceedance metric**: For pollen, total pollen concentration or influx; for charcoal, charcoal concentration or influx (whichever was available).
4. **Threshold**: A site was classified as showing exceedance if > 10% of test-period samples exceeded the 95th percentile of the baseline distribution.
5. **Minimum sample requirement**: At least 5 samples in both baseline and test periods.

This framework is deliberately simple and conservative. The 10% threshold ensures that isolated outliers do not trigger false positives. The 5,000 cal BP boundary provides a broadly applicable pre-anthropogenic baseline for most global regions, though we acknowledge that human impacts predate this boundary in some areas (e.g., the Fertile Crescent, early Holocene Australia).

### 2.4 2 x 2 contingency analysis

Each of the 111 sites was classified into one of four categories based on its exceedance status in each proxy:

- **Both exceedance**: Pollen exceedance AND charcoal exceedance
- **Charcoal only**: Charcoal exceedance AND pollen no signal
- **Pollen only**: Pollen exceedance AND charcoal no signal
- **Neither**: No exceedance in either proxy

We tested for non-independence using Fisher's exact test (appropriate for the modest sample sizes in some cells) and computed the odds ratio as a measure of association strength.

### 2.5 Concordance metrics

We computed two concordance metrics:

- **Overall concordance**: The proportion of sites where both proxies agree (both exceedance + neither), regardless of direction.
- **Discordance rate**: The proportion of sites where the proxies disagree (charcoal only + pollen only).

We also computed region-specific concordance where sample sizes permitted.

---

## 3. Results

### 3.1 Same-site pair identification

Of 351 unique charcoal sites in Neotoma, 216 (61.5%) also contained pollen records from the same site. After temporal filtering (overlap requirement + baseline > 5,000 cal BP), 129 sites were feasible for exceedance analysis. Exclusion of sites with insufficient samples in either proxy reduced the final analytical sample to 111 sites.

The 111 sites span a broad geographic range, with representation from Europe, North and South America, high-latitude regions, and Oceania. The majority of excluded sites failed the temporal overlap criterion due to charcoal records lacking age-depth models (age ranges reported as infinite).

### 3.2 Dual exceedance: the 2 x 2 table

The central result of this study is the 2 x 2 contingency table for the 111 same-core sites (Table 1).

**Table 1.** Dual exceedance classification of 111 same-core pollen-charcoal sites.

|  | Pollen exceedance | Pollen no signal | Row total |
|---|---|---|---|
| **Charcoal exceedance** | 23 (20.7%) | 22 (19.8%) | 45 (40.5%) |
| **Charcoal no signal** | 15 (13.5%) | 51 (45.9%) | 66 (59.5%) |
| **Column total** | 38 (34.2%) | 73 (65.8%) | 111 (100%) |

Fisher's exact test: p = 0.002
Odds ratio: 3.51 (95% CI: 1.53--8.21)
Overall concordance: 74/111 = 66.7%
Discordance rate: 37/111 = 33.3%

### 3.3 Hypothesis testing

**H1 (Non-independence): Supported.** The p-value of 0.002 rejects the null hypothesis of independence. The odds ratio of 3.51 indicates that a site showing charcoal exceedance is 3.5 times more likely to also show pollen exceedance than a site without charcoal exceedance. This confirms that many human (and natural) disturbances produce correlated fire and vegetation signals.

**H2 (Incomplete concordance): Supported.** Despite the significant association, concordance is only 66.7%. One-third of sites show discordant signals---exceedance in one proxy but not the other. If the proxies were redundant, concordance would approach 100% (limited only by noise). The 33.3% discordance demonstrates that each proxy captures a substantial fraction of signals invisible to the other.

**H3 (Domain specificity): Supported.** The discordant cells are not symmetric. Charcoal-only exceedance (22 sites, 19.8%) slightly exceeds pollen-only exceedance (15 sites, 13.5%), suggesting that fire regime changes without lasting vegetation impact are somewhat more common than vegetation changes without fire intensification. Both discordant categories are large enough to be ecologically meaningful and cannot be dismissed as measurement noise.

### 3.4 Marginal detection rates

Charcoal detected exceedance at 45 of 111 sites (40.5%). Pollen detected exceedance at 38 of 111 sites (34.2%). Neither proxy alone captured the full set of disturbed sites. The union of both proxies detected exceedance at 60 sites (54.1%)---a 33% increase over charcoal alone and a 58% increase over pollen alone.

This result quantifies the detection gain from multi-proxy analysis. A study relying solely on pollen would miss 22 sites with fire-regime changes (37% of all charcoal-exceedance sites). A study relying solely on charcoal would miss 15 sites with vegetation changes (39% of all pollen-exceedance sites).

### 3.5 Exceedance magnitude in concordant vs. discordant sites

Among the 23 dual-exceedance sites, charcoal exceedance percentages ranged from 11.1% to 93.3% (median: 24.2%), and pollen exceedance percentages ranged from 14.1% to 92.3% (median: 45.5%). Among charcoal-only sites, charcoal exceedance ranged from 11.8% to 85.7% (median: 23.3%), with pollen values uniformly below the 10% threshold. Among pollen-only sites, pollen exceedance ranged from 11.0% to 84.0% (median: 35.1%), with charcoal values below threshold.

The comparable exceedance magnitudes across categories indicate that discordant detection is not driven by weak signals narrowly missing the threshold. Rather, many discordant sites show strong exceedance in one proxy and near-zero signal in the other, consistent with genuinely different impact mechanisms.

### 3.6 Site-level examples

Several site-level patterns illustrate the ecological interpretation of each cell:

**Dual exceedance.** Tišice (Czech Republic): charcoal 75.0%, pollen 84.8%. This European lowland site shows the classic agricultural signature---simultaneous forest clearance and increased burning associated with Neolithic and later agricultural expansion. Flotatjønn (Norway): charcoal 52.0%, pollen 84.0%. High-latitude pastoral expansion with concurrent burning.

**Charcoal only.** Lake Flåfattjønna: charcoal 84.8%, pollen 0%. An extreme case where fire regime intensified dramatically without detectable vegetation change in the pollen record. Vestre Øykjamyrtjørn (Norway): charcoal 73.8%, pollen 0%. Fire management in a landscape where vegetation composition remained stable despite altered burning regimes. Grizzly Lake (Alaska): charcoal 23.3%, pollen 2.6%. Subarctic fire regime change without vegetation structural transformation.

**Pollen only.** Sa Curcurica (Sardinia): charcoal 0%, pollen 42.4%. Mediterranean vegetation transformation without fire intensification---consistent with arboriculture or grazing-driven change. Lago del Segrino (Italy): charcoal 2.7%, pollen 32.8%. Vegetation reorganisation through non-fire mechanisms. Kinnshaugen (Norway): charcoal 0%, pollen 32.1%. Possible pastoral impact without burning.

**Neither.** Harberton (Tierra del Fuego): charcoal 0%, pollen 0%. Remote high-latitude site with minimal human impact in either domain. Lago Los Niños (Patagonia): charcoal 0%, pollen 0%. Pre-Columbian baseline maintained.

---

## 4. Discussion

### 4.1 Proxy-specific detection domains

The central finding of this study---that one-third of same-core sites show discordant pollen and charcoal signals---formalizes a principle that multi-proxy paleoecologists have long intuited but rarely quantified: different proxies have different detection domains, and human impacts that operate through one mechanism may be systematically invisible to a proxy that monitors a different mechanism.

Pollen monitors the vegetation domain. It detects changes in plant community composition, structure (forest vs. open land), and the presence/absence of specific taxa (including cultivated species). Its detection domain encompasses Types A, B, and C of the three-domain framework developed in Paper 8: structural transformation (forest clearance), crop-specific transformation (cultivated taxa), and compositional transformation (tree reorganisation).

Charcoal monitors the fire domain. It detects changes in fire frequency, intensity, and the spatial extent of burning. Its detection domain is Type D: fire regime change. Charcoal cannot distinguish between natural and anthropogenic fire, but when combined with other evidence (archaeological context, pollen-charcoal discordance), anthropogenic fire management can be inferred.

The key insight is that these domains are partially independent. Some human impacts produce correlated signals in both domains (Neolithic slash-and-burn agriculture: Type A + D). Others produce signals in only one domain (fire management without permanent clearance: Type D only; pastoral deforestation without fire: Type A only; arboriculture replacing wild forest: Type C only).

### 4.2 The 22 charcoal-only sites: fire invisible to pollen

The 22 sites (19.8%) showing charcoal exceedance without pollen exceedance represent the most striking finding. These are locations where fire regime intensified---potentially through anthropogenic burning---but the regional vegetation detected by pollen analysis did not change beyond baseline variability.

Several mechanisms could produce this pattern:

**Low-intensity fire management.** Deliberate burning of understorey vegetation to promote game habitat, facilitate travel, or encourage specific plant resources. Such burning alters fire frequency and charcoal deposition without converting forest to open land. This is well documented ethnographically for Indigenous peoples of North America (Pyne, 1982; Stewart, 2002), Australia (Bowman, 1998), and northern Europe (Hicks, 1985).

**Fire-resistant vegetation.** In some ecosystems, increased burning does not produce lasting vegetation change because the dominant species are fire-adapted. Boreal forests dominated by fire-tolerant conifers can experience altered fire regimes without compositional change detectable in pollen records (Carcaillet et al., 2001).

**Scale mismatch.** Pollen records integrate vegetation signals over a catchment area of several km², while charcoal (especially macroscopic charcoal) reflects more local fire events. Localized anthropogenic burning may appear as charcoal exceedance without producing pollen change at the catchment scale.

**Temporal offset.** Fire regime change may precede vegetation change by decades to centuries. If burning intensifies but has not yet produced permanent vegetation transformation by the end of the record, the site would show charcoal-only exceedance. This interpretation is difficult to test without high-resolution chronological control but is ecologically plausible.

The charcoal-only category is not an artefact of the detection framework. These sites show strong charcoal exceedance (median 23.3%, maximum 85.7%) with near-zero pollen exceedance. The signal is real; it is simply in a domain that pollen does not monitor.

### 4.3 The 15 pollen-only sites: deforestation without fire

The 15 sites (13.5%) showing pollen exceedance without charcoal exceedance represent the converse: vegetation transformation through non-fire mechanisms. Several processes could generate this pattern:

**Pastoral deforestation.** Gradual forest conversion through livestock grazing, which prevents tree regeneration without requiring burning. This is documented in highland and Mediterranean environments (Grove & Rackham, 2001).

**Selective logging and ring-barking.** Mechanical forest removal that deposits no charcoal. Archaeological evidence from various periods demonstrates forest clearance through girdling and felling rather than burning (Iversen, 1956).

**Arboriculture.** Replacement of wild forest with cultivated tree crops (olive, chestnut, walnut) produces strong pollen compositional change without fire. The Mediterranean sites in our dataset (Sa Curcurica, Lago del Segrino) may exemplify this mechanism.

**Climate-driven vegetation change.** Some pollen-only exceedance sites may reflect climatic rather than anthropogenic vegetation change. Without independent archaeological evidence, we cannot exclude this possibility. However, the use of a 5,000 cal BP baseline means that mid-Holocene climatic shifts are captured in the baseline, and late-Holocene exceedance is more likely to reflect human impact.

### 4.4 Extending the identifiability framework: Type D

Paper 8 proposed a three-domain framework for agricultural impact detection within pollen records (Types A, B, C). The present results motivate adding a fourth domain:

**Type D: Fire regime transformation.** Alteration of burning frequency, intensity, or spatial extent. Detected by charcoal concentration/influx exceedance beyond pre-anthropogenic baselines. Type D is independent of vegetation domain changes: fire can intensify without altering vegetation composition (charcoal-only sites), and vegetation can change without fire intensification (pollen-only sites).

The extended four-domain framework provides a more complete classification of human impact mechanisms:

| Domain | Mechanism | Primary proxy | Example |
|---|---|---|---|
| Type A | Forest clearance | Pollen (AP decline) | European Neolithic |
| Type B | Crop cultivation | Pollen (crop indicators) | Eastern Agricultural Complex |
| Type C | Tree reorganisation | Pollen (compositional) | Indigenous fire management |
| Type D | Fire regime change | Charcoal | Deliberate landscape burning |

Critically, Types A and D can co-occur (slash-and-burn agriculture) or occur independently (pastoral clearance without fire; fire management without clearance). The 2 x 2 table quantifies this independence: only 51% of charcoal-exceedance sites also show pollen exceedance, and only 61% of pollen-exceedance sites also show charcoal exceedance. Neither proxy subsumes the other.

### 4.5 Implications for global syntheses

Large-scale syntheses of human impact increasingly rely on single proxies applied globally. The Global Pollen Project (Mottl et al., 2021) and the Global Charcoal Database (Marlon et al., 2016) both provide continent-to-global reconstructions, but they are rarely combined site-by-site because of the analytical challenges of matching records across databases.

Our results suggest that single-proxy syntheses systematically underestimate the spatial extent of human impact. If the 33.3% discordance rate observed in our 111-site sample is representative, then:

- Pollen-only syntheses miss approximately 20% of fire-regime impacts (the charcoal-only fraction).
- Charcoal-only syntheses miss approximately 14% of vegetation impacts (the pollen-only fraction).
- Only multi-proxy approaches capture the full 54% of sites showing exceedance in at least one domain.

These are conservative estimates because our exceedance framework uses a high threshold (> 10% of samples exceeding the 95th percentile). More sensitive detection methods would likely find additional discordant sites.

### 4.6 The odds ratio as a measure of proxy coupling

The odds ratio of 3.51 quantifies the degree of coupling between the fire and vegetation domains. An OR of 1.0 would indicate complete independence; an OR approaching infinity would indicate perfect coupling. The observed value of 3.51 indicates moderate coupling: fire and vegetation change are associated but far from deterministically linked.

This intermediate coupling is ecologically sensible. Fire is one of several mechanisms that produce vegetation change, and vegetation change is one of several consequences of fire. The OR of 3.51 may be interpreted as the "cross-domain leakage"---the degree to which impact in one domain tends to produce detectable impact in the other. The substantial residual discordance (33.3%) confirms that this leakage is far from complete.

### 4.7 Limitations

**Baseline assumption.** The 5,000 cal BP baseline may be inappropriate for regions with earlier human impact (Near East, East Asia) or later colonization (Pacific Islands, New Zealand). Sites in these regions may show misclassified exceedance or missed early impacts. Future analyses should employ region-specific baselines.

**Charcoal data quality.** Neotoma's charcoal records vary in measurement type (concentration, influx, counts, area), preparation method, and particle size fraction. We applied a uniform framework across all measurement types, which may introduce heterogeneity. Records lacking age-depth models (which produced the infinite age ranges visible in our temporal filtering) further reduced the available sample.

**Pollen metric.** We used total pollen concentration/influx rather than compositional metrics (AP/NAP ratio, specific indicator taxa). This captures gross vegetation productivity changes but may miss compositional transformations (Type B, Type C) that do not alter total pollen output. A future extension could apply the full three-domain pollen framework of Paper 8 to the same-core sites.

**Sample size.** The 111-site sample, while the largest same-core pollen-charcoal comparison to date, remains modest for region-specific analyses. The geographic distribution is heavily weighted toward Europe and the Americas, with minimal representation from Asia, Africa, and Oceania.

**Natural vs. anthropogenic drivers.** Neither pollen nor charcoal exceedance can be attributed to human impact without independent evidence. Climate-driven fire regime changes and vegetation shifts could produce exceedance patterns indistinguishable from anthropogenic impacts in our framework. The 2 x 2 table characterizes proxy discordance, not the human vs. natural attribution of that discordance.

### 4.8 Recommendations for multi-proxy practice

Based on these results, we offer four methodological recommendations:

1. **Report proxy-specific detection domains.** Every paleoecological study of human impact should explicitly state which impact domains its chosen proxy can and cannot detect. A pollen study should acknowledge its blindness to Type D (fire regime) changes; a charcoal study should acknowledge its blindness to Types A, B, and C (vegetation domain) changes.

2. **Quantify the detection gap.** Where same-core multi-proxy data exist, studies should report the concordance and discordance rates between proxies, providing an empirical estimate of what each proxy misses.

3. **Prioritize same-core multi-proxy records.** For questions about the total magnitude of human impact, same-core pollen-charcoal records provide the strongest evidence because they eliminate inter-site confounds. Database efforts should prioritize archiving and linking multi-proxy records from the same cores.

4. **Extend to additional proxies.** The framework developed here for pollen vs. charcoal extends naturally to other proxy combinations: pollen vs. non-pollen palynomorphs (coprophilous fungi for pastoralism), charcoal vs. geochemical proxies (lead, phosphorus for metallurgy and agriculture), or any combination where proxies monitor different ecological domains.

---

## 5. Conclusions

1. **Significant but incomplete concordance.** Pollen and charcoal exceedance at 111 same-core sites are significantly associated (Fisher's exact p = 0.002, OR = 3.51) but show only 66.7% concordance. One-third of sites show exceedance in one proxy but not the other.

2. **Charcoal detects fire-domain impacts invisible to pollen.** Twenty-two sites (19.8%) show charcoal exceedance without pollen exceedance, representing fire regime changes---potentially anthropogenic burning---that did not produce lasting vegetation transformation detectable in pollen records.

3. **Pollen detects vegetation-domain impacts invisible to charcoal.** Fifteen sites (13.5%) show pollen exceedance without charcoal exceedance, representing vegetation transformations through non-fire mechanisms such as pastoralism, selective logging, or arboriculture.

4. **Multi-proxy studies are structurally necessary.** The union of both proxies detects exceedance at 54.1% of sites, compared to 40.5% for charcoal alone and 34.2% for pollen alone. Single-proxy studies systematically underestimate the spatial extent of human impact by missing impacts in domains they do not monitor.

5. **The identifiability framework extends to inter-proxy detection.** We propose adding Type D (fire regime domain, monitored by charcoal) to the three-domain pollen framework (Types A--C), creating a four-domain classification of human impact mechanisms. Complete identifiability requires proxy coverage of all four domains.

---

## Data availability

All pollen and charcoal data used in this study are publicly available through the Neotoma Paleoecology Database (https://www.neotomadb.org/). Site lists, exceedance results, and analysis scripts are archived in the project repository.

---

## References

Abrams, M. D., & Nowacki, G. J. (2008). Native Americans as active and passive promoters of mast and fruit trees in the eastern USA. *The Holocene*, 18(7), 1065--1077.

Behre, K.-E. (1981). The interpretation of anthropogenic indicators in pollen diagrams. *Pollen et Spores*, 23, 225--245.

Bowman, D. M. J. S. (1998). The impact of Aboriginal landscape burning on the Australian biota. *New Phytologist*, 140(3), 385--410.

Bowman, D. M. J. S., Balch, J., Artaxo, P., Bond, W. J., Cochrane, M. A., D'Antonio, C. M., ... Pyne, S. J. (2011). The human dimension of fire regimes on Earth. *Journal of Biogeography*, 38(12), 2223--2236.

Carcaillet, C., Bergeron, Y., Richard, P. J. H., Frechette, B., Gauthier, S., & Prairie, Y. T. (2001). Change of fire frequency in the eastern Canadian boreal forests during the Holocene: Does vegetation composition or climate trigger the fire regime? *Journal of Ecology*, 89(6), 930--946.

Colombaroli, D., Tinner, W., van Leeuwen, J., Noti, R., Vescovi, E., Vanniere, B., ... Saurer, M. (2009). Response of broadleaved evergreen Mediterranean forest vegetation to fire disturbance during the Holocene: Insights from the peri-Adriatic region. *Journal of Biogeography*, 36(2), 314--326.

Ellis, E. C., Gauthier, N., Goldewijk, K. K., Bird, R. B., Boivin, N., Diaz, S., ... Watson, J. E. M. (2021). People have shaped most of terrestrial nature for at least 12,000 years. *Proceedings of the National Academy of Sciences*, 118(17), e2023483118.

Fyfe, R. M., Woodbridge, J., & Roberts, N. (2015). From forest to farmland: Pollen-inferred land cover change across Europe using the pseudobiomization approach. *Global Change Biology*, 21(3), 1197--1212.

Grove, A. T., & Rackham, O. (2001). *The Nature of Mediterranean Europe: An Ecological History*. Yale University Press.

Hicks, S. (1985). Modern pollen deposition records from Kuusamo, Finland. *Grana*, 24(3), 167--184.

Iversen, J. (1956). Forest clearance in the Stone Age. *Scientific American*, 194(3), 36--41.

Marlon, J. R., Bartlein, P. J., Carcaillet, C., Gavin, D. G., Harrison, S. P., Higuera, P. E., ... Power, M. J. (2008). Climate and human influences on global biomass burning over the past two millennia. *Nature Geoscience*, 1(10), 697--702.

Marlon, J. R., Kelly, R., Daniau, A.-L., Vanniere, B., Power, M. J., Bartlein, P., ... Tinner, W. (2016). Reconstructions of biomass burning from sediment-charcoal records to improve data-model comparisons. *Biogeosciences*, 13, 3225--3244.

Mottl, O., Flantua, S. G. A., Bhatt, S., Felde, V. A., Giesecke, T., Goring, S., ... Williams, J. W. (2021). Global acceleration in rates of vegetation change over the past 18,000 years. *Science*, 372(6544), 860--864.

Pyne, S. J. (1982). *Fire in America: A Cultural History of Wildland and Rural Fire*. Princeton University Press.

Roberts, N., Fyfe, R. M., Woodbridge, J., Gaillard, M.-J., Davis, B. A. S., Kaplan, J. O., ... Leydet, M. (2018). Europe's lost forests: A pollen-based synthesis for the last 11,000 years. *Scientific Reports*, 8, 716.

Stewart, O. C. (2002). *Forgotten Fires: Native Americans and the Transient Wilderness*. University of Oklahoma Press.

Vanniere, B., Power, M. J., Roberts, N., Tinner, W., Carrion, J., Magny, M., ... Sadori, L. (2011). Circum-Mediterranean fire activity and climate changes during the mid-Holocene environmental transition (8500--2500 cal. BP). *The Holocene*, 21(1), 53--73.

Whitlock, C., & Larsen, C. (2001). Charcoal as a fire proxy. In J. P. Smol, H. J. B. Birks, & W. M. Last (Eds.), *Tracking Environmental Change Using Lake Sediments*, Vol. 3: Terrestrial, Algal, and Siliceous Indicators (pp. 75--97). Kluwer.

Whitlock, C., Higuera, P. E., McWethy, D. B., & Briles, C. E. (2007). Paleoecological perspectives on fire ecology: Revisiting the fire-regime concept. *The Open Ecology Journal*, 3, 6--23.

Williams, J. W., Grimm, E. C., Blois, J. L., Charles, D. F., Davis, E. B., Goring, S. J., ... Takahara, H. (2018). The Neotoma Paleoecology Database, a multiproxy, international, community-curated data resource. *Quaternary Research*, 89(1), 156--177.

---

**Word count**: ~5,800

---

*Manuscript prepared: 2026-03-28*
