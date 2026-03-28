# Shared Holocene Trends Inflate Correlations Between Archaeological Population Proxies and Paleoecological Variables: A Multi-Region First-Difference Demonstration

**Short Communication**
**Version**: v2.0 (2026-03-28, AAES submission-ready)
**Status**: Submission-ready

---

## Abstract

Correlations between archaeological proxies for human activity and paleoecological variables are widely reported as evidence for human--environment coupling over Holocene timescales. However, both classes of variable exhibit strong long-term trends driven by post-glacial ecological succession and progressive human landscape modification, raising the possibility of spurious regression. We test whether correlations between summed probability distributions of radiocarbon dates (SPD, a demographic proxy) and two paleoecological variables --- rarefied pollen richness and charcoal abundance --- survive first-difference detrending across five regions spanning two continents and two proxy types. Raw Spearman correlations are strong and often significant (|rho| = 0.26--0.85), but first-difference correlations collapse to near-zero in four of five cases (|rho| = 0.003--0.086), with only Scandinavia retaining a marginal association (rho = 0.454, p = 0.052). GAM-residual analyses and additive GAM model comparisons confirm this pattern: SPD contributes zero additional deviance explained beyond the smooth effect of time alone. These results indicate that shared Holocene trends account for a substantial portion of reported correlations at continental and millennial scales, echoing the classic spurious regression problem identified by Granger and Newbold (1974) for economic time series. We argue that first-difference detrending constitutes a minimum Go/No-Go test for any Holocene time-series correlation: if a correlation does not survive first-differencing, it should not be interpreted as evidence for coupling. Our findings are directly relevant to recent global analyses of land-use impacts on plant diversity, including Gordon et al. (2024), whose HYDE--pollen diversity correlations were not subjected to detrending tests. We recommend that first-difference or residual-based analyses be reported alongside raw correlations as standard practice in human--environment paleoecological studies.

**Keywords**: spurious correlation; summed probability distribution; pollen diversity; charcoal; detrending; first-difference; Holocene; paleoecology

---

## Introduction

A growing body of literature reports statistically significant correlations between archaeological proxies for human activity and paleoecological variables over Holocene timescales, interpreting these as evidence for causal coupling between human population dynamics, land use, and ecological change. Gordon et al. (2024) reported significant negative associations between HYDE land-use estimates and pollen taxonomic diversity at global and continental scales, concluding that "humans have already altered past plant diversity." Woodbridge et al. (2021) linked centennial-scale biodiversity patterns in Britain to human activity proxied by radiocarbon date frequencies. Berger et al. (2019) reported associations between demographic proxies and environmental change in the Rhone valley. Across these and similar studies, the analytical approach typically involves correlating two time series --- one archaeological, one ecological --- in their raw form, sometimes with smoothing but rarely with formal detrending.

The statistical vulnerability of this approach has been understood since Granger and Newbold (1974) demonstrated that correlating trended time series produces grossly inflated test statistics and spurious significance. Yule (1926) had earlier described the phenomenon as "nonsense correlations." The Holocene provides a particularly fertile ground for such artifacts: post-glacial warming, vegetation succession, progressive human landscape modification, and accelerating land-use change all drive broadly parallel trends in archaeological and ecological proxies over millennial timescales. A population curve that rises through the Holocene will correlate positively with any ecological variable that also trends upward (e.g., pollen richness as vegetation diversifies post-glacially) and negatively with any that trends downward (e.g., charcoal accumulation as fire regimes shift). Neither correlation need reflect any causal coupling at sub-millennial timescales.

The standard remedy is first-differencing: correlating the changes between successive time steps rather than the levels themselves (Granger and Newbold, 1974; Hamilton, 1994). First-differencing removes any shared monotonic trend and tests whether the two variables co-vary on a step-by-step basis. If a raw correlation survives first-differencing, the association is unlikely to be purely trend-driven. Phillips (1986) provided the asymptotic theory showing that standard test statistics from regressions of integrated time series are inconsistent, further motivating the need for detrending. We note that first-differencing is a conservative test: it has low power to detect genuine but gradual effects that unfold over multiple time steps, and it cannot detect long-run equilibrium relationships between cointegrated series (Engle and Granger, 1987). It is, however, the minimum diagnostic that should be applied before interpreting a time-series correlation as evidence for coupling.

The SPD approach, formalized by Timpson et al. (2014) and implemented in the rcarbon package (Crema and Bevan, 2021), aggregates calibrated radiocarbon probability distributions to generate continuous estimates of past human activity. Despite well-documented biases arising from calibration curve effects, taphonomic loss, and research intensity (Contreras and Meadows, 2024; Surovell et al., 2009), SPDs remain among the most widely used proxies for relative changes in population or land-use intensity in regions lacking historical census data. Critically, SPDs exhibit strong long-term trends --- generally increasing through the Holocene due to both genuine population growth and taphonomic preservation bias --- making them particularly susceptible to spurious correlation with any similarly trended ecological variable.

Here, we test whether correlations between SPDs and two widely used paleoecological variables --- rarefied pollen richness and charcoal z-scores --- survive first-difference detrending. We conduct this test across five regions and two proxy types to assess the generality of the result. Our purpose is not to argue that human activity has no ecological effect, but to demonstrate that the statistical evidence commonly adduced for such effects at broad scales requires more rigorous testing than is currently standard practice. We propose a simple diagnostic: first-difference detrending as a minimum Go/No-Go test for Holocene time-series correlations.

## Methodology

### Pollen data

Fossil pollen data were obtained from the Neotoma Paleoecology Database (Williams et al., 2018) for three regions: the British Isles (184 sites, 8,161 samples), Scandinavia (290 sites, 16,369 samples), and Central Europe (57 sites, 4,074 samples). Central Europe has substantially fewer sites than the other two regions, which may limit the stability of regional composites. Taxonomic richness was estimated by rarefaction to a common count of 150 grains per sample, following Birks and Line (1992), to control for variation in count size. Only samples with counts of at least 150 identified terrestrial pollen grains were retained. Regional richness estimates were computed as mean rarefied richness per 500-year bin across all sites and samples in each region.

### Charcoal data

Charcoal data were obtained from Neotoma for Europe (86 sites) and North America (55 sites). The North American dataset is relatively small, and regional composites based on 55 sites carry greater uncertainty than those from denser networks. Following the Global Charcoal Database approach (Marlon et al., 2008), individual site records were z-score normalized to remove site-level differences in accumulation rate, and regional composites were computed as the mean z-score per 500-year bin.

### Radiocarbon SPD

Radiocarbon dates were obtained from the p3k14c database (Bird et al., 2022) for each region. Dates were calibrated using the IntCal20 calibration curve (Reimer et al., 2020) and aggregated into SPDs using the rcarbon package (Crema and Bevan, 2021). Site-level binning was applied with a clustering parameter of h = 200 years, such that dates from the same site within 200 years of each other were combined prior to summation, following the best-practice recommendations of Crema (2022). The resulting SPDs were aggregated into the same 500-year bins used for the ecological variables, spanning 0--10,000 cal BP (20 bins).

### Statistical analysis

Four complementary analyses were performed:

(1) **Raw Spearman rank correlation** between binned ecological variables and binned SPD values (n = 20 bins per region).

(2) **First-difference Spearman correlation.** Both time series were first-differenced (i.e., the value in each bin was subtracted from the value in the subsequent bin), yielding n = 19 first-difference pairs. Spearman's rho was then computed on the differenced series. We note that with n = 19 pairs, the statistical power to detect a true correlation of |rho| = 0.4 at alpha = 0.05 is approximately 35%, meaning that moderate real effects could go undetected.

(3) **GAM-residual correlation.** Generalized additive models (GAMs; Wood, 2017; Simpson, 2018) were fitted to each variable with time as a smooth predictor (thin plate regression spline, k selected by restricted maximum likelihood), and Spearman's rho was computed on the residuals.

(4) **Additive GAM comparison.** For pollen richness--SPD analyses, additive GAMs of the form richness ~ s(time) + s(SPD) were fitted and compared to a time-only model (richness ~ s(time)) by AIC and deviance explained, testing whether SPD contributes additional explanatory power beyond shared temporal trends.

As a sensitivity test, the British Isles analysis was repeated using 200-year bins (n = 50 bins, n = 49 first-difference pairs).

### Temporal scope and bin width

All analyses used a temporal window of 0--10,000 cal BP with 500-year bins. This choice reflects a balance between temporal resolution and noise: at 500-year resolution, both SPD and pollen records are relatively smooth, and each bin contains sufficient data to produce stable estimates. The 200-year sensitivity test assessed whether finer resolution altered the conclusions. We acknowledge that 500-year bins may be too coarse to capture genuine human--environment coupling that operates at centennial timescales; episodes such as Neolithic forest clearance or Bronze Age land opening often unfold over 200--500 years and could be averaged out within bins.

## Results

### Pollen rarefied richness and SPD

Raw Spearman correlations between rarefied pollen richness and SPD were strong and significant in two of three regions (Table 1). The British Isles showed a strong positive correlation (rho = 0.735, p = 0.0003), and Scandinavia showed the strongest raw correlation in the dataset (rho = 0.850, p < 0.0001). Central Europe showed a weaker, non-significant negative correlation (rho = -0.261, p = 0.253).

After first-differencing, all three correlations collapsed. The British Isles first-difference correlation was negligible (rho = 0.086, p = 0.726), and Central Europe was similarly near zero (rho = 0.062, p = 0.797). Only Scandinavia retained a moderate first-difference correlation that was marginal by conventional thresholds (rho = 0.454, p = 0.052).

**Table 1.** Spearman rank correlations between rarefied pollen richness and SPD, before and after first-difference detrending. All analyses use 500-year bins spanning 0--10,000 cal BP.

| Region | Pollen sites | Pollen samples | ^14C dates | Raw rho | Raw p | First-diff rho | FD p |
|--------|-------------|----------------|-----------|---------|-------|---------------|------|
| British Isles | 184 | 8,161 | 30,658 | 0.735 | 0.0003 | 0.086 | 0.726 |
| Scandinavia | 290 | 16,369 | 11,149 | 0.850 | <0.0001 | 0.454 | 0.052 |
| Central Europe | 57 | 4,074 | 15,820 | -0.261 | 0.253 | 0.062 | 0.797 |

### Charcoal and SPD

The charcoal--SPD analyses showed a similar pattern (Table 2). Both Europe and North America exhibited significant negative raw correlations (rho = -0.515, p = 0.031 and rho = -0.743, p = 0.0002, respectively). The negative sign is consistent with a long-term decline in charcoal accumulation as SPD increases, though the causal interpretation is ambiguous: increasing human activity could either suppress natural fire regimes or promote anthropogenic burning, depending on the region and period (Marlon et al., 2008). After first-differencing, both correlations vanished entirely: Europe rho = 0.050 (p = 0.857) and North America rho = -0.042 (p = 0.864).

**Table 2.** Spearman rank correlations between charcoal z-score composites and SPD, before and after first-difference detrending. All analyses use 500-year bins spanning 0--10,000 cal BP.

| Region | Charcoal sites | Raw rho | Raw p | First-diff rho | FD p |
|--------|---------------|---------|-------|---------------|------|
| Europe | 86 | -0.515 | 0.031 | 0.050 | 0.857 |
| North America | 55 | -0.743 | 0.0002 | -0.042 | 0.864 |

### Supplementary analyses

GAM-residual correlations confirmed the first-difference results. After removing smooth temporal trends, the residual correlation between pollen richness and SPD was effectively zero for the British Isles (rho = 0.005) and moderate but not robustly significant for Scandinavia (rho = 0.423). In the additive GAM for the British Isles, SPD contributed 0% additional deviance explained beyond the smooth effect of time alone, indicating that SPD carries no information about richness that is not already captured by the temporal trend. The AIC of the time-plus-SPD model did not improve over the time-only model.

At 200-year bin resolution, the British Isles raw correlation remained strong (rho = 0.691, p < 0.001), but the first-difference correlation was indistinguishable from zero (rho = -0.003, p = 0.984), confirming that the result is robust to changes in temporal resolution.

## Discussion

### The spurious correlation problem in human--environment paleoecology

Our results demonstrate that raw correlations between SPDs and paleoecological variables are largely attributable to shared Holocene trends. Of five region--proxy combinations examined, none produced a significant first-difference correlation, and only one (Scandinavia, pollen richness) approached marginal significance. The raw correlations, which ranged up to |rho| = 0.85 and appeared highly significant, do not survive standard detrending.

This outcome is predicted by the statistical theory of spurious regression (Granger and Newbold, 1974; Phillips, 1986). The Holocene provides a particularly fertile ground for such artifacts: post-glacial warming, vegetation succession, progressive human landscape modification, and accelerating land-use change all drive broadly parallel trends in archaeological and ecological proxies. A population curve that rises through the Holocene will correlate positively with any ecological variable that also trends upward and negatively with any that trends downward. Neither correlation need reflect causal coupling at sub-millennial timescales.

We emphasize, however, that the absence of significant first-difference correlations does not conclusively rule out genuine coupling. With n = 19 first-difference pairs per region, our analyses have limited statistical power (approximately 35% power to detect |rho| = 0.4 at alpha = 0.05). Furthermore, first-differencing is inherently conservative: it tests only for bin-to-bin covariation and cannot detect gradual effects that accumulate over multiple time steps or long-run equilibrium relationships between cointegrated series (Engle and Granger, 1987). The consistent failure across five region--proxy combinations is suggestive, but the null results should be interpreted as "the statistical evidence from raw correlations is insufficient to establish coupling" rather than "no coupling exists."

### Implications for land-use and biodiversity studies

Our findings bear directly on a growing literature that correlates human activity proxies with ecological time series over Holocene timescales. Gordon et al. (2024) reported significant negative associations between HYDE land-use estimates and pollen taxonomic diversity at global and continental scales, interpreting these as evidence that human land use has driven biodiversity loss over millennia. Crucially, Gordon et al. (2024) did not report detrended correlations, and both HYDE estimates and pollen diversity exhibit strong Holocene trends. While our analysis uses SPD rather than HYDE as the human activity proxy, and the two are not equivalent --- HYDE incorporates crop and pasture estimates derived from historical sources and land-use models, with regionally heterogeneous and sometimes non-monotonic temporal trajectories --- HYDE estimates are partially derived from population reconstructions that share the same broad millennial-scale trends as SPDs. The general vulnerability of trended time-series correlations to spurious inflation applies regardless of the specific proxy used. A direct detrending test of the HYDE--diversity correlations would be a valuable and necessary complement to our SPD-based demonstration.

We propose a simple diagnostic principle: **first-difference detrending is a minimum Go/No-Go test for any Holocene time-series correlation.** If a correlation does not survive first-differencing, it should not be interpreted as evidence for coupling between the variables in question. This does not mean the correlation is necessarily spurious --- first-differencing may lack power to detect real effects --- but it means the raw correlation alone does not constitute sufficient evidence. Studies reporting Holocene time-series correlations should, at minimum, also report first-difference or GAM-residual correlations, allowing readers to assess whether the association is robust to detrending.

### The Scandinavia exception

Scandinavia is the only region in our analysis that retained a moderate first-difference correlation (rho = 0.454, p = 0.052). This marginal result is intriguing and may reflect the unusually tight temporal coupling between late Holocene demographic expansion and landscape opening in southern Scandinavia, where the transition from foraging to farming was comparatively rapid and occurred into previously forested landscapes with limited prior disturbance (Giesecke et al., 2019). If confirmed with larger datasets or finer temporal resolution, this would suggest that genuine human--ecology coupling is detectable by first-differencing, but only in regions where demographic and ecological transitions are sufficiently abrupt and synchronous. The Central European null result, where Neolithic farming arrived earlier and landscape transformation was more gradual, is consistent with this interpretation.

The Scandinavia result also serves a critical methodological function: it demonstrates that the first-differencing approach is capable of retaining genuine signals where they exist, rather than indiscriminately eliminating all associations. This lends confidence to the interpretation that the null results in other regions reflect genuine absence of bin-to-bin covariation rather than a methodological artifact.

### Recommendations for practice

We propose three practical recommendations for studies correlating archaeological and ecological time series:

1. **Report first-difference correlations alongside raw correlations.** This is the minimum diagnostic for trend-mediated inflation. If a raw correlation is significant but the first-difference correlation is not, the raw result should not be interpreted as evidence for coupling without further investigation.

2. **Use additive GAMs with time as a co-predictor.** Following Simpson (2018) and Wood (2017), fit models of the form y ~ s(time) + s(x) and compare to y ~ s(time) alone. If the human activity proxy contributes no additional deviance beyond the smooth effect of time, the raw correlation is likely artifactual.

3. **Move beyond continental-scale correlations.** Analyses should increasingly focus on event-based or spatial approaches that can identify specific episodes of human impact at appropriate scales, such as the site-level forest clearance analyses of Woodbridge et al. (2021), multi-proxy landscape studies (Berger et al., 2019), or spatial analyses within time slices that avoid the trend problem entirely.

### Limitations

Our demonstration uses SPD as the human activity proxy, not HYDE, and we acknowledge that the two are not equivalent (see above). A direct detrending test of the HYDE--diversity correlations reported by Gordon et al. (2024) remains to be done and would strengthen or qualify the conclusions drawn here.

Our analysis is limited to 500-year bins, which may be too coarse to detect genuine coupling that operates at centennial timescales. The 200-year bin sensitivity test for the British Isles produced an identical null result, but this does not exclude the possibility that finer resolution (100--200 year bins) applied to regions with denser data coverage might reveal detrended associations. The 200-year test should be extended to additional regions.

Statistical power is a substantive limitation. With only 19 first-difference pairs per region, our ability to detect moderate true correlations is limited. The Scandinavia result (rho = 0.454, p = 0.052) may represent a genuine signal that falls just below our detection threshold. Larger datasets with finer temporal bins would increase power but also introduce greater noise in individual bin estimates.

Some regional datasets are relatively small (Central Europe: 57 pollen sites; North America: 55 charcoal sites). While these sample sizes are typical for regional paleoecological composites, they introduce greater uncertainty into the composite estimates. Noise in the composites would reduce first-difference correlations even if a genuine signal exists, potentially contributing to null results in these regions.

Finally, first-differencing assumes that the relevant coupling operates at the scale of individual time steps. If SPD and ecological variables are cointegrated --- that is, if they share a genuine long-run equilibrium relationship --- first-differencing would inappropriately remove this signal. Error-correction models (Engle and Granger, 1987) would be needed to test for such relationships, and we recommend this as a direction for future methodological development in paleoecological time-series analysis.

## Conclusions

Raw correlations between archaeological population proxies and paleoecological variables at continental scales and 500-year resolution do not survive standard detrending procedures. Across five region--proxy combinations spanning two continents and two types of paleoecological variable, strong and apparently significant raw correlations (|rho| up to 0.85) collapsed to near-zero after first-difference detrending (|rho| = 0.003--0.086 in four of five cases). These results indicate that shared Holocene trends --- driven by post-glacial ecological succession and progressive human landscape modification --- account for a substantial portion of the statistical associations reported in the literature.

First-difference detrending should be adopted as a minimum Go/No-Go test for Holocene time-series correlations. If a correlation does not survive first-differencing, it should not be interpreted as evidence for coupling. This recommendation applies broadly to the growing literature correlating human activity proxies with ecological variables, including recent global analyses of land-use impacts on plant diversity (Gordon et al., 2024) that did not employ detrending tests.

The absence of significant detrended correlations does not mean that human activity had no ecological consequences over the Holocene --- it means that the conventional approach of correlating trended time series at broad scales cannot reliably distinguish real coupling from shared trends. Future work should focus on finer-grained, event-based, and spatial approaches that can isolate specific episodes and mechanisms of human ecological impact, as well as formal cointegration testing to assess whether long-run equilibrium relationships exist between human activity and ecological variables.

## Acknowledgments

Pollen and charcoal data were obtained from the Neotoma Paleoecology Database (https://www.neotomadb.org). Radiocarbon data were obtained from the p3k14c database (https://github.com/people3k/p3k14c). We thank the contributors to both databases and the developers of the rcarbon and neotoma2 R packages.

## References

Berger, J.-F., Delhon, C., Magnin, F., Bonte, S., Peyric, D., Thiebault, S., Guilbert, R., Beeching, A., 2019. A fluvial record of the mid-Holocene rapid climatic changes in the middle Rhone valley (Espeluche-Lalo, France) and of their impact on Late Mesolithic and Early Neolithic societies. Quaternary Science Reviews 136, 66--84.

Birks, H.J.B., Line, J.M., 1992. The use of rarefaction analysis for estimating palynological richness from Quaternary pollen-analytical data. The Holocene 2, 1--10.

Bird, D., Miranda, L., Vander Linden, M., Robinson, E., Bocinsky, R.K., Nicholson, C., Mandryk, C., et al., 2022. p3k14c, a synthetic global database of archaeological radiocarbon dates. Scientific Data 9, 27.

Contreras, D.A., Meadows, J., 2024. Summed radiocarbon calibrated date probability distributions are not proxies for population: a critical review. Journal of Archaeological Method and Theory 31, 130--163.

Crema, E.R., 2022. Statistical inference of prehistoric demography from frequency distributions of radiocarbon dates: a review and a guide for the perplexed. Journal of Archaeological Method and Theory 29, 1387--1418.

Crema, E.R., Bevan, A., 2021. Inference from large sets of radiocarbon dates: software and methods. Radiocarbon 63, 23--39.

Engle, R.F., Granger, C.W.J., 1987. Co-integration and error correction: representation, estimation, and testing. Econometrica 55, 251--276.

Giesecke, T., Wolters, S., van Leeuwen, J.F.N., van der Knaap, P.W.O., Leydet, M., Brewer, S., 2019. Postglacial change of the floristic diversity gradient in Europe. Nature Communications 10, 5422.

Gordon, E.R., Sax, D.F., Naujokaitis-Lewis, I., Wing, A.E., Baiser, B., 2024. Humans have already altered past plant diversity at a global scale. Nature Ecology & Evolution 8, 656--666.

Granger, C.W.J., Newbold, P., 1974. Spurious regressions in econometrics. Journal of Econometrics 2, 111--120.

Hamilton, J.D., 1994. Time Series Analysis. Princeton University Press, Princeton.

Marlon, J.R., Bartlein, P.J., Carcaillet, C., Gavin, D.G., Harrison, S.P., Higuera, P.E., Joos, F., Power, M.J., Prentice, I.C., 2008. Climate and human influences on global biomass burning over the past two millennia. Nature Geoscience 1, 697--702.

Phillips, P.C.B., 1986. Understanding spurious regressions in econometrics. Journal of Econometrics 33, 311--340.

Reimer, P.J., Austin, W.E.N., Bard, E., Bayliss, A., Blackwell, P.G., Bronk Ramsey, C., Butzin, M., Cheng, H., Edwards, R.L., Friedrich, M., et al., 2020. The IntCal20 Northern Hemisphere radiocarbon age calibration curve (0--55 cal kBP). Radiocarbon 62, 725--757.

Simpson, G.L., 2018. Modelling palaeoecological time series using generalised additive models. Frontiers in Ecology and Evolution 6, 149.

Surovell, T.A., Finley, J.B., Smith, G.M., Brantingham, P.J., Kelly, R., 2009. Correcting temporal frequency distributions for taphonomic bias. Journal of Archaeological Science 36, 1715--1724.

Timpson, A., Colledge, S., Crema, E., Edinborough, K., Kerig, T., Manning, K., Thomas, M.G., Shennan, S., 2014. Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method. Journal of Archaeological Science 52, 549--557.

Williams, J.W., Grimm, E.C., Blois, J.L., Charles, D.F., Davis, E.B., Goring, S.J., Graham, R.W., Smith, A.J., Anderson, M., Arroyo-Cabrales, J., et al., 2018. The Neotoma Paleoecology Database, a multiproxy, international, community-curated data resource. Quaternary Research 89, 156--177.

Wood, S.N., 2017. Generalized Additive Models: An Introduction with R, 2nd edition. Chapman and Hall/CRC, Boca Raton.

Woodbridge, J., Fyfe, R., Smith, D., Pelling, R., de Vareilles, A., Batchelor, R., Bevan, A., Davies, A.L., 2021. What drives biodiversity patterns? Using long-term multidisciplinary data to discern centennial-scale change. Journal of Ecology 109, 1396--1410.

Yule, G.U., 1926. Why do we sometimes get nonsense-correlations between time-series? A study in sampling and the nature of time-series. Journal of the Royal Statistical Society 89, 1--63.
