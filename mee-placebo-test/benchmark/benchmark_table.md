# Benchmark Summary: Diagnostic Outcomes Across Analyses

**Table 3.** Summary of diagnostic outcomes across all analyses in this study. Each row represents an independent analysis of citizen science data subjected to our diagnostic framework. "Collapse" indicates the diagnostic revealed the main result as likely artifactual; "Survive" indicates the main result withstood diagnostic scrutiny; "Partial" indicates some components survived while others collapsed.

| # | Analysis | Platform | Taxon | Data use (risk level) | Diagnostic type | Main effect | Diagnostic effect | Ratio | Outcome |
|---|----------|----------|-------|-----------------------|-----------------|-------------|-------------------|-------|---------|
| 1 | ALAN → dawn chorus timing | Xeno-canto | Birds (5 spp.) | Timestamps (Very high) | Type A: midday window | r = −0.142*** | r = −0.208*** | 1.46 | **Collapse** |
| 2 | ALAN → nocturnal insect richness | GBIF | Lepidoptera | Richness along gradient (High) | Type A: diurnal taxon | β = 1.32* | β = 3.21*** | 2.43 | **Collapse** |
| 3a | Forest fragmentation → bird richness (uncorrected) | GBIF | Birds | Richness along gradient (High) | Type A+C: diagnostic taxa | Significant** | Significant** | ~1.0 | **Collapse** (pre-correction) |
| 3b | Forest fragmentation → bird richness (corrected) | GBIF | Birds | Richness along gradient (High) | Type C: effort correction | p = 0.018* | p = 0.075–0.092 ns | <0.5 | **Survive** |
| 4 | Latitude → pollinator phenological shift | GBIF | Apidae | Disaggregated temporal (Mod–High) | Type B: effort correction | R² = 0.79** | R² < 0.05 ns | <0.1 | **Collapse** (gradient) |
| 4b | Overall pollinator phenological shift | GBIF | Apidae | Aggregated temporal (Moderate) | Type B: effort correction | −3.4 d/yr*** | Survives correction | — | **Survive** |
| 5 | ALAN → bird richness (GBIF) | GBIF | Birds | Richness along gradient (High) | Type A: diurnal birds | β_main | β = 2.61*** | >1.0 | **Collapse** |
| 5b | ALAN → bird richness (eBird) | eBird | Birds | Richness along gradient (High) | Type A: diurnal birds | β = 0.66*** | β = 2.61*** | 3.95 | **Collapse** (reduced) |
| 6 | Spatial occurrence → SDM | GBIF | Multiple | Spatial occurrence (Low) | — (boundary) | AUC 0.67–0.98 | Not required | — | **Survive** |

\* p < 0.05; \*\* p < 0.01; \*\*\* p < 0.001; ns = not significant

## Summary statistics

- **Total independent analyses:** 9
- **Collapse (artifact detected):** 5 (56%)
- **Survive (signal validated):** 3 (33%)
- **Partial (component-specific):** 1 (11%)

## Key patterns

1. **Risk level predicts outcome.** All analyses using "Very high" or "High" risk data showed at least partial collapse. The two analyses that fully survived used either "Moderate" risk (aggregated temporal) or "Low" risk (spatial occurrence) data.

2. **Type A diagnostics are most decisive.** In all three Type A applications, the diagnostic produced a clear verdict (ratio > 1.0 in all cases). Type B and C diagnostics produced more nuanced results, distinguishing components that survived from those that did not.

3. **Structured protocols reduce but do not eliminate bias.** The eBird analysis (row 5b) showed a smaller absolute diagnostic effect than the all-GBIF analysis, but the ratio remained well above 1.0, indicating that bias persists even in semi-structured platforms.

4. **Diagnostic ratios cluster bimodally.** Genuine signals produced ratios < 0.5 (after correction), while artifacts produced ratios > 1.0. The gap between 0.5 and 1.0 was sparsely populated in our empirical analyses, consistent with the simulation finding that the 0.5 threshold provides reasonable discrimination.
