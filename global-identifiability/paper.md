# Global Indicator Dependency in Pollen-Based Agriculture Detection: Evidence from 776 Sites Across Five Continents

**Version**: 1.0
**Date**: 2026-03-28
**Target journal**: Quaternary Science Reviews / Global Ecology and Biogeography

---

## Abstract

Pollen-based detection of past agricultural impact depends on the indicator taxa selected, yet the magnitude and generality of this dependency have not been tested at global scale. A companion study demonstrated that detection rates at the same sites vary from 0% to 94% depending on indicator choice across Europe and eastern North America, a phenomenon termed the indicator dependency gap (Gordon et al., 2026a). Here we test whether this dependency is universal by applying a standardised exceedance framework to 776 pollen records from five continents: Europe (n = 331), eastern North America (n = 111), China (n = 45), Latin America (n = 179), and sub-Saharan Africa (n = 110). All records were obtained from the Neotoma Paleoecology Database and analysed with identical threshold criteria (baseline mean + 2 SD). European pastoral indicators (Behre, 1981) detect 87% of European sites but only 16% in China, 16% in eastern North America, 28% in Africa, and 41% in Latin America. Region-specific crop indicators recover additional signal in every non-European continent but at highly variable rates (13--58%). Tree taxa show near-universal exceedance (87--94%) on all five continents, confirming their status as non-diagnostic change detectors. The indicator dependency gap --- the range of detection outcomes produced by different indicator sets applied to the same records --- exceeds 60 percentage points on every continent except Europe. These results demonstrate that indicator dependency is a global structural property of pollen-based agricultural detection, not a peculiarity of any single cross-cultural comparison. We propose that all pollen-based assessments of agricultural impact report indicator dependency gaps alongside detection rates, and that region-specific indicator sets be developed and validated before cross-cultural synthesis is attempted.

**Keywords**: pollen analysis, indicator taxa, agricultural impact, indicator dependency, exceedance threshold, Neotoma, cross-cultural comparison, global paleoecology

---

## 1. Introduction

### 1.1 The indicator dependency problem

The detection of past agricultural impact from pollen records is foundational to our understanding of human-environment interaction over millennia. Continental and global syntheses increasingly rely on pollen-derived metrics to reconstruct the timing and magnitude of agricultural transformation (Mottl et al., 2021; Gordon et al., 2024; Nogué et al., 2025). These reconstructions inform debates about Anthropocene onset, baseline vegetation states, and the deep-time relationship between agricultural societies and biodiversity (Ellis et al., 2021; Stephens et al., 2019).

A companion study (Gordon et al., 2026a) demonstrated that agricultural impact detection is not a neutral measurement but a model-dependent operation. Using 442 pollen records from Europe and eastern North America, that study showed that detection rates at the same sites varied from 0% to 94% depending solely on which indicator taxa were selected. European pastoral indicators (*Plantago lanceolata*, *Rumex*, Cerealia-type) detected 87% of European sites but 0% of eastern North American sites, despite the latter representing an independent centre of plant domestication with a well-documented agricultural tradition (Smith, 2006; Fritz, 2019). Region-specific indicators (Amaranthaceae) recovered 39% of those invisible sites. Tree taxa exceeded baseline at approximately 94% of sites on both continents, demonstrating universal compositional change unrelated to any specific agricultural regime.

Gordon et al. (2026a) formalised this phenomenon through a three-domain framework, decomposing agricultural transformation into structural (Type A: forest clearance), crop-specific (Type B: cultivated taxa), and compositional (Type C: tree reorganisation) domains. They demonstrated that when a transformation operates in a domain orthogonal to the one monitored by an indicator set, the transformation is structurally unidentifiable: detection probability is invariant to signal strength (dD/ds approximately 0). The study concluded by proposing that indicator dependency is a general property of pollen-based detection, not limited to the Europe-North America comparison.

### 1.2 From proposal to test

That proposal, however, was based on only two continents. The European-ENA comparison, while revealing, could reflect something particular about the ENA case --- the low pollen productivity of EAC cultigens, the family-level taxonomic resolution of Amaranthaceae, or the specific ecological context of eastern deciduous forests. To claim that indicator dependency is a universal structural property of pollen-based agricultural detection, the framework must be tested against the full diversity of agricultural traditions and ecological settings.

This paper provides that test. We extend the exceedance framework to five continents, adding China (n = 45), Latin America (n = 179), and sub-Saharan Africa (n = 110) to the original European and eastern North American datasets. These three additions were chosen because each represents a distinct agricultural tradition with different crop taxa, different ecological settings, and different relationships between cultivation and forest structure:

**China.** Rice (*Oryza sativa*), millet (*Setaria italica*, *Panicum miliaceum*), and hemp (*Cannabis sativa*) were independently domesticated in the Yangtze and Yellow River valleys (Fuller et al., 2009; Zuo et al., 2017). Rice pollen is morphologically similar to wild Poaceae, creating a taxonomic resolution problem analogous to --- but distinct from --- the Amaranthaceae problem in ENA.

**Latin America.** Maize (*Zea mays*), squash (*Cucurbita*), manioc (*Manihot esculenta*), and a suite of Andean crops represent multiple independent domestication events across Mesoamerica and South America (Piperno, 2011). *Zea mays* is distinctive: wind-pollinated with large, morphologically unambiguous grains, it represents the best-case scenario for crop-specific detection. If indicator dependency persists even when the crop produces abundant, diagnostic pollen, the phenomenon cannot be attributed solely to poor pollen representation.

**Sub-Saharan Africa.** Sorghum (*Sorghum bicolor*), pearl millet (*Pennisetum glaucum*), yam (*Dioscorea*), and oil palm (*Elaeis guineensis*) represent yet another independent agricultural tradition (Marshall & Hildebrand, 2002; Neumann, 2005). African staple crops pose acute pollen identification challenges: sorghum and millet pollen cannot be reliably distinguished from wild grasses at the family level (Poaceae), and yam produces negligible pollen. Oil palm is the primary morphologically diagnostic crop indicator for the continent.

Together, these five regions encompass the major independent centres of agricultural origin, a range of ecological biomes from temperate deciduous forest to tropical rainforest to savanna, and the full spectrum of crop pollen diagnosticity from unambiguous (*Zea*) to cryptic (sorghum/millet as Poaceae). If indicator dependency is observed across all five, it constitutes strong evidence that the phenomenon is structural rather than contingent.

### 1.3 Study design

We apply three indicator categories to all five continental datasets using identical analytical procedures:

1. **European indicators** (Behre, 1981): *Plantago*, *Rumex*, Cerealia-type --- the default framework in global syntheses.
2. **Region-specific indicators**: taxa diagnostic of each continent's agricultural tradition (Amaranthaceae for ENA, Oryza-type/Cannabis for China, Zea/Cucurbita/Manihot for Latin America, Elaeis/Poaceae for Africa).
3. **Tree taxa**: continent-appropriate arboreal genera as universal change detectors.

For each indicator category and continent, we calculate the exceedance rate (proportion of sites where any taxon exceeds its pre-agricultural baseline) and the indicator dependency gap (difference between maximum and minimum detection rates). The five-continent comparison table is the paper's central result.

---

## 2. Methods

### 2.1 Pollen data

All pollen records were obtained from the Neotoma Paleoecology Database (Williams et al., 2018) using the `neotoma2` R package. Sites were selected to maximise geographic coverage within each region while requiring a minimum of 10 Holocene samples and chronological control from at least two radiocarbon dates. The five regional datasets comprise:

- **Europe** (n = 331): Central Europe, British Isles, and Scandinavia. Agricultural baseline set at >5,500 cal BP, corresponding to the pre-Neolithic period in most of northwestern Europe.
- **Eastern North America** (n = 111): 30--50 degrees N, 65--95 degrees W. Agricultural baseline set at >3,000 cal BP, corresponding to the pre-EAC period.
- **China** (n = 45): East and Central Asia (25--56 degrees N, 75--135 degrees E), including sites in Japan and Korea where comparable agricultural traditions existed. Agricultural baseline set at >5,000 cal BP, pre-dating rice intensification in most regions.
- **Latin America** (n = 179): Central and South America (55 degrees S to 25 degrees N, 120 degrees W to 35 degrees W), spanning Mesoamerica, northern South America, central South America, and southern South America. Agricultural baseline set at >5,000 cal BP.
- **Sub-Saharan Africa** (n = 110): 35 degrees S to 15 degrees N, 20 degrees W to 55 degrees E, spanning West Africa, Central Africa, East Africa, and southern Africa. Agricultural baseline set at >5,000 cal BP.

The combined dataset of 776 sites represents one of the largest multi-indicator, multi-continental pollen analyses for agricultural impact detection.

### 2.2 Exceedance framework

The exceedance framework follows Gordon et al. (2026a). For each taxon at each site:

1. A baseline period was defined from all samples predating the regional agricultural chronology.
2. The exceedance threshold was set at the baseline mean plus 2 standard deviations (97.5th percentile).
3. A taxon was classified as exceeding baseline if any post-baseline sample surpassed this threshold.
4. For each indicator set, a site was scored as "exceedance" if any taxon in the set exceeded its baseline.

Detection rates are reported both as proportions of testable sites (sites where at least one taxon in the set was recorded) and as proportions of total sites. The latter is the more conservative metric and is used for cross-continental comparison because it accounts for the possibility that absence of a taxon from the record reflects genuine absence from the landscape rather than sampling artefact.

### 2.3 Indicator sets

**European pastoral indicators (Set A).** *Plantago lanceolata*, *Plantago* spp., *Rumex*, Cerealia-type, *Secale*. Applied identically across all five continents. This is the indicator set most commonly used in global syntheses and serves as the test case for framework transferability.

**Region-specific indicators (Set B).** Defined separately for each non-European continent based on the archaeobotanical literature:

- *ENA*: Amaranthaceae (including Chenopodiaceae), targeting EAC cultigens (*Chenopodium berlandieri*, *Iva annua*).
- *China*: Oryza-type, *Cannabis*/*Humulus*, *Fagopyrum* --- markers of East Asian rice, hemp, and buckwheat agriculture.
- *Latin America*: *Zea mays*, *Cucurbita*, *Manihot*, *Ipomoea*, *Phaseolus*, *Capsicum*, *Chenopodium*, *Amaranthus* --- the Mesoamerican and Andean crop suite.
- *Africa*: *Elaeis guineensis*, Poaceae (including size-sorted categories as sorghum/millet proxies), *Dioscorea*, *Vigna*, *Musanga* (secondary forest indicator of forest disturbance).

**Tree taxa (Set E).** Continent-appropriate arboreal genera selected to capture the dominant forest composition:

- *Europe*: *Quercus*, *Fagus*, *Betula*, *Pinus*, *Corylus*, *Ulmus*, *Alnus*, *Tilia*, *Fraxinus*, *Picea*, *Carpinus*, *Acer*.
- *ENA*: *Quercus*, *Pinus*, *Betula*, *Acer*, *Carya*, *Picea*, *Fagus*, *Tsuga*, *Castanea*.
- *China*: *Quercus*, *Pinus*, *Betula*, *Picea*, *Abies*, *Ulmus*, *Fagus*, *Tsuga*, *Castanopsis*, *Podocarpus*, *Cryptomeria*, *Liquidambar*.
- *Latin America*: *Nothofagus*, *Podocarpus*, *Araucaria*, *Quercus*, *Alnus*, *Hedyosmum*, *Weinmannia*, *Polylepis*, *Cecropia*, *Moraceae*, *Myrtaceae*, *Melastomataceae*, and additional tropical genera.
- *Africa*: *Olea*, *Podocarpus*, *Juniperus*, *Hagenia*, *Macaranga*, *Celtis*, *Combretaceae*, *Uapaca*, *Brachystegia*, *Syzygium*, and additional forest and savanna genera.

### 2.4 Onset timing

For each exceedance event, the onset age was recorded as the calibrated age of the first post-baseline sample exceeding the threshold. Median and mean onset ages are reported per indicator set per continent to assess whether exceedance timing is consistent with known agricultural chronologies.

### 2.5 Regional sub-analyses

Latin America and Africa were subdivided into ecological-cultural sub-regions to test for within-continent heterogeneity:

- *Latin America*: Mesoamerica (n = 20), northern South America (n = 47), central South America (n = 42), southern South America (n = 70).
- *Africa*: West Africa (n = 9), Central Africa (n = 47), East Africa (n = 14), southern Africa (n = 40).

### 2.6 Statistical framework

Confidence intervals use the Clopper-Pearson exact method. The indicator dependency gap is defined as the difference between the maximum and minimum detection rates (total-site basis) across indicator sets within a single region. All analyses were conducted in R 4.3.x.

---

## 3. Results

### 3.1 Europe: the calibration continent

The European results from Gordon et al. (2026a) serve as the calibration benchmark. European pastoral indicators (Set A) detect exceedance at 87% of 331 sites (95% CI [83%, 90%]). The signal is dominated by pastoral taxa (*Plantago lanceolata*, *Rumex*, Poaceae), which exceed pre-5,500 BP baselines well before cereal pollen in most records. Tree taxa (Set E) show 93% exceedance, driven by late-Holocene genus migrations (*Picea* expansion, *Fagus* migration) and anthropogenic deforestation (*Quercus* decline).

In Europe, the indicator dependency gap is narrow (6 percentage points) because the dominant indicator framework was designed for and calibrated against exactly this agricultural tradition. Europe is the continent where framework and reality are aligned. The question is what happens when this alignment is broken.

### 3.2 Eastern North America: the indicator dependency gap discovered

Eastern North America provides the extreme case. European pastoral indicators detect 0% of 111 sites (Gordon et al., 2026a) --- later revised to 16% with expanded indicator definitions including broader *Plantago* and *Rumex* species concepts. Region-specific indicators (Amaranthaceae) detect 39% of testable sites. Tree taxa detect 94%.

The ENA indicator dependency gap is 78 percentage points (94% tree minus 16% European). The gap is driven by the fundamental mismatch between European pastoral agriculture (structural transformation) and the Eastern Agricultural Complex (crop-specific transformation within an existing forest matrix). This result motivated the present global extension.

### 3.3 China: rice agriculture and the Poaceae problem

China provides the most data-constrained test. Of 45 analysed sites, only 12 had sufficient European indicator taxa and only 8 had Chinese agricultural indicator taxa for exceedance testing. This reflects two factors: the geographic bias of the Chinese Neotoma holdings toward Siberia, Mongolia, and the Tibetan Plateau (rather than the rice-growing Yangtze and Yellow River valleys), and the low taxonomic resolution of Poaceae-family crops.

**European indicators.** Seven of 12 testable sites (58% of testable, 16% of total) show exceedance. Six of the seven are located in Central Asia (Kazakhstan, western Siberia), where *Plantago* and *Rumex* are native components of steppe vegetation. The Zoige Basin site on the Tibetan Plateau margin (4,558 BP onset) is the only exceedance in core China. European indicators detect little meaningful agricultural signal in the Chinese dataset.

**Region-specific indicators.** Six of 8 testable sites (75% of testable, 13% of total) show exceedance. The Chinese agricultural indicator set (Oryza-type, *Cannabis*, *Fagopyrum*) recovers signal predominantly from *Cannabis*/*Humulus* (found at Kanas Lake, Lake Baikal, Wenquan wetland) and Oryza-type (found at Lake Tianchi, Yunnan). The absence of Oryza-type pollen from most records reflects both the geographic bias of available sites (away from the Yangtze delta where rice was domesticated) and the difficulty of distinguishing Oryza-type pollen from wild Poaceae at standard taxonomic resolution.

**Tree taxa.** Forty-one of 42 testable sites (98% of testable, 91% of total) show exceedance. Tree exceedance in China is driven by Holocene forest dynamics including *Betula* and *Pinus* fluctuations across the forest-steppe ecotone and *Quercus* compositional shifts in the temperate mixed forest zone.

**Indicator dependency gap.** The gap is 78 percentage points (91% tree minus 13% Chinese agricultural). This is virtually identical to the ENA gap (78 pp), despite the entirely different agricultural tradition, ecological setting, and pollen taxonomic challenges involved.

**Onset timing.** Chinese agricultural exceedances show median onset at 4,421 BP (mean 4,064 BP), consistent with the mid-Holocene intensification of rice and hemp cultivation. European indicator exceedances show later median onset (3,276 BP), consistent with Central Asian agropastoral expansion rather than indigenous Chinese agriculture.

### 3.4 Latin America: Zea mays and the cosmopolitan indicator problem

Latin America provides the largest non-European dataset (179 sites) and the most internally heterogeneous agricultural landscape.

**European indicators.** Seventy-four of 105 testable sites (71% of testable, 41% of total) show exceedance --- the highest non-European rate. This requires explanation. *Plantago* and *Rumex* are cosmopolitan genera with native species in the Americas: *Plantago rigida*, *P. australis*, *P. sericea*, and *P. barbata* are native Andean and Patagonian species; *Rumex acetosella* was likely introduced but *Rumex* spp. are present in pre-Columbian records. The 41% European indicator exceedance rate in Latin America thus conflates two signals: (a) genuine agricultural disturbance detected through cosmopolitan weedy taxa and (b) natural exceedance of native *Plantago* and *Rumex* species in Andean grasslands and Patagonian steppe. This cosmopolitan false-positive problem does not arise in ENA or China, where *Plantago lanceolata* is unambiguously introduced.

**Region-specific indicators.** Forty-eight of 55 testable sites (87% of testable, 27% of total) show exceedance. The Mesoamerican indicator set (*Zea mays*, *Cucurbita*, *Manihot*, *Ipomoea*, *Phaseolus*, *Capsicum*, *Chenopodium*, *Amaranthus*) detects agricultural impact most effectively in its region of origin (Mesoamerica: 55% of all sites) and in northern South America (43%), declining sharply to 6% in southern South America.

**The Zea gradient.** *Zea mays* pollen provides the sharpest test of geographic specificity. It is testable at 30 sites (present in at least one sample) and shows 100% exceedance at all testable sites. But testability itself is geographically structured: Zea is present at 50% of Mesoamerican sites, 26% of northern South American sites, 17% of central South American sites, and only 1% of southern South American sites. The Zea exceedance gradient --- from 50% of all Mesoamerican sites to 1% of all southern South American sites --- traces the geographic attenuation of maize agriculture from its Mesoamerican centre of origin into southern temperate latitudes where it was never a staple. This gradient is precisely what the identifiability framework predicts: where the crop was not cultivated, its pollen is absent, and detection under a crop-specific indicator is zero regardless of other agricultural activities (such as *Araucaria* management or camelid pastoralism) that may have been occurring.

**Tree taxa.** One hundred fifty-seven of 170 testable sites (92% of testable, 88% of total) show exceedance. Tree exceedance in Latin America is driven by tropical forest dynamics (*Cecropia* as secondary forest indicator, *Moraceae*/*Urticaceae* fluctuations), Andean treeline shifts (*Polylepis*, *Weinmannia*), and southern *Nothofagus* forest dynamics.

**Indicator dependency gap.** The gap is 61 percentage points (88% tree minus 27% Mesoamerican). If the inflated European indicator rate (41%) is used as the maximum herb indicator, the gap narrows to 47 pp --- still substantial, though complicated by the cosmopolitan false-positive issue.

**Regional breakdown.** Within-continent heterogeneity is pronounced:

| Sub-region | n | European | Mesoamerican | Zea only | Tree |
|---|---|---|---|---|---|
| Mesoamerica | 20 | 10% | 55% | 50% | 95% |
| Northern S. America | 47 | 53% | 43% | 26% | 96% |
| Central S. America | 42 | 38% | 31% | 17% | 88% |
| Southern S. America | 70 | 44% | 6% | 1% | 80% |

Mesoamerica shows the inverted dependency pattern: European indicators detect only 10% (because *Plantago lanceolata* is absent), while region-specific indicators detect 55%. Northern South America shows the highest European rate (53%) due to cosmopolitan *Plantago*/*Rumex* species in Andean environments. Southern South America shows minimal Mesoamerican detection (6%) because maize agriculture did not penetrate significantly south of 30 degrees S.

### 3.5 Sub-Saharan Africa: the Poaceae ceiling and the West African gap

Africa provides the most challenging pollen-taxonomic context for agricultural detection.

**European indicators.** Thirty-one of 45 testable sites (69% of testable, 28% of total) show exceedance. As in Latin America, *Plantago* has native African species (*P. africana*), and *Rumex* is present in highland East African vegetation. However, the geographic distribution of European indicator exceedance is informative: West Africa shows 0% (0/9 sites), while East Africa shows 50% (7/14) and Central Africa shows 36% (17/47). The complete absence of European indicator exceedance in West Africa --- the heartland of indigenous African agriculture --- is the purest example of the identifiability gap in our global dataset.

**Region-specific indicators.** Sixty-four of 99 testable sites (65% of testable, 58% of total) show exceedance --- the highest region-specific rate of any non-European continent. The African indicator set is dominated by Poaceae (wild grass family, which subsumes sorghum and millet pollen) and *Elaeis guineensis* (oil palm). The 58% rate must be interpreted with caution because Poaceae exceedance conflates agricultural (sorghum/millet cultivation) and non-agricultural (fire-driven savanna expansion, climate-driven grassland shifts) signals. This is the African analogue of the Oryza-type/wild Poaceae problem in China.

**The Elaeis test.** *Elaeis guineensis* provides the most diagnostic crop-specific indicator for Africa. It is testable at 23 sites (present in at least one sample) and shows exceedance at 15 (65% of testable, 14% of total). Elaeis exceedance is concentrated in West and Central Africa (West Africa: 67% of sites where testable, Central Africa: 19%), exactly corresponding to the oil palm cultivation zone. The 14% total-site rate for the most diagnostic African crop indicator illustrates the severity of the detection problem: even with the best available marker, over 85% of sites show no crop-specific exceedance.

**Tree taxa.** Ninety-seven of 102 testable sites (95% of testable, 88% of total) show exceedance. Tree exceedance in Africa is driven by forest-savanna boundary dynamics (Macaranga, Celtis, Moraceae fluctuations), highland forest shifts (Podocarpus, Juniperus, Hagenia in East Africa), and miombo woodland dynamics (Brachystegia, Julbernardia, Uapaca in southern Africa).

**Indicator dependency gap.** The gap is 60 percentage points (88% tree minus 28% European). If the Poaceae-inflated African indicator rate (58%) is used as maximum, the gap narrows to 30 pp, but this comparison is misleading because the African indicator rate is artificially elevated by non-diagnostic Poaceae.

**Regional breakdown:**

| Sub-region | n | European | African | Elaeis only | Tree |
|---|---|---|---|---|---|
| West Africa | 9 | 0% | 78% | 67% | 89% |
| Central Africa | 47 | 36% | 60% | 19% | 92% |
| East Africa | 14 | 50% | 21% | 0% | 86% |
| Southern Africa | 40 | 18% | 65% | 0% | 85% |

The West Africa result is the paper's most striking regional finding. European indicators detect 0% of sites in the region with the longest indigenous agricultural history in Africa. African indicators detect 78%, driven by oil palm and Poaceae. This 78-percentage-point within-region gap (0% European vs. 78% African) equals the most extreme gaps found anywhere in the global dataset and demonstrates that the indicator dependency problem is not merely a between-continent phenomenon --- it operates at sub-continental scales wherever agricultural traditions diverge from the European model.

### 3.6 Five-continent synthesis

**Table 1.** Global indicator exceedance comparison (776 sites, 5 continents). This is the paper's central result.

| Indicator set | Europe (331) | ENA (111) | China (45) | Latin America (179) | Africa (110) |
|---|---|---|---|---|---|
| European (Behre) | 87% | 16% | 16% | 41%* | 28% |
| Region-specific | --- | 39% (EAC) | 13% (rice) | 27% (Meso) | 58% (palm+Poac)** |
| Tree taxa | 93% | 94% | 91% | 88% | 88% |
| Dependency gap (pp) | 6 | 78 | 78 | 61 | 60 |

\* Inflated by cosmopolitan native *Plantago*/*Rumex* species.
\** Inflated by non-diagnostic Poaceae (sorghum/millet indistinguishable from wild grasses).

**Table 2.** Onset timing by indicator set (median, cal BP).

| Indicator set | Europe | ENA | China | Latin America | Africa |
|---|---|---|---|---|---|
| European | ~5,200 | post-contact | 3,276 | 3,504 | 3,903 |
| Region-specific | --- | ~2,700 | 4,421 | 3,298 | 3,573 |
| Tree | ~5,000 | ~5,000 | 4,569 | 4,604 | 4,680 |

Three patterns are universal across all five continents:

**Pattern 1: Tree exceedance universality.** Tree taxa exceed baseline at 87--94% of sites on every continent. The range is narrow (7 percentage points) despite enormous differences in forest type, climate, and human history. Tree compositional change is a structural feature of Holocene forests worldwide, driven by climate-vegetation feedbacks, species migrations, fire regime changes, and anthropogenic disturbance in varying proportions. It detects everything and therefore diagnoses nothing specific.

**Pattern 2: European indicator attenuation.** The European pastoral indicator set shows monotonic attenuation away from its calibration region: 87% (Europe) > 41% (Latin America, inflated) > 28% (Africa) > 16% (China) approximately equals 16% (ENA). The non-European rates reflect (a) presence of cosmopolitan congeners (*Plantago*, *Rumex* in LatAm and Africa), (b) late introduction of European weeds via trade or colonialism, and (c) coincidental presence of these taxa in non-agricultural contexts. In no case does the European indicator set function as a reliable detector of indigenous non-European agriculture.

**Pattern 3: Region-specific indicators always add detection but never match the European rate in Europe.** Every non-European continent shows some positive detection when region-specific indicators are applied, but the rates range from 13% (China, rice) to 58% (Africa, palm+Poaceae). The variation reflects the pollen diagnosticity of the regional crop suite: where the signature crop produces distinctive, abundant pollen (*Zea mays* in Mesoamerica: 50% detection), detection is higher; where the crop pollen is taxonomically cryptic (Oryza in China: 13%), detection is lower.

---

## 4. Discussion

### 4.1 Tree exceedance is universal and non-diagnostic

The most robust finding is the convergence of tree exceedance rates across five continents: 93% (Europe), 94% (ENA), 91% (China), 88% (Latin America), 88% (Africa). This near-universality confirms the interpretation advanced by Gordon et al. (2026a): tree compositional change is a baseline property of Holocene forest ecosystems, not a diagnostic marker of agricultural impact. An investigator who relies solely on tree-based metrics (AP decline, forest compositional turnover) will detect "impact" at approximately 90% of sites regardless of the agricultural history of the region.

This has immediate practical consequences. Global compilations that use AP decline or tree compositional turnover as their primary impact metric (e.g., Mottl et al., 2021) are measuring something real --- Holocene vegetation change --- but cannot attribute it to agriculture without independent corroboration. The 88--94% rate sets the baseline false-positive ceiling for tree-based agricultural detection.

### 4.2 European indicators: false negatives and false positives

The European indicator set was designed to detect European pastoral-arable agriculture and does so effectively (87%). When transferred to other continents, it produces two classes of error.

**False negatives** occur where agricultural traditions do not produce the structural transformation (forest-to-open conversion) that European indicators detect. The clearest cases are ENA (16%, driven by post-contact European weed introduction rather than indigenous agriculture), China (16%, driven by Central Asian steppe vegetation rather than rice cultivation), and West Africa (0%, the purest case of complete detection failure in our dataset). These false negatives are not failures of sensitivity --- the exceedance framework detects European pastoral indicators at 87% in their calibration region. They are failures of specificity to the transformation domain: the indicator set monitors the wrong ecological dimension.

**False positives** (or more precisely, detection of non-target processes) occur where cosmopolitan congeners of European indicator taxa are present for reasons unrelated to European-style agriculture. The 41% rate in Latin America is partly or largely attributable to native *Plantago* (*P. rigida*, *P. australis*) and *Rumex* species in Andean and Patagonian environments. The 28% rate in Africa includes native *Plantago africana* in highland East Africa. These detections are not necessarily "false" in the strict sense --- the taxa may increase under agricultural disturbance --- but they conflate European-framework signal with indigenous ecological processes, making the detection rate uninterpretable without additional context.

The false-positive problem is asymmetric across continents. In ENA, where *Plantago lanceolata* is unambiguously introduced, European indicator exceedance can be cleanly attributed to post-contact ecological disruption. In Latin America and Africa, where native congeners exist, the signal is irreducibly ambiguous.

### 4.3 Region-specific indicators: variable but essential

Region-specific indicator sets recover agricultural signal that European indicators miss on every non-European continent. But the recovered fraction varies enormously:

- Africa: 58% (palm + Poaceae) --- but Poaceae is non-diagnostic
- ENA: 39% (Amaranthaceae) --- family-level, with zero false positives
- Latin America: 27% (Mesoamerican crop suite) --- geographically structured
- China: 13% (rice/hemp/buckwheat) --- data gap-limited

The variation is not random. It reflects three independent factors:

**Pollen diagnosticity.** Where the crop pollen is morphologically distinctive and abundantly produced, detection is higher. *Zea mays* (large, diagnostic grains, wind-pollinated) achieves 100% exceedance where testable. *Elaeis guineensis* (distinctive triporate grains) achieves 65% where testable. By contrast, rice (Oryza-type, morphologically similar to wild Poaceae) and sorghum/millet (literally classified as Poaceae) have low effective diagnosticity.

**Geographic coverage of data.** The Chinese dataset (45 sites) is heavily biased toward Siberia, Mongolia, and the Tibetan Plateau rather than the Yangtze and Yellow River valleys where rice agriculture originated. This data gap may be the primary explanation for the low Chinese agricultural detection rate (13%). Sites that are in the core rice zone (e.g., Lake Tianchi, Yunnan) do show Oryza-type exceedance.

**Agricultural intensity and extent.** Southern South America shows only 6% Mesoamerican indicator exceedance because maize agriculture genuinely did not penetrate significantly into southern temperate latitudes. This is a true negative, not a detection failure --- the indicator set correctly reports the absence of a crop that was not cultivated.

### 4.4 The Zea gradient as proof of concept

The *Zea mays* detection gradient across Latin American sub-regions constitutes the strongest evidence that the exceedance framework detects real agricultural geography rather than methodological artefact:

| Sub-region | Zea exceedance (% of all sites) |
|---|---|
| Mesoamerica | 50% |
| Northern South America | 26% |
| Central South America | 17% |
| Southern South America | 1% |

This gradient mirrors the known archaeological distribution of maize agriculture: intensive in Mesoamerica, present but less dominant in northern South America, patchy in central South America, and minimal in southern South America (Piperno, 2011). The gradient is not an artefact of data coverage --- all sub-regions have reasonable site density. It demonstrates that crop-specific indicators, when the crop produces diagnostic pollen, can recover meaningful agricultural geography from pollen records.

The Zea case also highlights why the identifiability problem is structural. Southern South American sites that show no Zea exceedance may still have experienced significant agricultural transformation through other pathways (camelid pastoralism, *Araucaria* nut management, potato cultivation). These transformations are invisible to the Mesoamerican indicator set because they operate in a different ecological domain. A truly comprehensive detection framework for Latin American agriculture would need to include indicators for each of these distinct agricultural traditions --- but for many (potato cultivation, nut management), no palynological marker exists.

### 4.5 West Africa: the purest identifiability gap

West Africa (n = 9 sites) provides the most unambiguous demonstration of the indicator dependency gap in our global dataset:

- European indicators: 0%
- African indicators: 78%
- Tree taxa: 89%

The 0% European indicator rate in the region with some of Africa's oldest agricultural traditions (oil palm management dating to at least 4,000 BP; Neumann, 2005; Sowunmi, 1999) is the purest case of structural unidentifiability. *Plantago lanceolata* is not native to West Africa. *Rumex* is absent from lowland tropical forest pollen records. Cerealia-type pollen is not produced by African cereal crops (sorghum and millet pollen falls within the general Poaceae morphotype). Every assumption embedded in the European framework fails simultaneously.

The 78% African indicator rate demonstrates that agricultural signal is recoverable when the framework matches the agricultural tradition. The signal is dominated by *Elaeis guineensis* (oil palm) and Poaceae exceedance, consistent with the documented expansion of oil palm cultivation and savannisation in the West African forest-savanna mosaic during the late Holocene (Sowunmi, 1999; Maley, 2002).

### 4.6 The African Poaceae ceiling

The high African region-specific rate (58%) must be tempered by a fundamental limitation: Poaceae dominates the African indicator set, and Poaceae exceedance is non-specific. Sorghum (*Sorghum bicolor*), pearl millet (*Pennisetum glaucum*), finger millet (*Eleusine coracana*), and fonio (*Digitaria exilis*) --- the four most important indigenous African cereals --- all produce pollen classified as Poaceae at standard light-microscopy resolution. Their cultivation cannot be distinguished from natural grassland expansion driven by fire, climate change, or herbivore dynamics.

This creates the African analogue of the Chinese Oryza problem, but with greater severity. In China, Oryza-type can at least be tentatively distinguished from wild Poaceae at the sub-family level in some records. In Africa, the cereal pollen is taxonomically invisible. The true diagnostic power of the African indicator set depends almost entirely on *Elaeis guineensis*, which is geographically restricted to the West and Central African forest zone.

If we isolate Elaeis as the sole diagnostic African crop indicator, the detection rate drops to 14% of total sites --- comparable to the Chinese rice rate (13%). This convergence suggests that when crop pollen taxonomic diagnosticity is low, detection rates cluster around 10--15% regardless of agricultural tradition. The higher aggregate African rate (58%) reflects the inclusion of taxonomically non-diagnostic Poaceae, not superior diagnostic power.

### 4.7 The China data gap

The Chinese dataset (45 sites) is the smallest in our analysis and the most geographically biased. The Neotoma Paleoecology Database's Chinese holdings are concentrated in three zones: (1) the Russian Far East (Primorye), (2) southern Siberia and Mongolia, and (3) the Tibetan Plateau margin. The Yangtze River delta --- where rice was domesticated and where the densest early agricultural populations lived --- is represented by zero sites. The Yellow River basin, homeland of millet agriculture, is represented by one marginal site (Zoige Basin, on the Tibetan Plateau edge).

This geographic bias has two consequences. First, it depresses both European and Chinese agricultural detection rates because most sites are in regions where agriculture arrived late (Siberia) or never (Tibetan Plateau). Second, it makes the Chinese results provisional: they demonstrate that indicator dependency exists in the available Chinese data but cannot characterise the full magnitude of the phenomenon in the core agricultural zones. A future analysis incorporating the Chinese Pollen Database (Zheng et al., 2014) and the Global Pollen Database's Chinese holdings would likely show higher region-specific detection rates while maintaining or widening the indicator dependency gap with European indicators.

### 4.8 The indicator dependency gap as a universal property

The most important result of this study is the demonstration that the indicator dependency gap is large (>= 60 pp) on every non-European continent:

| Continent | Gap (pp) | Driving contrast |
|---|---|---|
| Europe | 6 | Framework matches reality |
| ENA | 78 | European 16% vs. Tree 94% |
| China | 78 | Chinese 13% vs. Tree 91% |
| Latin America | 61 | Mesoamerican 27% vs. Tree 88% |
| Africa | 60 | European 28% vs. Tree 88% |

The gap exceeds 60 percentage points everywhere except in Europe, where the framework was developed. This is precisely the pattern predicted by the identifiability framework (Gordon et al., 2026a): when the indicator set's ecological domain aligns with the transformation's ecological domain, detection succeeds and the gap is small. When alignment fails, the gap opens to 60+ percentage points regardless of continent, biome, or agricultural tradition.

The consistency of the gap magnitude across such different settings --- temperate deciduous forest (ENA), continental steppe (China), tropical forest (Africa, Latin America), and Andean grassland (Latin America) --- argues against explanations rooted in any specific ecological or taphonomic context. The indicator dependency gap is a structural property of indicator-based detection systems applied across agricultural traditions, not an empirical accident.

### 4.9 Implications for global synthesis

These results have direct implications for how pollen-based agricultural impact assessments are conducted and reported.

**First**, any cross-cultural comparison of agricultural impact that uses a single indicator set (especially one derived from European traditions) will produce systematically biased results. Regions with European-style agriculture will appear more impacted; regions with other agricultural traditions will appear less impacted or unimpacted. This bias is not correctable post hoc --- it is built into the framework.

**Second**, tree-based metrics (AP decline, compositional turnover) cannot serve as universal indicators of agricultural impact because they detect Holocene forest change at approximately 90% of sites regardless of agricultural history. They are necessary for detecting structural transformation (Type A) but insufficient for attributing change to agriculture.

**Third**, the path forward requires developing, validating, and reporting region-specific indicator sets alongside any European-derived framework. The validation protocol should include: (a) temporal consistency with known agricultural chronologies, (b) geographic consistency with known agricultural distributions (the Zea gradient test), and (c) explicit reporting of the indicator dependency gap.

**Fourth**, the indicator dependency gap should be reported as a standard metric in all pollen-based agricultural impact studies, analogous to how confidence intervals are reported for effect sizes. A study that reports only a single detection rate without acknowledging the range of rates achievable under alternative indicator sets provides an incomplete picture that is structurally prone to interpretive error.

---

## 5. Limitations

**Sample size heterogeneity.** The five continental datasets range from 45 (China) to 331 (Europe) sites. The Chinese and African datasets are particularly small relative to the size and agricultural diversity of their continents. Results for these regions are provisional.

**Geographic bias.** Neotoma's global coverage is uneven. The Chinese dataset lacks core agricultural zones (Yangtze, Yellow River). The African dataset is biased toward the East African highlands and Central African forests, with minimal coverage of the Sahel and Sudan zones where cereal agriculture originated.

**Taxonomic resolution.** The Poaceae problem (sorghum/millet/rice indistinguishable from wild grasses) affects both China and Africa. Our region-specific detection rates for these continents are inflated by non-diagnostic Poaceae exceedance. The true crop-specific detection rate is better represented by the Elaeis-only rate (14%) for Africa and the Oryza-type-only rate (data insufficient) for China.

**Baseline definition.** A uniform >5,000 BP baseline was used for China, Latin America, and Africa, but agricultural chronologies differ substantially across regions. Early oil palm management in West Africa may overlap with our baseline period, producing conservative (under-)estimates of agricultural exceedance.

**Cosmopolitan taxa.** The *Plantago*/*Rumex* false-positive problem in Latin America and Africa cannot be resolved without species-level pollen identification, which is not routinely achieved in these records.

---

## 6. Conclusions

1. **Indicator dependency is global.** The indicator dependency gap exceeds 60 percentage points on every non-European continent tested (ENA: 78 pp; China: 78 pp; Latin America: 61 pp; Africa: 60 pp). This confirms the prediction of Gordon et al. (2026a) that indicator dependency is a structural property of pollen-based detection, not a contingent finding of the Europe-ENA comparison.

2. **Tree exceedance is universal and non-diagnostic.** Tree taxa exceed pre-agricultural baselines at 87--94% of sites on all five continents, confirming that tree compositional change is a structural feature of Holocene forests worldwide, not a marker of any specific agricultural tradition.

3. **European indicators produce systematic false negatives outside Europe.** The European pastoral indicator set detects 87% of sites in its calibration region but only 16--41% elsewhere, with the detection rate reflecting cosmopolitan congener presence rather than agricultural impact. The 0% rate in West Africa represents the purest identifiability failure.

4. **Region-specific indicators are essential but taxonomically constrained.** Every non-European continent shows additional detection when region-specific indicators are applied (13--58%), but the magnitude depends critically on the pollen diagnosticity of the regional crop suite. Where diagnostic markers exist (*Zea mays*, *Elaeis guineensis*), detection is geographically meaningful. Where crop pollen is taxonomically cryptic (Poaceae-family cereals), the detection rate conflates agricultural and non-agricultural signals.

5. **Cross-cultural synthesis requires indicator dependency reporting.** We propose that all pollen-based assessments of agricultural impact report the indicator dependency gap alongside detection rates. Studies that apply a single indicator framework to multiple agricultural traditions without reporting the range of detection outcomes achievable under alternative frameworks provide an incomplete and potentially misleading picture of past human-environment interaction.

---

## Acknowledgements

All pollen data were obtained from the Neotoma Paleoecology Database (www.neotomadb.org). We thank the hundreds of data contributors whose records make continental- and global-scale synthesis possible. The p3k14c radiocarbon database was used for archaeological cross-validation.

---

## References

Behre, K.-E. (1981). The interpretation of anthropogenic indicators in pollen diagrams. *Pollen et Spores*, 23, 225--245.

Bird, D., et al. (2022). p3k14c, a synthetic global database of archaeological radiocarbon dates. *Scientific Data*, 9, 27.

Chevalier, M., et al. (2020). Pollen-based climate reconstruction techniques for late Quaternary studies. *Earth-Science Reviews*, 210, 103384.

Ellis, E. C., et al. (2021). People have shaped most of terrestrial nature for at least 12,000 years. *Proceedings of the National Academy of Sciences*, 118, e2023483118.

Fritz, G. J. (2019). Feeding Cahokia: Early Agriculture in the North American Heartland. University of Alabama Press.

Fuller, D. Q., et al. (2009). The domestication process and domestication rate in rice. *Science*, 323, 1607--1610.

Fyfe, R. M., et al. (2015). The European Pollen Database: past efforts and current activities. *Vegetation History and Archaeobotany*, 24, 243--248.

Gordon, V., et al. (2024). Continental-scale pollen synthesis of agricultural impact. *Quaternary Science Reviews*.

Gordon, V., et al. (2026a). Indicator-dependent detection of agricultural transformation: how pollen frameworks render some agricultural systems visible and others invisible. *Nature Ecology & Evolution* [companion paper].

Maley, J. (2002). A catastrophic destruction of African forests about 2,500 years ago still exerts a major influence on present vegetation formations. *IDS Bulletin*, 33, 13--30.

Marshall, F., & Hildebrand, E. (2002). Cattle before crops: the beginnings of food production in Africa. *Journal of World Prehistory*, 16, 99--143.

Mercuri, A. M., et al. (2013). Pollen and macroremains from Holocene archaeological sites: a dataset for the understanding of the bio-cultural diversity of the Italian landscape. *Review of Palaeobotany and Palynology*, 198, 44--68.

Mottl, O., et al. (2021). Global acceleration in rates of vegetation change over the past 18,000 years. *Science*, 372, 860--864.

Mueller, N. G. (2017). Documenting domestication in a lost crop (*Polygonum erectum* L.): evolutionary bet-hedging. *Vegetation History and Archaeobotany*, 26, 313--327.

Neumann, K. (2005). The romance of farming: plant cultivation and domestication in Africa. In *African Archaeology* (pp. 249--275). Springer.

Nogué, S., et al. (2025). Long-term human impacts on biodiversity. *Nature Ecology & Evolution*.

Piperno, D. R. (2011). The origins of plant cultivation and domestication in the New World tropics. *Current Anthropology*, 52, S453--S470.

Smith, B. D. (2006). Eastern North America as an independent center of plant domestication. *Proceedings of the National Academy of Sciences*, 103, 12223--12228.

Sowunmi, M. A. (1999). The significance of the oil palm (*Elaeis guineensis* Jacq.) in the late Holocene environments of west and west central Africa. *Vegetation History and Archaeobotany*, 8, 199--210.

Stephens, L., et al. (2019). Archaeological assessment reveals Earth's early transformation through land use. *Science*, 365, 897--902.

Williams, J. W., et al. (2018). The Neotoma Paleoecology Database, a multiproxy, international, community-curated data resource. *Quaternary Research*, 89, 156--177.

Zheng, Z., et al. (2014). East Asian pollen database: modern pollen distribution and its quantitative relationship with vegetation and climate. *Journal of Biogeography*, 41, 1819--1832.

Zuo, X., et al. (2017). Dating rice remains through phytolith carbon-14 study reveals domestication at the beginning of the Holocene. *Proceedings of the National Academy of Sciences*, 114, 6486--6491.
