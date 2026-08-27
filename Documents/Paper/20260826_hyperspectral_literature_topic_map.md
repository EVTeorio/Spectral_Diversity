# Hyperspectral Literature Topic Map For Introduction Development

Created: 2026-08-26

Purpose: organize the papers in the Zotero Hyperspectral Paper folder around the topics that support the updated research direction: scale-dependent spectral-biodiversity relationships, PCA-based spectral heterogeneity, phylogenetic signal in spectra, and differences among biodiversity metrics.

## 1. Spectral Variation Hypothesis And Forest Biodiversity Across Spatial Scales

### Crofts et al. 2024 - Linking aerial hyperspectral data to canopy tree biodiversity

Relevance: One of the closest conceptual matches. The paper directly tests the spectral variation hypothesis with airborne imaging spectroscopy and forest biodiversity measurements, and reports that spectral-biodiversity relationships depend on the biodiversity dimension and the spectral method used.

Difference from this study: Their framework compares broader biodiversity dimensions and environmental controls, while this study focuses more narrowly on PCA-derived canopy spectral heterogeneity metrics across 10, 20, and 50 m quadrat scales, with emphasis on vector-normalized PCA and phylogenetic diversity.

Source: https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1605

### Hayden et al. - Scale dependence in remotely sensed biodiversity

Relevance: Directly supports the scale-centered framing. The paper uses NEON imaging spectroscopy to evaluate how spectral diversity changes across spatial scales and shows that scale can alter inferred biodiversity patterns.

Difference from this study: Hayden et al. use continental-scale NEON data, plot-specific PCA workflows, and broad scaling relationships. This study uses a single forest system and compares biodiversity correlations across fixed quadrat grains using a common PCA basis and vector-normalized PCA metrics.

Source: https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.70068

### Crofts et al. 2026 - Testing the scale dependence of plant community assembly processes using imaging spectroscopy

Relevance: Supports the argument that scale changes ecological interpretation. The study links imaging spectroscopy to community assembly and shows that increasing grain can alter the inferred balance of ecological processes.

Difference from this study: Crofts et al. 2026 focus on assembly processes, while this paper focuses on how spectral heterogeneity metrics correlate with observed biodiversity metrics, especially phylogenetic diversity.

Sources: https://www.researchgate.net/publication/403374956_Testing_the_scale_dependence_of_plant_community_assembly_processes_using_imaging_spectroscopy ; https://visualize.jove.com/41916596-Testing-the-scale-dependence-of-plant-community-assembly-processes-using-imaging-spectroscopy

### Donnini, Kross, and Alejo 2024 - Spectral Diversity as a Predictor of Tree Diversity

Relevance: Connects spectral diversity to forest tree diversity and explicitly examines challenges and opportunities across forest ecosystems.

Difference from this study: Donnini et al. use Sentinel-2 scale observations across broader forest types. This study uses finer spatial hyperspectral information and asks which PCA-based spectral metric and quadrat scale best tracks phylogenetic diversity.

Sources: https://doi.org/10.1080/07038992.2024.2403495 ; https://doi.org/10.5683/SP3/U8TI1E

### Rocchini et al. 2010 - Remotely sensed spectral heterogeneity as a proxy of species diversity

Relevance: Foundational review for the spectral variation hypothesis. It highlights scale matching, spectral heterogeneity measurement, sensor choice, taxonomic diversity measures, and modeling challenges.

Difference from this study: Rocchini et al. synthesize the early SVH literature, while this study contributes a focused empirical comparison among PCA-based metrics, vector normalization, and biodiversity transformations.

Source: https://cris.unibo.it/handle/11585/715338

### Torresani et al. 2024 - Reviewing the Spectral Variation Hypothesis

Relevance: Current review that situates this study within two decades of SVH work, including metric choice, ecosystem context, sensor characteristics, and open limitations.

Difference from this study: The review summarizes the field broadly; this study provides a local, metric-explicit test of PCA-derived canopy spectral heterogeneity against biodiversity and phylogenetic metrics.

Source: https://www.sciencedirect.com/science/article/pii/S1574954124002449

## 2. PCA-Based, Distance-Based, And Geometry-Based Spectral Measures

### Hayden et al. - Scale dependence in remotely sensed biodiversity

Relevance: Uses PCA-based spectral diversity approaches and tests sensitivity to the number of retained PCs. This is important background for using reduced spectral feature spaces instead of all original bands.

Difference from this study: Their PCA is plot-specific and continental in scope. This study uses global PCA scores built from the sampled quadrat pixels, then compares raw and vector-normalized PCA metrics across quadrat scales.

Source: https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.70068

### Rocchini, Marcantonio, and Ricotta 2017 - Measuring Rao's Q diversity index from remote sensing

Relevance: Provides the conceptual basis for using Rao's Q with remote sensing data: spectral diversity can incorporate distances among pixel values rather than only frequency-like abundance.

Difference from this study: Rao's Q is one of the spectral metrics compared here, but this study shows why Rao's Q can behave differently from standardized PCA mean distance and alpha hull area because it weights pairwise pixel distances rather than centroid dispersion or occupied spectral space.

Source: https://www.sciencedirect.com/science/article/pii/S1470160X16304319

### Messinger, Ziemann, and Schlamm 2010 - Spectral image complexity estimated through local convex hull volume

Relevance: Directly relevant to convex hull spectral volume as a geometry-based measure of spectral complexity.

Difference from this study: Convex hull style metrics are considered here but are not carried forward for future analysis because they are highly sensitive to boundary points and extreme pixels compared with the main PCA metrics retained for interpretation.

### Schneider et al. 2017 - Mapping functional diversity from remotely sensed forest traits

Relevance: Shows how remotely sensed multivariate canopy information can be translated into functional diversity metrics and evaluated across spatial scale.

Difference from this study: Schneider et al. map trait-based functional richness, divergence, and evenness from modeled forest traits. This study uses spectral PCA space directly and asks how spectral geometry relates to taxonomic, abundance, and phylogenetic biodiversity metrics.

Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC5682291/

## 3. Spectral Shape, Plant Function, And Phylogenetic Signal

### Cavender-Bares et al. 2016 - Associations of leaf spectra with genetic and phylogenetic variation in oaks

Relevance: Supports the biological mechanism that spectra can contain phylogenetic information because leaf optical properties reflect coordinated chemistry, structure, and evolutionary history.

Difference from this study: Cavender-Bares et al. work at the leaf and oak lineage scale using 400-2400 nm spectra. This study tests whether canopy-level spectral heterogeneity, after image preprocessing and vector normalization, tracks plot-level phylogenetic diversity.

Source: https://www.mdpi.com/2072-4292/8/3/221

### Meireles, O'Meara, and Cavender-Bares 2020 - Linking Leaf Spectra to the Plant Tree of Life

Relevance: Useful for framing why spectra may relate to phylogeny but should be interpreted as phenotype rather than as phylogeny itself. A key idea is that spectra summarize evolved leaf traits, not evolutionary history directly.

Difference from this study: The chapter focuses on leaf spectra and macroevolutionary interpretation. This study evaluates whether canopy spectral heterogeneity can serve as a remote proxy for phylogenetic diversity at quadrat scales.

Source: https://link.springer.com/chapter/10.1007/978-3-030-33157-3_7

### Meireles et al. 2020 - Leaf reflectance spectra capture the evolutionary history of seed plants

Relevance: Strong support for phylogenetic signal in reflectance spectra across broad plant lineages.

Difference from this study: The source uses leaf-level spectra from many seed plants. This study works with canopy pixels, where shadow, illumination, crown mixing, and spatial grain complicate the biological signal.

Source: https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.16771

### Blanchard, Bruneau, and Laliberte 2024 - Foliar spectra distinguish temperate tree species and show phylogenetic signal

Relevance: Especially relevant because it focuses on temperate tree species. It supports the expectation that close relatives can be spectrally similar and more distant relatives can be more separable.

Difference from this study: Their evidence comes from foliar spectra and species discrimination. This study tests whether those leaf-level and species-level differences are strong enough to appear in canopy-scale spectral heterogeneity and correlate with quadrat-level phylogenetic diversity.

Source: https://www.citedrive.com/en/discovery/foliar-spectra-accurately-distinguish-most-temperate-tree-species-and-show-strong-phylogenetic-signal/

### Villa et al. 2025 - Spectral and phylogenetic diversity links with functional structure

Relevance: Supports integrating spectral and phylogenetic information to understand biodiversity and functional diversity.

Difference from this study: Villa et al. examine aquatic plant communities and functional structure. This study focuses on terrestrial forest canopy hyperspectral data and abundance-weighted phylogenetic diversity.

Source: https://iris.cnr.it/handle/20.500.14243/524865

## 4. Environmental, Phenological, And Canopy Drivers Of Spectral Heterogeneity

### Thornley et al. 2022 - Taxonomic and phenological drivers of spectral variance in grasslands

Relevance: Shows that spectral variance is not purely a species-diversity signal. Phenology and site context can strongly shape spectral variance.

Difference from this study: Thornley et al. focus on grasslands and intra-annual phenology. This study uses forest canopy imagery at a fixed acquisition context, but the same caution applies: illumination, shadow, phenological condition, and canopy structure can all influence spectral heterogeneity.

Source: https://nora.nerc.ac.uk/id/eprint/532304

### Davrinche and Haider 2024 - Soil conditions modify species diversity effects on tree traits

Relevance: Useful for explaining why species diversity does not translate into one fixed spectral response. Environmental conditions can modify trait expression, and spectral reflectance is partly a trait-expression signal.

Difference from this study: Davrinche and Haider study soil and functional trait expression experimentally, not remote sensing. This study observes canopy spectral patterns and treats environmental and structural effects as possible contributors to residual variation.

Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC11269567/

### Root and Nelson 2011 - Phylogenetic distance and environmental gradients

Relevance: Helps frame why phylogenetic diversity can sometimes align with environmental and spectral gradients if relevant traits are phylogenetically conserved.

Difference from this study: Root and Nelson evaluate phylogenetic distance in relation to environmental gradients and species composition, whereas this study evaluates phylogenetic diversity against hyperspectral heterogeneity.

Source: https://www.researchgate.net/publication/230547866_Does_phylogenetic_distance_aid_in_detecting_environmental_gradients_related_to_species_composition

### Schneider et al. 2017 - Mapping functional diversity from remotely sensed forest traits

Relevance: Shows that terrain, soil, and canopy conditions can structure remotely sensed functional diversity patterns.

Difference from this study: Schneider et al. model explicit functional traits; this study interprets PCA-based spectral heterogeneity as an integrated canopy signal without direct trait inversion.

Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC5682291/

## 5. Biodiversity Metrics: Taxonomic, Abundance-Weighted, Functional, And Phylogenetic Dimensions

### Faith 1992 - Conservation evaluation and phylogenetic diversity

Relevance: Foundational source for phylogenetic diversity as accumulated evolutionary history rather than species count alone.

Difference from this study: Faith establishes the conceptual and conservation basis of phylogenetic diversity. This study uses phylogenetic metrics as response variables in a remote-sensing test of spectral heterogeneity.

Sources: https://publications.australian.museum/1992-conservation-evaluation-and-phylogenetic-diversity/ ; https://bishtref.com/articles/10.1016/0006-3207%2892%2991201-3

### Scherson and Faith / Faith 2018 - Phylogenetic Diversity and Conservation Evaluation

Relevance: Useful for discussing multiple values, indices, and spatial scales of phylogenetic diversity.

Difference from this study: The chapter is conceptual and conservation-oriented. This study operationalizes phylogenetic diversity in a remote-sensing correlation framework.

Source: https://ouci.dntb.gov.ua/en/works/lxqwZbE9/

### Crofts et al. 2024 - Linking aerial hyperspectral data to canopy tree biodiversity

Relevance: Important because it compares biodiversity dimensions and shows that spectral relationships may be stronger for some dimensions than others.

Difference from this study: This paper's central contrast is narrower: phylogenetic diversity, Shannon diversity, richness, abundance, and related metrics are compared against three PCA-based spectral heterogeneity measures across scale.

Source: https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1605

### Villa et al. 2025 - Spectral and phylogenetic diversity links with functional structure

Relevance: Demonstrates the value of combining spectral and phylogenetic information when interpreting biodiversity patterns.

Difference from this study: Their endpoint is functional structure in aquatic plants; this paper emphasizes why phylogenetic diversity is more detectable than Shannon-like diversity in forest canopy spectra.

Source: https://iris.cnr.it/handle/20.500.14243/524865

## 6. Spectral Species And Community-Level Spectral Diversity

### Rocchini et al. 2021 - From local spectral species to global spectral communities

Relevance: Useful conceptual background for treating spectral variation as a community-level property that can be summarized across space.

Difference from this study: Rocchini et al. derive spectral species and spectral communities from broader remote-sensing data. This study calculates continuous PCA-based heterogeneity metrics within fixed forest quadrats rather than classifying spectral entities.

Source: https://www.sciencedirect.com/science/article/abs/pii/S157495412030145X

### Clark and Roberts 2012 - Species-level differences in hyperspectral metrics among tropical rainforest trees

Relevance: Supports the premise that tree species can differ spectrally across tissue, pixel, and crown scales.

Difference from this study: Clark and Roberts focus on species classification using hyperspectral metrics and machine learning. This study does not classify species from imagery; it asks whether aggregate canopy spectral heterogeneity reflects plot-level biodiversity.

Source: https://doaj.org/article/026245cb8a5c449289e86339a24d03ee

## 7. Comparison Of This Analytical Framework To The Literature

| Analytical element | This study | Closest literature parallels | Main difference or contribution |
| --- | --- | --- | --- |
| Spatial scale | 10, 20, and 50 m quadrats | Hayden et al.; Crofts et al. 2026; Schneider et al. | Fixed quadrat grains are used to test how correlations change as local canopy variation is aggregated. |
| Spectral feature space | Raw PCA and vector-normalized PCA from canopy pixels | Hayden et al.; Crofts et al. 2024 | Vector normalization is used to reduce brightness dominance and emphasize spectral shape before calculating heterogeneity. |
| Main spectral metrics | Standardized PCA alpha hull area, standardized PCA mean Euclidean distance, standardized PCA Rao's Q | Rocchini et al. 2017; Messinger et al.; Schneider et al.; Hayden et al. | Metrics are compared directly to ask how centroid dispersion, pairwise distance, and occupied spectral area differ in biodiversity sensitivity. |
| Biodiversity dimensions | Richness, Shannon-type diversity, abundance, and phylogenetic metrics | Faith; Crofts et al. 2024; Villa et al.; Schneider et al. | The main interpretive focus is why phylogenetic diversity shows stronger relationships than Shannon-like diversity. |
| Transformations | Biodiversity-squared axis transformation emphasized for reporting | Common biodiversity modeling practice; scale literature | The transformation is interpreted as changing sensitivity to high-diversity quadrats, especially at smaller grain sizes. |
| Metric interpretation | Pixel reflectance values drive spectral metrics, while biodiversity metrics differ by species identity, abundance, and evolutionary history | Rocchini 2010; Torresani 2024; Thornley et al. 2022 | The paper explicitly links metric construction to why correlations differ rather than treating all diversity indices as interchangeable. |
| Methods-informing analyses | Shadow mask AUC and brightness-PC1 correlation checks | Imaging spectroscopy preprocessing literature broadly | Intermediate analyses are used to justify downstream spectral metric reliability and the choice to emphasize vector-normalized PCA. |

## 8. Intro-Relevant Synthesis

The literature supports three points that should shape the introduction. First, the spectral variation hypothesis is well established but remains metric-, scale-, and ecosystem-dependent. Second, plant spectra can contain phylogenetic information because leaf optical properties are tied to conserved and convergent traits related to photosynthesis, structure, water, pigments, and chemistry. Third, translating that biological signal from leaves to canopy pixels requires careful attention to image processing, illumination, shadow, spatial grain, and the mathematical behavior of the spectral heterogeneity metric.

This study fits within that literature by asking a narrower and more operational question: among PCA-based canopy spectral heterogeneity metrics, which ones best track biodiversity across quadrat scales, and why do those relationships differ among biodiversity metrics? The framework is most similar to studies that use imaging spectroscopy to test the SVH across scale, but it differs by emphasizing vector-normalized PCA, standardized PCA metric comparisons, and the specific contrast between phylogenetic diversity and Shannon-like biodiversity metrics.
