# Final Research Narrative

Created: 2026-08-19

Working title: **Assessing the Spectral Variation Hypothesis Across Spatial Scales Using PCA-based Canopy Spectral Heterogeneity**

## Narrative Pivot

The paper should now be framed as a scale-dependent test of the Spectral Variation Hypothesis using UAV hyperspectral imagery and field-derived biodiversity metrics in the Paint Rock Forest Dynamics Plot.

The previous research track remains important: spectral heterogeneity still appears more strongly related to phylogenetic diversity than to standard species diversity, and abundance-weighted phylogenetic metrics remain central. The revised direction broadens that track by asking why spectral and biodiversity metrics differ from one another, how those differences change with quadrat grain, and which PCA-based spectral heterogeneity metrics best capture biodiversity structure.

The manuscript should no longer read as a single-metric paper about abundance-weighted Faith's PD. It should read as a paper about scale, metric behavior, and correlation structure:

- How do spectral heterogeneity metrics correlate with one another?
- How do biodiversity metrics correlate with one another?
- How do spectral-biodiversity correlations change from 10 m to 20 m to 50 m?
- Why do phylogenetic metrics show weak-to-moderate spectral relationships when Shannon diversity and related species-diversity metrics do not?

## Revised Central Question

How does spatial scale affect relationships among PCA-based canopy spectral heterogeneity, species diversity, and phylogenetic diversity, and which metric combinations provide the strongest evidence for the Spectral Variation Hypothesis?

## Updated Research Questions

1. What drives differences among spectral heterogeneity values?
2. What drives differences among biodiversity metrics?
3. How does quadrat scale affect spectral metrics, biodiversity metrics, and the correlations between them?
4. What quadrat scale produces the strongest phylogenetic diversity versus spectral heterogeneity correlation?
5. Do transformations such as log, square-root, Box-Cox, or power transforms strengthen spectral-biodiversity correlations?
6. What properties allow one spectral heterogeneity metric to capture more biodiversity signal than another?
7. Why does phylogenetic diversity show weak-to-moderate relationships with spectral diversity when Shannon diversity, richness, Simpson diversity, or evenness often do not?

In shorthand, the manuscript should answer:

- How do these metrics correlate?
- Why do they correlate?
- How does scale affect these correlations?

## Objectives To Preserve From The Original Direction

The original project objectives remain valid, but their emphasis should be updated.

1. Quantify spectral heterogeneity from drone-acquired hyperspectral imagery using quadrat-based sampling at 10 m, 20 m, and 50 m.
2. Quantify species and phylogenetic biodiversity from forest census data using canopy-weighted tree-crown overlap.
3. Evaluate relationships between spectral heterogeneity and biodiversity metrics at each spatial scale.
4. Assess whether environmental and illumination-related variables help explain differences among spectral heterogeneity metrics.
5. Evaluate scale dependence as a central result rather than a supporting observation.

## Updated Working Thesis

PCA-based canopy spectral heterogeneity provides scale-dependent support for the Spectral Variation Hypothesis. The strongest spectral-biodiversity relationships occur at broader quadrat grains, especially 50 m, and are most evident for vector-normalized PCA metrics compared with standard species-diversity metrics. This suggests that spatial aggregation reduces fine-scale canopy and illumination noise while allowing spectral summaries to better reflect canopy composition, abundance structure, and phylogenetic dissimilarity.

The paper should remain cautious. These are mostly correlation-based results from one forest plot, 50 m has the fewest quadrats, and spatial autocorrelation has not yet been finalized. The strongest claim is not that hyperspectral imagery directly measures biodiversity, but that PCA-based spectral heterogeneity captures a scale-dependent component of canopy biodiversity that is more visible in phylogenetic metrics than in conventional species-diversity metrics.

## Metric Calculation Summary

This section should be used as a concise methods anchor. It excludes sample-size sensitivity runs and lists the primary calculated spectral and biodiversity metrics, including metrics that may only be supporting or supplemental.

### Spectral Metrics

- `spec_sa`: spectral angle entropy from sunlit, shadow-masked, smoothed 5 nm spectra. Most quadrats use the mean of 70 bootstrap iterations with up to 5,000 retained pixels per iteration; small quadrats use exact all-pixel pairwise spectral angles.
- `spec_pca_mean`: mean Euclidean distance of retained pixels from the quadrat centroid in raw global PCA PC1-PC3 space.
- `spec_pca_med`: median Euclidean distance of retained pixels from the quadrat centroid in raw global PCA PC1-PC3 space.
- `spec_pca_sd`: standard deviation of pixel-to-centroid Euclidean distances in raw global PCA PC1-PC3 space.
- `spec_rao_q`: raw PCA spectral Rao's Q using equal pixel weights and squared Euclidean distance in global PC1-PC3 space; implemented with the equivalent centroid formula for computational efficiency.
- `spec_alpha`: alpha-hull area of retained pixels in raw global PCA PC1-PC2 space.
- `spec_convex`: convex-hull area of retained pixels in raw global PCA PC1-PC2 space. Archived diagnostic only; do not carry forward for future analysis.
- `spec_hull3d_v`: convex-hull volume of retained pixels in raw global PCA PC1-PC3 space. Archived diagnostic only; do not carry forward for future analysis.
- `spec_hull3d_a`: convex-hull surface area of retained pixels in raw global PCA PC1-PC3 space. Archived diagnostic only; do not carry forward for future analysis.
- `spec_spca_mean`: mean Euclidean distance of retained vector-normalized spectra from the quadrat centroid in standardized PCA PC1-PC3 space.
- `spec_spca_med`: median Euclidean distance of retained vector-normalized spectra from the quadrat centroid in standardized PCA PC1-PC3 space.
- `spec_spca_sd`: standard deviation of standardized PCA pixel-to-centroid Euclidean distances.
- `spec_spca_rao`: standardized PCA spectral Rao's Q using equal pixel weights and squared Euclidean distance in standardized PCA PC1-PC3 space; implemented with the equivalent centroid formula.
- `spec_spca_alpha`: alpha-hull area of retained vector-normalized spectra in standardized PCA PC1-PC2 space.
- `spec_spca_convex`: convex-hull area of retained vector-normalized spectra in standardized PCA PC1-PC2 space. Archived diagnostic only; do not carry forward for future analysis.
- `spec_spca_hull3d_v`: convex-hull volume of retained vector-normalized spectra in standardized PCA PC1-PC3 space. Archived diagnostic only; do not carry forward for future analysis.
- `spec_spca_hull3d_a`: convex-hull surface area of retained vector-normalized spectra in standardized PCA PC1-PC3 space. Archived diagnostic only; do not carry forward for future analysis.

### Biodiversity Metrics

- `sp_rich`: species richness, calculated as the number of species with positive crown-overlap abundance in a quadrat.
- `sp_shannon`: Shannon diversity, calculated from proportional canopy-weighted species abundances.
- `sp_simpson`: Simpson diversity, calculated as one minus the sum of squared proportional species abundances.
- `sp_even`: species evenness, calculated as Shannon diversity divided by log species richness where richness is greater than one.
- `phy_faith`: Faith's phylogenetic diversity, calculated as the total phylogenetic branch length represented by species present in a quadrat.
- `phy_rao`: phylogenetic Rao's Q, calculated from canopy-weighted species abundances and phylogenetic distances among species.
- `phy_afaith`: abundance-weighted Faith's PD, calculated by weighting phylogenetic branch-length contributions by species crown-overlap abundance.

### Context Variables

- `env_elev`: mean quadrat elevation from the leaf-off digital terrain model.
- `env_tri5`, `env_tri11`, `env_tri21`: mean Riley topographic roughness index values from 5 x 5, 11 x 11, and 21 x 21 moving-window DTM products.
- Pixel brightness summaries: retained sunlit pixel brightness at 563 nm and in blue, green, red, and near-infrared wavelength regions.

## PCA Emphasis

The revised manuscript should focus primarily on PCA-based spectral heterogeneity, especially vector-normalized PCA values.

Vector-normalized PCA metrics should lead because they reduce broad brightness dominance and better isolate spectral shape. Raw PCA metrics should still be discussed because they show how much illumination and overall reflectance magnitude can influence spectral heterogeneity. The contrast between raw and vector-normalized PCA is part of the story: preprocessing choices affect what the spectral heterogeneity metric is actually measuring.

Recommended hierarchy:

1. Lead with vector-normalized PCA metrics, especially `spec_spca_alpha` and `spec_spca_mean`.
2. Use `spec_spca_rao` to discuss how squared-distance and pairwise-diversity formulations behave differently.
3. Use raw PCA metrics as a brightness-sensitive comparison.
4. Use spectral angle entropy as an important non-PCA reference metric, not as the center of the new title.

## Scale As The Center Of The Paper

Scale should be the organizing principle for the Results and Discussion.

At 10 m, spectral heterogeneity may be dominated by individual crowns, local illumination, crown geometry, geolocation mismatch, and within-crown spectral variation. At 20 m, quadrats begin to integrate community structure while retaining enough sample size for robust comparison. At 50 m, correlations are strongest, suggesting that broader aggregation may better align spectral summaries with canopy composition and phylogenetic structure, though the reduced number of quadrats requires caution.

The paper should explicitly evaluate the optimal quadrat scale for phylogenetic versus spectral correlation. Current evidence points to 50 m for correlation strength, but 20 m may offer a better balance between stronger relationships and higher sample size.

## Metric Behavior Story

The paper should explain metric strengths and weaknesses rather than only reporting which metric wins.

Mean PCA distance is stable and interpretable as average spectral departure from the quadrat centroid. It may capture dominant canopy spectral dispersion while being less controlled by rare spectral outliers.

PCA Rao's Q represents pairwise squared spectral dissimilarity. Its strength is conceptual alignment with diversity theory, but its weakness is greater sensitivity to extreme pixels and residual illumination or canopy-structure effects.

PCA alpha hull captures the occupied area of the quadrat's spectral cloud in PC1-PC2 space. Its current strength is high correlation with phylogenetic metrics, especially at 20 m and 50 m. Its weakness is potential sensitivity to hull geometry, outlying points, and the chosen PCA axes.

Convex hull and 3D hull volume/area metrics were evaluated as hull-family diagnostics, but the current decision is not to move them forward. They should remain archived screening outputs because they primarily summarize broad envelope size, can include unoccupied spectral space, and add complexity without becoming primary manuscript metrics.

Raw PCA metrics are useful because they show the influence of brightness and overall reflectance magnitude. Their weakness is that PC1 can be strongly illumination dominated.

Vector-normalized PCA metrics are stronger candidates for biological interpretation because they reduce brightness magnitude effects and emphasize spectral shape. Their weakness is that vector normalization does not remove all illumination, topographic, canopy-structural, or sensor-related variation.

Species richness and Faith's PD are presence-based and therefore tend to converge. Shannon, Simpson, and evenness emphasize abundance distribution and dominance. Phylogenetic Rao's Q and abundance-weighted Faith's PD can diverge from species metrics because they combine abundance structure with evolutionary distances.

## Why Phylogenetic Metrics Can Outperform Shannon Diversity

Shannon diversity responds to the number of species and the evenness of their canopy-weighted abundances, but it treats species identities as interchangeable. A quadrat with two closely related canopy species and a quadrat with two distantly related canopy species can have similar Shannon values.

Phylogenetic diversity metrics retain information about evolutionary relatedness. If evolutionary history is linked to leaf traits, canopy architecture, phenology, pigments, water content, or other reflectance-relevant traits, phylogenetic diversity can track spectral differences that Shannon diversity cannot see.

This explains the central pattern: phylogenetic diversity can show weak-to-moderate relationships with spectral heterogeneity while Shannon diversity remains weak, near zero, or inconsistent. The spectral signal is probably not just "more species equals more spectral variation." It is more likely a mixture of species identity, abundance dominance, phylogenetic dissimilarity, canopy structure, illumination, and scale.

## Transformation Question

The transformation question should become a targeted analysis rather than a speculative statement.

Test whether log, square-root, Box-Cox, or simple power transforms improve the relationship between each PCA spectral heterogeneity metric and each biodiversity metric at each scale. This should be assessed with Pearson r, Spearman r, R2, residual diagnostics, and interpretability. The goal is not simply to maximize correlation, but to determine whether metric distributions and scale-dependent relationships become more linear and ecologically interpretable.

Candidate transformations:

- Log or `log1p` for right-skewed positive metrics such as hull areas, Rao's Q values, and some phylogenetic metrics.
- Square-root for positive skew when zeros or near-zeros make log transformation awkward.
- Box-Cox or Yeo-Johnson for systematic distribution correction.
- Standardization within scale for comparing effect sizes across metrics without changing correlation structure.

## Final Manuscript Argument

The final paper should argue that support for the Spectral Variation Hypothesis depends on:

1. Spatial grain.
2. Spectral metric construction.
3. Biodiversity metric construction.
4. Preprocessing choices such as vector normalization.
5. Environmental and illumination context.

The strongest current story is that PCA-based spectral heterogeneity, especially vector-normalized metrics, captures a modest but meaningful component of canopy biodiversity at broader spatial scales. Phylogenetic metrics show stronger relationships than Shannon diversity because they retain information about lineage dissimilarity and abundance-weighted evolutionary structure, while Shannon diversity mainly summarizes abundance evenness without species identity or phylogenetic distance.

## Output Direction

The next manuscript outputs should be:

1. A revised outline using the working title and scale-centered questions.
2. A Results section that begins with metric behavior and scale, then moves into spectral-biodiversity correlations.
3. A Methods section that clearly documents crown weighting, abundance-weighted Faith's PD, PCA basis construction, vector normalization, raw versus standardized PCA, and metric formulas.
4. A transformation analysis report.
5. A spatial autocorrelation report using Moran's I, variograms, or spatial residual diagnostics.
6. A figure/table inventory that separates main-text outputs from supplemental metric screens.

## Source Context

This narrative incorporates the original research objectives, the current final direction file, the current `SVH_Paper.docx` draft, and Evan's 2026-08-15 next-step suggestions. Evan's notes particularly support the scale-centered framing, all-metric biodiversity comparison, biodiversity-metric correlation layer, spatial autocorrelation check, and fuller documentation of crown weighting and abundance-weighted Faith's PD.
