# Working Manuscript Outline

Created: 2026-07-30

Working title: Evaluating relationships between drone hyperspectral spectral heterogeneity and phylogenetic diversity in the Paint Rock Forest Dynamics Plot

## Manuscript Purpose

Develop a publication manuscript around the current reproducible workflow linking quadrat-scale hyperspectral spectral variation to field-derived species and phylogenetic diversity across 10 m, 20 m, and 50 m spatial grains.

## Central Question

Can drone-based hyperspectral spectral heterogeneity serve as a proxy for phylogenetic diversity in a temperate forest dynamics plot, and does that relationship depend on spatial grain?

## Working Thesis

Current results support a weak-to-moderate positive relationship between hyperspectral spectral heterogeneity and phylogenetic diversity, with stronger relationships at coarser spatial grains and clearer support for phylogenetic metrics than for standard species-diversity metrics. Standardized PCA mean Euclidean distance and spectral angle entropy both show positive relationships with phylogenetic diversity, but standardized PCA distance generally produces stronger scale-dependent relationships.

## Candidate Abstract Structure

1. Context: Biodiversity monitoring needs spatially continuous indicators, and hyperspectral spectral heterogeneity may capture canopy functional and phylogenetic variation.
2. Objective: Test relationships between drone hyperspectral spectral heterogeneity and field-derived biodiversity metrics across 10 m, 20 m, and 50 m quadrats.
3. Methods: Use sunlit, shadow-masked, smoothed 5 nm spectra; calculate spectral angle entropy and standardized PCA mean Euclidean distance; compare with species and phylogenetic diversity metrics derived from filtered tree census data.
4. Results: Relationships are consistently positive for phylogenetic metrics, strongest for phylogenetic Rao's Q and abundance-weighted Faith's PD, and weaker or absent for Shannon and other species-diversity metrics.
5. Implication: Hyperspectral spectral heterogeneity appears more aligned with phylogenetic structure than with species-diversity summaries, but spatial and environmental sensitivity checks are needed before final inference.

## Introduction

### Opening Problem

- Forest biodiversity is spatially complex and expensive to monitor solely from field surveys.
- Remote sensing can provide spatially continuous measurements, but links between spectral heterogeneity and biodiversity depend on scale, canopy structure, and the biological diversity metric being tested.

### Conceptual Background

- Spectral variation hypothesis: variation in reflectance may track variation in plant traits, canopy structure, species composition, or phylogenetic diversity.
- Phylogenetic diversity may be more mechanistically connected to spectral heterogeneity than species richness alone if related species share spectral or canopy trait structure.
- Spatial grain matters because 10 m, 20 m, and 50 m quadrats integrate different amounts of canopy, species mixing, and environmental heterogeneity.

### Study Gap

- Need direct multiscale tests of hyperspectral spectral heterogeneity against species and phylogenetic diversity in a mapped forest plot.
- Need transparent handling of shadow masking, spectral metric uncertainty, edge quadrats, and environmental context.

### Objectives And Hypotheses

- Objective 1: Quantify spectral heterogeneity from drone hyperspectral imagery at 10 m, 20 m, and 50 m grains.
- Objective 2: Quantify species and phylogenetic diversity from forest census data at matching grains.
- Objective 3: Test direct pairwise relationships between spectral heterogeneity and biodiversity metrics.
- Objective 4: Evaluate whether relationships are sensitive to edge quadrats, bootstrap variation, and environmental covariates.
- Hypothesis 1: Spectral heterogeneity increases with phylogenetic diversity.
- Hypothesis 2: Relationships strengthen with coarser spatial grain.
- Hypothesis 3: Phylogenetic diversity metrics have stronger positive relationships with spectral heterogeneity than species-diversity metrics.

## Methods

### Study Site

- Paint Rock Forest Dynamics Plot.
- Include location, forest type, plot area, terrain context, and census background once final site text is assembled.

### Field Biodiversity Data

- Source: `PR_tree_DL.csv` and supporting taxonomy/phylogeny files.
- Current filtered analysis population: upper-canopy individuals based on DBH/crown-position logic documented in the project state and plant-diversity workflow.
- Diversity metrics:
  - Species richness
  - Shannon diversity
  - Simpson diversity
  - Species evenness
  - Faith's phylogenetic diversity
  - Phylogenetic Rao's Q
  - Abundance-weighted Faith's phylogenetic diversity

### Hyperspectral Data

- Sensor: Headwall Micro A-Series VNIR.
- Spectral range: approximately 396-1002 nm.
- Processing basis: confirmed quadrat spectra in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m`.
- Current spectral inputs: Savitzky-Golay smoothed and 5 nm resampled spectra in `Quad_Spectra/*_smooth_5nm`.
- Shadow masking: 563 nm threshold of 0.0305476, retaining sunlit pixels greater than the threshold.

### Spectral Heterogeneity Metrics

- Primary metric 1: `spec_spca_mean`, standardized PCA mean Euclidean distance in PC1-PC3 space after vector-normalizing spectra.
- Primary metric 2: `spec_sa`, spectral angle entropy from sunlit, shadow-masked, smoothed 5 nm spectra.
- Supporting metrics: raw PCA distance, standardized PCA median/sd, spectral Rao's Q, and alpha-hull area. Convex-hull and 3D hull volume/area summaries are archived diagnostics only and should not be carried forward for future analysis.

### Spatial Grains

- 10 m quadrats
- 20 m quadrats
- 50 m quadrats
- Current analysis-ready tables:
  - `quadrat_analysis_10m.csv`
  - `quadrat_analysis_20m.csv`
  - `quadrat_analysis_50m.csv`

### Statistical Analysis

- First-layer analysis: simple pairwise models within scale:
  - `spectral_variation ~ diversity_measure`
- Report Pearson `r`, `R2`, F statistic, F-test p-value, slope, intercept, and Spearman rank correlation.
- Primary inference currently excludes documented 10 m and 20 m edge quadrats; all 50 m quadrats are retained because no 50 m edge rule is documented.
- Sensitivity checks:
  - Edge-quadrat inclusion/exclusion
  - Equal-sample-size non-edge resampling
  - SA entropy bootstrap CV subsets
  - Environmental adjustment with elevation and TRI variables

## Results

### Data Coverage

- Primary non-edge quadrats:
  - 10 m: 1,656 primary quadrats
  - 20 m: 414 primary quadrats
  - 50 m: 80 quadrats
- Complete standardized PCA distance values:
  - 10 m: 1,440
  - 20 m: 360
  - 50 m: 74
- Complete SA entropy values:
  - 10 m: 1,587
  - 20 m: 405
  - 50 m: 80

### Main Pairwise Findings

- Strongest 10 m standardized PCA pairing: phylogenetic Rao's Q, r = 0.121, R2 = 0.015.
- Strongest 10 m SA entropy pairing: phylogenetic Rao's Q, r = 0.106, R2 = 0.011.
- Strongest 20 m standardized PCA pairing: abundance-weighted Faith's PD, r = 0.291, R2 = 0.085.
- Strongest 20 m SA entropy pairing: phylogenetic Rao's Q, r = 0.180, R2 = 0.033.
- Strongest 50 m standardized PCA pairing: phylogenetic Rao's Q, r = 0.444, R2 = 0.197.
- Strongest 50 m SA entropy pairing: phylogenetic Rao's Q, r = 0.340, R2 = 0.115.

### Scale Dependence

- Relationships are weakest at 10 m and strongest at 50 m.
- The 50 m standardized PCA relationship with phylogenetic Rao's Q is the current strongest primary result.
- Coarser grains may better integrate canopy mixtures and phylogenetic composition, but 50 m has only 80 quadrats and should be interpreted with sample-size caution.

### Phylogenetic Versus Species Diversity

- Phylogenetic metrics show the most consistent positive relationships with both spectral heterogeneity measures.
- Species Shannon, Simpson, evenness, and richness are weakly related, near zero, or negative in most primary comparisons.

### Sensitivity Findings

- Removing edge quadrats generally weakens all-quadrat correlations slightly at 10 m and 20 m, but the primary phylogenetic patterns remain.
- Edge-only correlations can appear stronger, especially at 20 m, but range restriction and environmental context need careful interpretation.
- Elevation is negatively correlated with spectral variation, especially at 20 m.
- Environmental-adjusted screening models suggest phylogenetic diversity retains some incremental explanatory value, particularly for phylogenetic Rao's Q and abundance-weighted Faith's PD.
- For SA entropy, low bootstrap-CV subsets retain or strengthen phylogenetic correlations, while high-CV subsets show weak or absent relationships.

## Discussion

### Interpretation Of Main Result

- The current evidence suggests hyperspectral spectral heterogeneity is more closely aligned with phylogenetic diversity than with conventional species-diversity metrics.
- Standardized PCA mean distance may capture biologically relevant spectral structure after reducing brightness dominance.
- Spectral angle entropy provides an independent angular spectral heterogeneity metric that supports the same broad phylogenetic pattern.

### Why Phylogenetic Metrics May Perform Better

- Phylogenetic structure may integrate trait differences that are spectrally expressed through canopy chemistry, leaf structure, crown architecture, or phenology.
- Species counts and Shannon diversity may miss functional or evolutionary dissimilarity among co-occurring canopy species.

### Why Scale Matters

- 10 m quadrats may be too fine to capture stable canopy-community relationships.
- 20 m provides intermediate support with enough sample size for sensitivity checks.
- 50 m shows the strongest associations but has fewer quadrats, so spatial non-independence and leverage need final checks.

### Limitations

- Spectral metrics are affected by shadow masking, atmospheric/cloud exclusions, edge quadrats, and retained-pixel counts.
- Current direct pairwise models do not yet fully handle spatial autocorrelation.
- Environmental covariates are screening adjustments, not final spatial inference.
- Bootstrap uncertainty is fully developed for SA entropy but not for standardized PCA metrics.

### Implications

- Drone hyperspectral spectral heterogeneity may be useful for mapping relative phylogenetic diversity patterns in forest canopies.
- The most defensible manuscript claim should emphasize scale-dependent, modest, phylogenetically stronger relationships rather than a universal biodiversity proxy.

## Figures To Assemble

- Study site and quadrat-scale map.
- Workflow diagram from hyperspectral imagery and tree census to combined quadrat tables.
- Spectral heterogeneity maps for `spec_spca_mean` and `spec_sa`.
- Phylogenetic diversity maps for `phy_rao`, `phy_faith`, and `phy_afaith`.
- Pairwise correlation heatmap across scales.
- Scatterplots for primary pairings:
  - `spec_spca_mean ~ phy_rao`
  - `spec_spca_mean ~ phy_afaith`
  - `spec_sa ~ phy_rao`
- Edge/sensitivity figure as supplement or robustness panel.
- PCA loading spectral-region figure as methods/results support.

## Tables To Assemble

- Table 1: Study data and quadrat coverage by scale.
- Table 2: Diversity and spectral metric definitions.
- Table 3: Strongest primary pairwise relationships by scale and spectral metric.
- Table 4: Full primary pairwise correlation results, or supplement if too large.
- Supplementary table: Edge, environmental, and bootstrap sensitivity summaries.

## Current Source Reports

- `reports/analysis/20260710_sv_diversity_pairwise_correlations.md`
- `reports/analysis/20260725_edge_bootstrap_sensitivity_sv_diversity.md`
- `reports/analysis/20260707_pca_loading_spectral_regions.md`
- `reports/analysis/20260707_standardized_pca_metric_sample_size_effects.md`
- `reports/analysis/20260618_bootstrap_variation_analysis.md`
- `reports/combined_quadrat_variable_guide.md`
- `reports/data_dictionary.md`
- `reports/project_state.md`

## Open Decisions Before Drafting

- Target journal and required manuscript format.
- Whether to lead with standardized PCA mean distance, SA entropy, or a paired two-metric framing.
- Whether 50 m results are framed as primary despite lower sample size, or as a scale-dependent endpoint.
- Whether edge quadrats are excluded from all main figures or shown but visually flagged.
- Which spatial model or spatial autocorrelation diagnostic becomes the final inferential layer.
- Whether environmental models stay in supplemental sensitivity or enter the main results.
