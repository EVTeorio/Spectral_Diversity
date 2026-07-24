# SV-Diversity Pairwise Correlation Analysis

Date: 2026-07-10

## Research Question

What is the relationship between spectral variation and phylogenetic/species diversity across 10 m, 20 m, and 50 m quadrats?

This analysis is intentionally direct: each primary spectral variation measure is independently paired with each phylogenetic and species diversity measure. No environmental covariates, interaction terms, or multivariable model-selection steps are included in this first relationship layer.

## Primary Spectral Variation Measures

- `spec_spca_mean`: standardized PCA mean Euclidean distance in PC1-PC3 space after vector-normalizing spectra.
- `spec_sa`: spectral angle entropy mean from sunlit, shadow-masked smoothed 5 nm spectra.

## Diversity Measures

- Phylogenetic diversity: `phy_faith`, `phy_rao`, `phy_afaith`.
- Species diversity: `sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even`.

## Model Form

For each scale, spectral variation measure, and diversity measure, the workflow fits:

`SV_measure ~ diversity_measure`

Reported values include Pearson `r`, `R2`, simple-regression `F`, F-test p-value, slope, intercept, and Spearman rank correlation. In these one-predictor models, `R2` is the squared Pearson correlation.

`sv_diversity_analysis_dataset.csv` also includes sampling/provenance flags from the upstream spectral heterogeneity summaries. `sa_all_pixels_sampled` is TRUE when the SA entropy workflow sampled all retained pixels for that quadrat and FALSE when the 5,000-pixel cap was used. `spca_euclidean_all_pixels_sampled` is TRUE when the standardized PCA Euclidean-distance metric used all retained pixels for that quadrat.

Documented 10 m and 20 m edge quadrats are excluded from primary analysis; all 50 m quadrats are retained because no separate 50 m edge rule is documented.

## Coverage

| Scale | Total quadrats | Primary quadrats | Edge flagged | Complete SA entropy | Complete standardized PCA distance | SA all pixels sampled | Std PCA all pixels sampled |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | 2000 | 1656 | 344 | 1587 | 1440 | 108 | 1440 |
| 20m |  500 |  414 |  86 |  405 |  360 |  12 |  360 |
| 50m |   80 |   80 |   0 |   80 |   74 |   1 |   74 |

## Strongest Pairing Per SV Measure And Scale

| Scale | SV measure | Strongest diversity pairing | Group | n | r | R2 | F | F p-value |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | Standardized PCA mean Euclidean distance | Phylogenetic Rao's Q | Phylogenetic diversity | 1440 | 0.121 | 0.015 | 21.54 | 3.77e-06 |
| 10m | SA entropy | Phylogenetic Rao's Q | Phylogenetic diversity | 1587 | 0.106 | 0.011 | 18.01 | 2.33e-05 |
| 20m | Standardized PCA mean Euclidean distance | Abundance-weighted Faith's PD | Phylogenetic diversity |  360 | 0.291 | 0.085 | 33.17 | 1.82e-08 |
| 20m | SA entropy | Phylogenetic Rao's Q | Phylogenetic diversity |  405 | 0.180 | 0.033 | 13.55 | 2.64e-04 |
| 50m | Standardized PCA mean Euclidean distance | Phylogenetic Rao's Q | Phylogenetic diversity |   74 | 0.444 | 0.197 | 17.65 | 7.49e-05 |
| 50m | SA entropy | Phylogenetic Rao's Q | Phylogenetic diversity |   80 | 0.340 | 0.115 | 10.17 | 0.002 |

## Full Pairwise Results

| Scale | SV measure | Diversity measure | Group | n | r | R2 | F | F p-value |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | Standardized PCA mean Euclidean distance | Faith's PD | Phylogenetic diversity | 1440 | 0.098 | 0.010 | 14.01 | 1.89e-04 |
| 10m | Standardized PCA mean Euclidean distance | Phylogenetic Rao's Q | Phylogenetic diversity | 1440 | 0.121 | 0.015 | 21.54 | 3.77e-06 |
| 10m | Standardized PCA mean Euclidean distance | Abundance-weighted Faith's PD | Phylogenetic diversity | 1440 | 0.106 | 0.011 | 16.36 | 5.51e-05 |
| 10m | Standardized PCA mean Euclidean distance | Species richness | Species diversity | 1440 | -0.073 | 0.005 | 7.77 | 0.005 |
| 10m | Standardized PCA mean Euclidean distance | Shannon diversity | Species diversity | 1440 | -0.049 | 0.002 | 3.49 | 0.062 |
| 10m | Standardized PCA mean Euclidean distance | Simpson diversity | Species diversity | 1440 | -0.039 | 0.002 | 2.19 | 0.139 |
| 10m | Standardized PCA mean Euclidean distance | Species evenness | Species diversity | 1440 | -0.016 | 0.000 | 0.35 | 0.556 |
| 10m | SA entropy | Faith's PD | Phylogenetic diversity | 1587 | 0.088 | 0.008 | 12.26 | 4.76e-04 |
| 10m | SA entropy | Phylogenetic Rao's Q | Phylogenetic diversity | 1587 | 0.106 | 0.011 | 18.01 | 2.33e-05 |
| 10m | SA entropy | Abundance-weighted Faith's PD | Phylogenetic diversity | 1587 | 0.102 | 0.010 | 16.75 | 4.47e-05 |
| 10m | SA entropy | Species richness | Species diversity | 1587 | -0.048 | 0.002 | 3.72 | 0.054 |
| 10m | SA entropy | Shannon diversity | Species diversity | 1587 | -0.051 | 0.003 | 4.10 | 0.043 |
| 10m | SA entropy | Simpson diversity | Species diversity | 1587 | -0.046 | 0.002 | 3.33 | 0.068 |
| 10m | SA entropy | Species evenness | Species diversity | 1587 | -0.035 | 0.001 | 1.93 | 0.165 |
| 20m | Standardized PCA mean Euclidean distance | Faith's PD | Phylogenetic diversity |  360 | 0.222 | 0.049 | 18.51 | 2.19e-05 |
| 20m | Standardized PCA mean Euclidean distance | Phylogenetic Rao's Q | Phylogenetic diversity |  360 | 0.249 | 0.062 | 23.57 | 1.81e-06 |
| 20m | Standardized PCA mean Euclidean distance | Abundance-weighted Faith's PD | Phylogenetic diversity |  360 | 0.291 | 0.085 | 33.17 | 1.82e-08 |
| 20m | Standardized PCA mean Euclidean distance | Species richness | Species diversity |  360 | -0.067 | 0.005 | 1.64 | 0.201 |
| 20m | Standardized PCA mean Euclidean distance | Shannon diversity | Species diversity |  360 | -0.075 | 0.006 | 2.03 | 0.155 |
| 20m | Standardized PCA mean Euclidean distance | Simpson diversity | Species diversity |  360 | -0.063 | 0.004 | 1.43 | 0.232 |
| 20m | Standardized PCA mean Euclidean distance | Species evenness | Species diversity |  360 | -0.020 | 0.000 | 0.14 | 0.707 |
| 20m | SA entropy | Faith's PD | Phylogenetic diversity |  405 | 0.140 | 0.020 | 8.10 | 0.005 |
| 20m | SA entropy | Phylogenetic Rao's Q | Phylogenetic diversity |  405 | 0.180 | 0.033 | 13.55 | 2.64e-04 |
| 20m | SA entropy | Abundance-weighted Faith's PD | Phylogenetic diversity |  405 | 0.166 | 0.027 | 11.39 | 8.12e-04 |
| 20m | SA entropy | Species richness | Species diversity |  405 | -0.029 | 0.001 | 0.33 | 0.566 |
| 20m | SA entropy | Shannon diversity | Species diversity |  405 | -0.065 | 0.004 | 1.69 | 0.195 |
| 20m | SA entropy | Simpson diversity | Species diversity |  405 | -0.086 | 0.007 | 3.00 | 0.084 |
| 20m | SA entropy | Species evenness | Species diversity |  405 | -0.058 | 0.003 | 1.36 | 0.244 |
| 50m | Standardized PCA mean Euclidean distance | Faith's PD | Phylogenetic diversity |   74 | 0.282 | 0.080 | 6.23 | 0.015 |
| 50m | Standardized PCA mean Euclidean distance | Phylogenetic Rao's Q | Phylogenetic diversity |   74 | 0.444 | 0.197 | 17.65 | 7.49e-05 |
| 50m | Standardized PCA mean Euclidean distance | Abundance-weighted Faith's PD | Phylogenetic diversity |   74 | 0.419 | 0.175 | 15.31 | 2.05e-04 |
| 50m | Standardized PCA mean Euclidean distance | Species richness | Species diversity |   74 | 0.092 | 0.008 | 0.61 | 0.436 |
| 50m | Standardized PCA mean Euclidean distance | Shannon diversity | Species diversity |   74 | 0.009 | 0.000 | 0.01 | 0.942 |
| 50m | Standardized PCA mean Euclidean distance | Simpson diversity | Species diversity |   74 | -0.009 | 0.000 | 0.01 | 0.940 |
| 50m | Standardized PCA mean Euclidean distance | Species evenness | Species diversity |   74 | -0.107 | 0.012 | 0.84 | 0.363 |
| 50m | SA entropy | Faith's PD | Phylogenetic diversity |   80 | 0.171 | 0.029 | 2.34 | 0.130 |
| 50m | SA entropy | Phylogenetic Rao's Q | Phylogenetic diversity |   80 | 0.340 | 0.115 | 10.17 | 0.002 |
| 50m | SA entropy | Abundance-weighted Faith's PD | Phylogenetic diversity |   80 | 0.228 | 0.052 | 4.27 | 0.042 |
| 50m | SA entropy | Species richness | Species diversity |   80 | 0.072 | 0.005 | 0.40 | 0.527 |
| 50m | SA entropy | Shannon diversity | Species diversity |   80 | 0.047 | 0.002 | 0.17 | 0.678 |
| 50m | SA entropy | Simpson diversity | Species diversity |   80 | 0.001 | 0.000 | 0.00 | 0.996 |
| 50m | SA entropy | Species evenness | Species diversity |   80 | -0.033 | 0.001 | 0.09 | 0.769 |

## Figures

- `reports/figures/multiscale_spectral_biodiversity/01_sv_diversity_pairwise_correlation_heatmap.png`
- `reports/figures/multiscale_spectral_biodiversity/02_sa_entropy_diversity_scatterplots.png`
- `reports/figures/multiscale_spectral_biodiversity/03_standardized_pca_distance_diversity_scatterplots.png`

## Output Tables

- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_pairwise_correlations.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_top_pairings.csv`

## Superseded Analysis Direction

Earlier candidate-model tables that ranked environmental and interaction models are retained only as historical context if present in the output folder. The active direction for this stage is the direct pairwise relationship between the two primary SV measures and each diversity measure.
