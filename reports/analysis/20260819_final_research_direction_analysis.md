# Final Research Direction Analysis

Date: 2026-08-19

## Purpose

This analysis completes the current `Checklist.md` analysis-needed items that can be completed from existing 10 m, 20 m, and 50 m outputs. It compiles PCA-centered spectral-biodiversity relationships, compares spectral and biodiversity metrics, tests transformations, summarizes metric drivers, and runs spatial autocorrelation diagnostics for priority variables and residuals.

The discarded 100 m optimal-scale task is not run because no 100 m products are currently part of the validated workflow.

## Top PCA Spectral-Biodiversity Relationships

| Scale | Spectral | Biodiversity | n | r | R2 | p |
| --- | --- | --- | --- | --- | --- | --- |
| 50m | Std PCA alpha hull | Abundance-weighted Faith's PD |  74 | 0.497 | 0.247 | 6.79e-06 |
| 50m | Std PCA alpha hull | Phylogenetic Rao's Q |  74 | 0.487 | 0.238 | 1.06e-05 |
| 50m | Std PCA mean distance | Phylogenetic Rao's Q |  74 | 0.444 | 0.197 | 7.49e-05 |
| 50m | Std PCA mean distance | Abundance-weighted Faith's PD |  74 | 0.419 | 0.175 | 2.05e-04 |
| 50m | Raw PCA alpha hull | Abundance-weighted Faith's PD |  74 | 0.412 | 0.170 | 2.62e-04 |
| 20m | Std PCA alpha hull | Abundance-weighted Faith's PD | 436 | 0.384 | 0.147 | 9.97e-17 |
| 50m | Raw PCA mean distance | Abundance-weighted Faith's PD |  74 | 0.375 | 0.141 | 9.80e-04 |
| 50m | Std PCA Rao's Q | Phylogenetic Rao's Q |  74 | 0.352 | 0.124 | 0.002 |
| 50m | Raw PCA Rao's Q | Abundance-weighted Faith's PD |  74 | 0.346 | 0.119 | 0.003 |
| 20m | Std PCA alpha hull | Phylogenetic Rao's Q | 436 | 0.336 | 0.113 | 6.08e-13 |
| 50m | Raw PCA alpha hull | Phylogenetic Rao's Q |  74 | 0.332 | 0.111 | 0.004 |
| 50m | Raw PCA Rao's Q | Phylogenetic Rao's Q |  74 | 0.319 | 0.102 | 0.006 |
| 50m | Std PCA alpha hull | Faith's PD |  74 | 0.312 | 0.097 | 0.007 |
| 50m | Raw PCA 3D hull volume | Phylogenetic Rao's Q |  74 | 0.304 | 0.093 | 0.008 |
| 20m | Std PCA mean distance | Abundance-weighted Faith's PD | 436 | 0.301 | 0.091 | 1.36e-10 |
| 50m | Raw PCA mean distance | Phylogenetic Rao's Q |  74 | 0.292 | 0.085 | 0.012 |
| 20m | Std PCA alpha hull | Faith's PD | 436 | 0.285 | 0.081 | 1.40e-09 |
| 50m | Std PCA mean distance | Faith's PD |  74 | 0.282 | 0.080 | 0.015 |
| 20m | Std PCA mean distance | Phylogenetic Rao's Q | 436 | 0.274 | 0.075 | 6.16e-09 |
| 50m | Std PCA Rao's Q | Abundance-weighted Faith's PD |  74 | 0.272 | 0.074 | 0.019 |

## Largest Transformation Gains

| Scale | Spectral | Biodiversity | Best | Identity_r | Best_r | Delta_abs_r |
| --- | --- | --- | --- | --- | --- | --- |
| 10m | Std PCA Rao's Q | Abundance-weighted Faith's PD | power_2 x / yeojohnson y | 0.025 | 0.170 | 0.145 |
| 20m | Std PCA Rao's Q | Abundance-weighted Faith's PD | identity x / boxcox y | 0.133 | 0.269 | 0.136 |
| 10m | Std PCA Rao's Q | Faith's PD | power_2 x / yeojohnson y | 0.026 | 0.160 | 0.135 |
| 10m | Std PCA Rao's Q | Phylogenetic Rao's Q | power_2 x / yeojohnson y | 0.044 | 0.168 | 0.124 |
| 20m | Std PCA Rao's Q | Phylogenetic Rao's Q | identity x / yeojohnson y | 0.132 | 0.237 | 0.104 |
| 20m | Std PCA Rao's Q | Faith's PD | identity x / boxcox y | 0.123 | 0.226 | 0.103 |
| 10m | Std PCA 3D hull volume | Abundance-weighted Faith's PD | power_2 x / boxcox y | -0.003 | 0.087 | 0.085 |
| 10m | Std PCA 3D hull area | Abundance-weighted Faith's PD | power_2 x / boxcox y | 0.001 | 0.081 | 0.079 |
| 10m | Std PCA mean distance | Abundance-weighted Faith's PD | power_2 x / yeojohnson y | 0.122 | 0.197 | 0.075 |
| 10m | Std PCA mean distance | Faith's PD | power_2 x / yeojohnson y | 0.119 | 0.193 | 0.074 |
| 10m | Std PCA convex hull | Abundance-weighted Faith's PD | power_2 x / boxcox y | -0.004 | 0.076 | 0.072 |
| 50m | Std PCA Rao's Q | Abundance-weighted Faith's PD | sqrt x / boxcox y | 0.272 | 0.342 | 0.070 |
| 20m | Std PCA 3D hull volume | Abundance-weighted Faith's PD | power_2 x / boxcox y | 0.003 | 0.068 | 0.065 |
| 10m | Std PCA mean distance | Phylogenetic Rao's Q | power_2 x / yeojohnson y | 0.140 | 0.205 | 0.065 |
| 50m | Std PCA Rao's Q | Species evenness | log x / boxcox y | -0.027 | -0.091 | 0.064 |

## Strongest Metric Driver Relationships

| Scale | Target | Driver | Group | n | r | R2 | p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20m | Species richness | Present species count | biodiversity |  500 | 1.000 | 1.000 | 0.00e+00 |
| 50m | Species richness | Present species count | biodiversity |   80 | 1.000 | 1.000 | 0.00e+00 |
| 10m | Species richness | Present species count | biodiversity | 2000 | 1.000 | 1.000 | 0.00e+00 |
| 10m | Shannon diversity | Simpson diversity | biodiversity | 2000 | 0.969 | 0.939 | 0.00e+00 |
| 20m | Shannon diversity | Simpson diversity | biodiversity |  500 | 0.949 | 0.901 | 7.24e-252 |
| 50m | Shannon diversity | Simpson diversity | biodiversity |   80 | 0.939 | 0.883 | 5.06e-38 |
| 50m | Simpson diversity | Species evenness | biodiversity |   80 | 0.902 | 0.814 | 3.20e-30 |
| 50m | Species evenness | Simpson diversity | biodiversity |   80 | 0.902 | 0.814 | 3.20e-30 |
| 10m | Simpson diversity | Species evenness | biodiversity | 2000 | 0.860 | 0.740 | 0.00e+00 |
| 10m | Species evenness | Simpson diversity | biodiversity | 2000 | 0.860 | 0.740 | 0.00e+00 |
| 20m | Shannon diversity | Present species count | biodiversity |  500 | 0.856 | 0.732 | 1.32e-144 |
| 50m | Shannon diversity | Present species count | biodiversity |   80 | 0.847 | 0.718 | 3.82e-23 |
| 20m | Simpson diversity | Species evenness | biodiversity |  500 | 0.834 | 0.695 | 1.83e-130 |
| 20m | Species evenness | Simpson diversity | biodiversity |  500 | 0.834 | 0.695 | 1.83e-130 |
| 50m | Shannon diversity | Species evenness | biodiversity |   80 | 0.829 | 0.688 | 2.13e-21 |
| 10m | Shannon diversity | Present species count | biodiversity | 2000 | 0.824 | 0.679 | 0.00e+00 |
| 10m | Shannon diversity | Species evenness | biodiversity | 2000 | 0.763 | 0.582 | 0.00e+00 |
| 50m | Faith's PD | Present species count | biodiversity |   80 | 0.761 | 0.580 | 2.45e-16 |
| 10m | Raw PCA alpha hull | Green brightness | spectral | 1741 | 0.758 | 0.575 | 0.00e+00 |
| 50m | Spectral angle entropy | NIR brightness | spectral |   80 | -0.753 | 0.566 | 8.36e-16 |

## Strongest Spatial Autocorrelation Diagnostics

| Scale | Type | Variable_or_model | n | Moran_I | p |
| --- | --- | --- | --- | --- | --- |
| 10m | variable | Std PCA alpha hull | 1744 | 0.417 | 0.005 |
| 10m | residual | spec_spca_alpha ~ sp_shannon | 1744 | 0.416 | 0.005 |
| 10m | residual | spec_spca_alpha ~ phy_rao | 1744 | 0.394 | 0.005 |
| 10m | residual | spec_spca_alpha ~ phy_afaith | 1744 | 0.385 | 0.005 |
| 50m | variable | Phylogenetic Rao's Q |   80 | 0.373 | 0.005 |
| 20m | variable | Phylogenetic Rao's Q |  500 | 0.368 | 0.005 |
| 20m | variable | Abundance-weighted Faith's PD |  500 | 0.367 | 0.005 |
| 10m | variable | Abundance-weighted Faith's PD | 2000 | 0.349 | 0.005 |
| 20m | variable | Shannon diversity |  500 | 0.333 | 0.005 |
| 20m | residual | spec_spca_alpha ~ sp_shannon |  436 | 0.318 | 0.005 |
| 20m | variable | Std PCA alpha hull |  436 | 0.318 | 0.005 |
| 10m | variable | Phylogenetic Rao's Q | 2000 | 0.295 | 0.005 |
| 10m | variable | Shannon diversity | 2000 | 0.284 | 0.005 |
| 50m | variable | Abundance-weighted Faith's PD |   80 | 0.283 | 0.005 |
| 10m | variable | Std PCA mean distance | 1744 | 0.281 | 0.005 |
| 10m | residual | spec_spca_mean ~ sp_shannon | 1744 | 0.281 | 0.005 |
| 20m | residual | spec_spca_mean ~ sp_shannon |  436 | 0.263 | 0.005 |
| 20m | variable | Std PCA mean distance |  436 | 0.263 | 0.005 |
| 10m | residual | spec_spca_mean ~ phy_afaith | 1744 | 0.262 | 0.005 |
| 10m | residual | spec_spca_mean ~ phy_rao | 1744 | 0.262 | 0.005 |

## PCA Metric Computation Method Check

| Scale | Field | Method | n |
| --- | --- | --- | --- |
| 10m | raw_pca_metric | all_pixels | 1744 |
| 10m | raw_pca_metric | manual_excluded |  165 |
| 10m | standardized_pca_metric | all_pixels | 1744 |
| 10m | standardized_pca_metric | manual_excluded |  165 |
| 10m | raw_alpha_hull | all_pixels | 1734 |
| 10m | raw_alpha_hull | fallback_sampled_pixels |    7 |
| 10m | raw_alpha_hull | fallback_sampled_pixels_alpha_failed |    3 |
| 10m | raw_alpha_hull | manual_excluded |  165 |
| 10m | standardized_alpha_hull | all_pixels | 1744 |
| 10m | standardized_alpha_hull | manual_excluded |  165 |
| 10m | raw_3d_hull | all_pixels | 1744 |
| 10m | raw_3d_hull | manual_excluded |  165 |
| 10m | standardized_3d_hull | all_pixels | 1744 |
| 10m | standardized_3d_hull | manual_excluded |  165 |
| 20m | raw_pca_metric | all_pixels |  436 |
| 20m | raw_pca_metric | manual_excluded |   49 |
| 20m | standardized_pca_metric | all_pixels |  436 |
| 20m | standardized_pca_metric | manual_excluded |   49 |
| 20m | raw_alpha_hull | all_pixels |    2 |
| 20m | raw_alpha_hull | manual_excluded |   49 |
| 20m | raw_alpha_hull | sampled_pixels |  434 |
| 20m | standardized_alpha_hull | all_pixels |    2 |
| 20m | standardized_alpha_hull | manual_excluded |   49 |
| 20m | standardized_alpha_hull | sampled_pixels |  434 |
| 20m | raw_3d_hull | all_pixels |    2 |
| 20m | raw_3d_hull | manual_excluded |   49 |
| 20m | raw_3d_hull | sampled_pixels |  434 |
| 20m | standardized_3d_hull | all_pixels |    2 |
| 20m | standardized_3d_hull | manual_excluded |   49 |
| 20m | standardized_3d_hull | sampled_pixels |  434 |
| 50m | raw_pca_metric | all_pixels |   74 |
| 50m | raw_pca_metric | manual_excluded |    6 |
| 50m | standardized_pca_metric | all_pixels |   74 |
| 50m | standardized_pca_metric | manual_excluded |    6 |
| 50m | raw_alpha_hull | manual_excluded |    6 |
| 50m | raw_alpha_hull | sampled_pixels |   74 |
| 50m | standardized_alpha_hull | manual_excluded |    6 |
| 50m | standardized_alpha_hull | sampled_pixels |   74 |
| 50m | raw_3d_hull | manual_excluded |    6 |
| 50m | raw_3d_hull | sampled_pixels |   74 |
| 50m | standardized_3d_hull | manual_excluded |    6 |
| 50m | standardized_3d_hull | sampled_pixels |   74 |

## Figures

- `reports/figures/final_research_direction/01_transformation_delta_abs_r_priority_relationships.png`
- `reports/figures/final_research_direction/02_spatial_moran_i_priority_variables_and_residuals.png`
- `reports/figures/final_research_direction/03_spectral_metric_correlation_heatmap.png`
- `reports/figures/final_research_direction/04_biodiversity_metric_correlation_heatmap.png`

## Output Tables

- `reports/tables/final_research_direction/all_spectral_biodiversity_relationships.csv`
- `reports/tables/final_research_direction/biodiversity_metric_pairwise_relationships.csv`
- `reports/tables/final_research_direction/metric_driver_relationships.csv`
- `reports/tables/final_research_direction/pca_metric_computation_method_summary.csv`
- `reports/tables/final_research_direction/spatial_moran_diagnostics.csv`
- `reports/tables/final_research_direction/spatial_variogram_summary.csv`
- `reports/tables/final_research_direction/spectral_metric_pairwise_relationships.csv`
- `reports/tables/final_research_direction/transformation_best_relationship_summary.csv`
- `reports/tables/final_research_direction/transformation_spectral_biodiversity_relationships.csv`

## Interpretation Notes

- Vector-normalized PCA metrics are reported as the primary PCA evidence layer, with raw PCA metrics retained as brightness-sensitive comparisons.
- Transformation gains should be interpreted as diagnostics, not as automatic replacements for untransformed metrics. A transform is useful only if it improves correlation and residual behavior while remaining interpretable.
- Moran's I diagnostics use eight-nearest-neighbor spatial weights and 199 permutations. Significant residual autocorrelation means coefficient p-values from ordinary least squares should remain cautious.
- Driver relationships are screening summaries. They identify plausible influences on metric values but do not by themselves establish causation.
