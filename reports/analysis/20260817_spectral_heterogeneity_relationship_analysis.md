# Spectral Heterogeneity Metric Relationship Analysis

Date: 2026-08-17

## Purpose

This analysis compares all seven focal spectral heterogeneity measures across 10 m, 20 m, and 50 m quadrats: spectral angle entropy, raw PCA alpha-hull area, raw PCA spectral Rao's Q, raw PCA mean Euclidean distance, standardized PCA alpha-hull area, standardized PCA spectral Rao's Q, and standardized PCA mean Euclidean distance.

## Methods

- Pairwise scatterplots and correlations were calculated for every unique pair among the seven focal spectral heterogeneity measures within each scale.
- Additional figure sets color the same pairwise relationships by mean elevation, species presence/absence composition type, crown-equivalent individuals, number of species present, mean retained-pixel brightness at 563 nm, and retained-pixel brightness in blue, green, red, and near-infrared spectral regions.
- Elevation panels report the incremental R2 and p-value for adding `env_elev` after the x-axis spectral metric.
- Composition-type panels report panel-specific mean silhouette width calculated from the two plotted spectral metrics.

## Strongest Spectral Metric Correlations

| Scale | Relationship | n | r | R2 | p-value |
| --- | --- | --- | --- | --- | --- |
| 10m | Std PCA Rao's Q vs Std PCA mean distance | 1744 | 0.907 | 0.823 | 0.00e+00 |
| 10m | Raw PCA Rao's Q vs Raw PCA mean distance | 1744 | 0.906 | 0.821 | 0.00e+00 |
| 10m | Raw PCA Rao's Q vs Std PCA Rao's Q | 1744 | 0.855 | 0.732 | 0.00e+00 |
| 20m | Raw PCA Rao's Q vs Raw PCA mean distance |  436 | 0.926 | 0.858 | 3.90e-186 |
| 20m | Std PCA alpha hull vs Std PCA mean distance |  436 | 0.904 | 0.817 | 3.48e-162 |
| 20m | Std PCA Rao's Q vs Std PCA mean distance |  436 | 0.898 | 0.806 | 1.72e-156 |
| 50m | Raw PCA Rao's Q vs Raw PCA mean distance |   74 | 0.940 | 0.883 | 2.91e-35 |
| 50m | Std PCA alpha hull vs Std PCA mean distance |   74 | 0.929 | 0.864 | 7.08e-33 |
| 50m | Std PCA Rao's Q vs Std PCA mean distance |   74 | 0.900 | 0.810 | 1.13e-27 |

## Largest Elevation Contributions

| Scale | Relationship | n | Elevation delta R2 | Elevation p-value |
| --- | --- | --- | --- | --- |
| 50m | Raw PCA mean distance vs Std PCA alpha hull |  74 | 0.109 | 2.35e-05 |
| 50m | Raw PCA Rao's Q vs Std PCA alpha hull |  74 | 0.106 | 4.57e-05 |
| 50m | Raw PCA mean distance vs Std PCA mean distance |  74 | 0.103 | 1.98e-05 |
| 50m | Raw PCA Rao's Q vs Std PCA mean distance |  74 | 0.088 | 5.70e-06 |
| 20m | Raw PCA mean distance vs Std PCA alpha hull | 436 | 0.081 | 6.10e-14 |
| 20m | Raw PCA Rao's Q vs Std PCA alpha hull | 436 | 0.075 | 1.22e-12 |
| 20m | Raw PCA mean distance vs Std PCA mean distance | 436 | 0.065 | 1.18e-13 |
| 50m | Raw PCA alpha hull vs Std PCA mean distance |  74 | 0.065 | 2.33e-04 |
| 20m | Raw PCA alpha hull vs Std PCA alpha hull | 436 | 0.063 | 2.96e-16 |
| 50m | Raw PCA alpha hull vs Std PCA alpha hull |  74 | 0.057 | 3.90e-05 |

## Strongest Composition-Type Separation In Spectral Metric Space

| Scale | Relationship | Overall composition silhouette |
| --- | --- | --- |
| 50m | Spectral angle entropy vs Std PCA mean distance | -0.006 |
| 50m | Spectral angle entropy vs Std PCA alpha hull | -0.011 |
| 50m | Spectral angle entropy vs Raw PCA mean distance | -0.019 |
| 50m | Std PCA alpha hull vs Std PCA mean distance | -0.022 |
| 50m | Spectral angle entropy vs Raw PCA alpha hull | -0.025 |
| 50m | Raw PCA mean distance vs Std PCA mean distance | -0.027 |
| 50m | Spectral angle entropy vs Raw PCA Rao's Q | -0.027 |
| 50m | Raw PCA alpha hull vs Std PCA mean distance | -0.028 |
| 50m | Raw PCA mean distance vs Std PCA alpha hull | -0.030 |
| 50m | Spectral angle entropy vs Std PCA Rao's Q | -0.033 |

## Figures

- `reports/figures/spectral_heterogeneity_relationships/01_spectral_metric_pairwise_scatter_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/01_spectral_metric_pairwise_scatter_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/01_spectral_metric_pairwise_scatter_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/02_spectral_metric_pairwise_scatter_elevation_gradient_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/02_spectral_metric_pairwise_scatter_elevation_gradient_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/02_spectral_metric_pairwise_scatter_elevation_gradient_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/03_spectral_metric_pairwise_scatter_composition_cluster_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/03_spectral_metric_pairwise_scatter_composition_cluster_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/03_spectral_metric_pairwise_scatter_composition_cluster_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/04_spectral_metric_pairwise_scatter_individual_ramp_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/04_spectral_metric_pairwise_scatter_individual_ramp_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/04_spectral_metric_pairwise_scatter_individual_ramp_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/05_spectral_metric_pairwise_scatter_species_count_ramp_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/05_spectral_metric_pairwise_scatter_species_count_ramp_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/05_spectral_metric_pairwise_scatter_species_count_ramp_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/06_spectral_metric_pairwise_scatter_pixel_brightness_ramp_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/06_spectral_metric_pairwise_scatter_pixel_brightness_ramp_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/06_spectral_metric_pairwise_scatter_pixel_brightness_ramp_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/07_spectral_metric_pairwise_scatter_blue_brightness_ramp_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/07_spectral_metric_pairwise_scatter_blue_brightness_ramp_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/07_spectral_metric_pairwise_scatter_blue_brightness_ramp_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/08_spectral_metric_pairwise_scatter_green_brightness_ramp_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/08_spectral_metric_pairwise_scatter_green_brightness_ramp_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/08_spectral_metric_pairwise_scatter_green_brightness_ramp_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/09_spectral_metric_pairwise_scatter_red_brightness_ramp_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/09_spectral_metric_pairwise_scatter_red_brightness_ramp_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/09_spectral_metric_pairwise_scatter_red_brightness_ramp_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/10_spectral_metric_pairwise_scatter_nir_brightness_ramp_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/10_spectral_metric_pairwise_scatter_nir_brightness_ramp_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/10_spectral_metric_pairwise_scatter_nir_brightness_ramp_50m.png`
- `reports/figures/spectral_heterogeneity_relationships/11_spectral_metrics_vs_regional_illumination_10m.png`
- `reports/figures/spectral_heterogeneity_relationships/11_spectral_metrics_vs_regional_illumination_20m.png`
- `reports/figures/spectral_heterogeneity_relationships/11_spectral_metrics_vs_regional_illumination_50m.png`

## Output Tables

- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_pairwise_correlations.csv`
- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_elevation_adjusted_models.csv`
- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_summary.csv`
- `reports/tables/spectral_heterogeneity_relationships/spectral_composition_scatterplot_silhouette.csv`
- `reports/tables/spectral_heterogeneity_relationships/quadrat_pixel_brightness_summary.csv`
- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_regional_illumination_correlations.csv`
