# Spectral Heterogeneity Metric Relationship Analysis

Date: 2026-08-17

## Purpose

This analysis compares all 13 focal spectral heterogeneity measures across 10 m, 20 m, and 50 m quadrats: spectral angle entropy plus raw and standardized PCA alpha-hull area, spectral Rao's Q, mean Euclidean distance, convex hull, 3D hull volume, and 3D hull area metrics.

## Methods

- Pairwise scatterplots and correlations were calculated for every unique pair among the 13 focal spectral heterogeneity measures within each scale.
- Additional figure sets color the same pairwise relationships by mean elevation, species presence/absence composition type, crown-equivalent individuals, number of species present, mean retained-pixel brightness at 563 nm, and retained-pixel brightness in blue, green, red, and near-infrared spectral regions.
- Elevation panels report the incremental R2 and p-value for adding `env_elev` after the x-axis spectral metric.
- Composition-type panels report panel-specific mean silhouette width calculated from the two plotted spectral metrics.

## Strongest Spectral Metric Correlations

| Scale | Relationship | n | r | R2 | p-value |
| --- | --- | --- | --- | --- | --- |
| 10m | Std PCA convex hull vs Std PCA 3D hull area | 1744 | 0.992 | 0.985 | 0.00e+00 |
| 10m | Raw PCA convex hull vs Raw PCA 3D hull area | 1744 | 0.991 | 0.983 | 0.00e+00 |
| 10m | Std PCA 3D hull volume vs Std PCA 3D hull area | 1744 | 0.980 | 0.961 | 0.00e+00 |
| 20m | Std PCA 3D hull volume vs Std PCA 3D hull area |  436 | 0.981 | 0.962 | 2.65e-311 |
| 20m | Raw PCA 3D hull volume vs Raw PCA 3D hull area |  436 | 0.980 | 0.961 | 1.74e-307 |
| 20m | Raw PCA convex hull vs Raw PCA 3D hull area |  436 | 0.976 | 0.952 | 4.66e-288 |
| 50m | Std PCA 3D hull volume vs Std PCA 3D hull area |   74 | 0.975 | 0.952 | 4.61e-49 |
| 50m | Raw PCA 3D hull volume vs Raw PCA 3D hull area |   74 | 0.975 | 0.950 | 1.45e-48 |
| 50m | Raw PCA Rao's Q vs Raw PCA mean distance |   74 | 0.940 | 0.883 | 2.91e-35 |

## Largest Elevation Contributions

| Scale | Relationship | n | Elevation delta R2 | Elevation p-value |
| --- | --- | --- | --- | --- |
| 50m | Raw PCA convex hull vs Std PCA alpha hull |  74 | 0.227 | 3.03e-06 |
| 50m | Raw PCA convex hull vs Std PCA mean distance |  74 | 0.220 | 2.38e-07 |
| 50m | Raw PCA 3D hull area vs Std PCA alpha hull |  74 | 0.203 | 3.02e-06 |
| 50m | Raw PCA 3D hull area vs Std PCA mean distance |  74 | 0.191 | 2.26e-07 |
| 50m | Raw PCA 3D hull volume vs Std PCA alpha hull |  74 | 0.176 | 9.04e-06 |
| 50m | Raw PCA 3D hull volume vs Std PCA mean distance |  74 | 0.158 | 1.11e-06 |
| 20m | Raw PCA convex hull vs Std PCA alpha hull | 436 | 0.133 | 1.73e-19 |
| 20m | Raw PCA convex hull vs Std PCA mean distance | 436 | 0.128 | 1.91e-24 |
| 50m | Std PCA mean distance vs Std PCA 3D hull volume |  74 | 0.126 | 2.66e-04 |
| 20m | Raw PCA 3D hull area vs Std PCA alpha hull | 436 | 0.117 | 1.33e-17 |

## Strongest Composition-Type Separation In Spectral Metric Space

| Scale | Relationship | Overall composition silhouette |
| --- | --- | --- |
| 50m | Std PCA mean distance vs Std PCA 3D hull volume | 0.014 |
| 50m | Std PCA mean distance vs Std PCA 3D hull area | 0.009 |
| 50m | Std PCA alpha hull vs Std PCA 3D hull volume | -0.003 |
| 50m | Spectral angle entropy vs Std PCA mean distance | -0.006 |
| 50m | Std PCA alpha hull vs Std PCA 3D hull area | -0.009 |
| 50m | Spectral angle entropy vs Std PCA alpha hull | -0.011 |
| 50m | Std PCA Rao's Q vs Std PCA 3D hull volume | -0.015 |
| 50m | Spectral angle entropy vs Std PCA 3D hull area | -0.018 |
| 50m | Spectral angle entropy vs Raw PCA mean distance | -0.019 |
| 50m | Spectral angle entropy vs Std PCA 3D hull volume | -0.021 |

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
