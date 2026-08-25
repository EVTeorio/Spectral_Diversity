# Spectral Biodiversity All-Metric Scatter Analysis

Date: 2026-08-19

## Purpose

This analysis creates one scatterplot figure for each of the seven focal spectral heterogeneity measures. In each figure, rows are the seven biodiversity metrics and columns are the 10 m, 20 m, and 50 m quadrat scales.

## Data

- Input tables are `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`.
- Each panel uses all complete quadrat records for that spectral-biodiversity metric pair; no edge filtering is applied.
- Panels report sample size, Pearson r, linear-model R2, and the simple-regression F-test p-value.

## Strongest Absolute Correlations

| Scale | Spectral | Biodiversity | n | r | R2 | p-value |
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
| 50m | Spectral angle entropy | Phylogenetic Rao's Q |  80 | 0.340 | 0.115 | 0.002 |
| 20m | Std PCA alpha hull | Phylogenetic Rao's Q | 436 | 0.336 | 0.113 | 6.08e-13 |
| 50m | Raw PCA alpha hull | Phylogenetic Rao's Q |  74 | 0.332 | 0.111 | 0.004 |
| 50m | Raw PCA Rao's Q | Phylogenetic Rao's Q |  74 | 0.319 | 0.102 | 0.006 |
| 50m | Std PCA alpha hull | Faith's PD |  74 | 0.312 | 0.097 | 0.007 |
| 20m | Std PCA mean distance | Abundance-weighted Faith's PD | 436 | 0.301 | 0.091 | 1.36e-10 |
| 50m | Raw PCA mean distance | Phylogenetic Rao's Q |  74 | 0.292 | 0.085 | 0.012 |
| 20m | Std PCA alpha hull | Faith's PD | 436 | 0.285 | 0.081 | 1.40e-09 |
| 50m | Std PCA mean distance | Faith's PD |  74 | 0.282 | 0.080 | 0.015 |
| 20m | Std PCA mean distance | Phylogenetic Rao's Q | 436 | 0.274 | 0.075 | 6.16e-09 |
| 50m | Std PCA Rao's Q | Abundance-weighted Faith's PD |  74 | 0.272 | 0.074 | 0.019 |

## Complete Case Ranges

| scale | rows | min_complete_n | max_complete_n |
| --- | --- | --- | --- |
| 10m | 2000 | 1741 | 1903 |
| 20m |  500 |  436 |  485 |
| 50m |   80 |   74 |   80 |

## Figures

- `reports/figures/spectral_biodiversity_all_metrics/01_sa_entropy_vs_all_biodiversity_metrics_by_scale.png`
- `reports/figures/spectral_biodiversity_all_metrics/01_raw_pca_alpha_hull_vs_all_biodiversity_metrics_by_scale.png`
- `reports/figures/spectral_biodiversity_all_metrics/01_raw_pca_rao_q_vs_all_biodiversity_metrics_by_scale.png`
- `reports/figures/spectral_biodiversity_all_metrics/01_raw_pca_mean_distance_vs_all_biodiversity_metrics_by_scale.png`
- `reports/figures/spectral_biodiversity_all_metrics/01_std_pca_alpha_hull_vs_all_biodiversity_metrics_by_scale.png`
- `reports/figures/spectral_biodiversity_all_metrics/01_std_pca_rao_q_vs_all_biodiversity_metrics_by_scale.png`
- `reports/figures/spectral_biodiversity_all_metrics/01_std_pca_mean_distance_vs_all_biodiversity_metrics_by_scale.png`

## Tables

- `reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_all_metric_dataset.csv`
- `reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_all_metric_relationships.csv`
