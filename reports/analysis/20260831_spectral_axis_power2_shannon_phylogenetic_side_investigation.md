# Side Investigation: Spectral-Axis Power-2 Transform

Date: 2026-08-31

## Purpose

This side investigation repeats the spectral-biodiversity scatter comparison with a power-2 transformation applied only to the spectral metric axis. Biodiversity metrics remain on their original scale.

## Scope

- Spectral metrics: spectral angle entropy plus raw and standardized PCA alpha hull, Rao's Q, mean distance, convex hull, 3D hull volume, and 3D hull area.
- Biodiversity metrics: Faith's PD, phylogenetic Rao's Q, abundance-weighted Faith's PD, and Shannon diversity.
- Spatial scales: 10 m, 20 m, and 50 m.
- Each panel reports sample size, Pearson r, linear-model R2, and the simple-regression F-test p-value.

## Strongest Absolute Pearson Correlations

| Scale | Spectral | Biodiversity | n | r | R2 | p-value |
| --- | --- | --- | --- | --- | --- | --- |
| 50m | Std PCA alpha hull | Abundance-weighted Faith's PD |  74 | 0.490 | 0.240 | 9.22e-06 |
| 50m | Std PCA alpha hull | Phylogenetic Rao's Q |  74 | 0.481 | 0.232 | 1.42e-05 |
| 50m | Std PCA mean distance | Phylogenetic Rao's Q |  74 | 0.434 | 0.188 | 1.14e-04 |
| 50m | Raw PCA alpha hull | Abundance-weighted Faith's PD |  74 | 0.424 | 0.179 | 1.69e-04 |
| 50m | Std PCA mean distance | Abundance-weighted Faith's PD |  74 | 0.406 | 0.165 | 3.28e-04 |
| 50m | Raw PCA mean distance | Abundance-weighted Faith's PD |  74 | 0.388 | 0.150 | 6.45e-04 |
| 20m | Std PCA alpha hull | Abundance-weighted Faith's PD | 436 | 0.380 | 0.144 | 2.05e-16 |
| 50m | Raw PCA Rao's Q | Abundance-weighted Faith's PD |  74 | 0.355 | 0.126 | 0.002 |
| 50m | Raw PCA alpha hull | Phylogenetic Rao's Q |  74 | 0.346 | 0.119 | 0.003 |
| 50m | Raw PCA Rao's Q | Phylogenetic Rao's Q |  74 | 0.341 | 0.116 | 0.003 |
| 50m | Spectral angle entropy | Phylogenetic Rao's Q |  80 | 0.332 | 0.110 | 0.003 |
| 20m | Std PCA alpha hull | Phylogenetic Rao's Q | 436 | 0.328 | 0.108 | 2.15e-12 |
| 50m | Raw PCA 3D hull volume | Phylogenetic Rao's Q |  74 | 0.315 | 0.099 | 0.006 |
| 50m | Std PCA Rao's Q | Phylogenetic Rao's Q |  74 | 0.313 | 0.098 | 0.007 |
| 50m | Raw PCA mean distance | Phylogenetic Rao's Q |  74 | 0.305 | 0.093 | 0.008 |
| 50m | Std PCA alpha hull | Faith's PD |  74 | 0.292 | 0.085 | 0.012 |
| 20m | Std PCA alpha hull | Faith's PD | 436 | 0.266 | 0.071 | 1.71e-08 |
| 50m | Std PCA mean distance | Faith's PD |  74 | 0.263 | 0.069 | 0.024 |
| 20m | Std PCA mean distance | Abundance-weighted Faith's PD | 436 | 0.262 | 0.069 | 2.77e-08 |
| 50m | Raw PCA 3D hull volume | Faith's PD |  74 | 0.254 | 0.065 | 0.029 |

## Complete Case Ranges

| scale | rows | min_complete_n | max_complete_n |
| --- | --- | --- | --- |
| 10m | 2000 | 1741 | 1903 |
| 20m |  500 |  436 |  485 |
| 50m |   80 |   74 |   80 |

## Figures

- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_sa_entropy_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_raw_pca_alpha_hull_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_raw_pca_rao_q_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_raw_pca_mean_distance_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_raw_pca_convex_hull_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_raw_pca_3d_hull_volume_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_raw_pca_3d_hull_area_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_std_pca_alpha_hull_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_std_pca_rao_q_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_std_pca_mean_distance_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_std_pca_convex_hull_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_std_pca_3d_hull_volume_vs_shannon_phylogenetic_by_scale.png`
- `reports/figures/spectral_axis_power2_shannon_phylogenetic/spectral_power2_std_pca_3d_hull_area_vs_shannon_phylogenetic_by_scale.png`

## Table

- `reports/tables/spectral_axis_power2_shannon_phylogenetic/spectral_axis_power2_shannon_phylogenetic_summary.csv`
