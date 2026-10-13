# Multiscale all-pixel alpha hull spectral variation

Date: 2026-10-12

This note records the completed all-pixel 2D alpha hull production runs for alpha = 0.5, 1, and 1.5 across the 10 m, 20 m, and 50 m quadrat scales.

## Settings

- PCA basis: `vector_normalized_pca`
- PCA dimensions: PC1 and PC2
- Pixel rule: all valid pixels per included quadrat
- Alpha values: 0.5, 1, 1.5
- Spatial scales: 10 m, 20 m, 50 m

## Output tables

- Values: `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values.csv`
- Status summary: `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_status_summary.csv`
- Runtime summary: `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_runtime_summary.csv`
- Faith's PD correlations: `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_phy_afaith_correlations.csv`

## Status Summary

| scale | alpha | ok | nonfinite area | skipped manual exclusions |
|:---|---:|---:|---:|---:|
| 10m | 0.5 | 1741 | 3 | 165 |
| 10m | 1.0 | 1744 | 0 | 165 |
| 10m | 1.5 | 1742 | 2 | 165 |
| 20m | 0.5 | 434 | 2 | 49 |
| 20m | 1.0 | 436 | 0 | 49 |
| 20m | 1.5 | 436 | 0 | 49 |
| 50m | 0.5 | 74 | 0 | 6 |
| 50m | 1.0 | 74 | 0 | 6 |
| 50m | 1.5 | 74 | 0 | 6 |

## Runtime Summary

Times are summed only across `ok` rows with completed alpha hull areas.

| scale | alpha | total seconds | average seconds per computed quadrat | timed rows |
|:---|---:|---:|---:|---:|
| 10m | 0.5 | 7003.27 | 4.02 | 1741 |
| 10m | 1.0 | 2031.50 | 1.16 | 1744 |
| 10m | 1.5 | 1209.34 | 0.69 | 1742 |
| 20m | 0.5 | 4763.90 | 10.98 | 434 |
| 20m | 1.0 | 2487.47 | 5.71 | 436 |
| 20m | 1.5 | 2015.74 | 4.62 | 436 |
| 50m | 0.5 | 9905.36 | 133.86 | 74 |
| 50m | 1.0 | 9256.36 | 125.09 | 74 |
| 50m | 1.5 | 8625.87 | 116.57 | 74 |

## Correlation With Abundance-Weighted Faith's PD

Correlations use only rows with `status == "ok"` and finite `alpha_area_all_pixels` and `phy_afaith`.

| scale | alpha | n | Pearson r | R2 | Pearson p | Spearman rho | Spearman p |
|:---|---:|---:|---:|---:|---:|---:|---:|
| 10m | 0.5 | 1741 | 0.2219583 | 0.0492655 | 7.133321e-21 | 0.1439839 | 1.590355e-09 |
| 10m | 1.0 | 1744 | 0.1882354 | 0.0354326 | 2.258438e-15 | 0.1324393 | 2.835020e-08 |
| 10m | 1.5 | 1742 | 0.1419426 | 0.0201477 | 2.678233e-09 | 0.1170107 | 9.729521e-07 |
| 20m | 0.5 | 434 | 0.3748696 | 0.1405272 | 6.284746e-16 | 0.2586578 | 4.595551e-08 |
| 20m | 1.0 | 436 | 0.2539370 | 0.0644840 | 7.636895e-08 | 0.2112733 | 8.620216e-06 |
| 20m | 1.5 | 436 | 0.1736628 | 0.0301588 | 2.689574e-04 | 0.1614064 | 7.175226e-04 |
| 50m | 0.5 | 74 | 0.3439137 | 0.1182767 | 2.698200e-03 | 0.2929137 | 1.132021e-02 |
| 50m | 1.0 | 74 | 0.1901393 | 0.0361529 | 1.046677e-01 | 0.1578230 | 1.792833e-01 |
| 50m | 1.5 | 74 | 0.1182831 | 0.0139909 | 3.155155e-01 | 0.0780007 | 5.088869e-01 |

Across all three spatial scales, alpha = 0.5 currently has the strongest Pearson relationship with abundance-weighted Faith's PD among the completed all-pixel alpha hull runs.
