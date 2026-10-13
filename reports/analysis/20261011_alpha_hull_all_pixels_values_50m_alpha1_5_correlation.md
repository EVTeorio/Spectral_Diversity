# 50 m all-pixel alpha hull spectral variation, alpha 1.5

Date: 2026-10-11

This note records the production 50 m all-pixel alpha hull spectral variation run for alpha = 1.5, using the same workflow previously used for alpha = 1.

## Inputs and settings

- Spatial scale: 50 m quadrats
- PCA basis: `vector_normalized_pca`
- PCA dimensions: PC1 and PC2
- Alpha value: 1.5
- Pixel rule: all valid pixels per included quadrat
- Manual atmospheric/cloud exclusions retained as skipped rows: `sub50_80`, `sub50_79`, `sub50_71`, `sub50_70`, `sub50_62`, `sub50_53`

## Output tables

- Values: `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values.csv`
- Summary: `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_summary.csv`
- Validation: `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_validation.csv`
- Correlations: `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_phy_afaith_correlations.csv`

The values table now contains both alpha = 1 and alpha = 1.5 records for the 50 m `vector_normalized_pca` basis.

## Run summary

| alpha | scale | computed rows | skipped manual exclusions | total alphahull seconds | average seconds per computed quadrat |
|---:|:---|---:|---:|---:|---:|
| 1.5 | 50m | 74 | 6 | 8625.87 | 116.57 |

## Correlation with abundance-weighted Faith's PD

Correlations use only computed 50 m rows with finite alpha hull area and finite `phy_afaith`.

| alpha | n | Pearson r | R2 | Pearson p | Spearman rho | Spearman p |
|---:|---:|---:|---:|---:|---:|---:|
| 1.0 | 74 | 0.1901393 | 0.0361529 | 0.1046677 | 0.1578230 | 0.1792833 |
| 1.5 | 74 | 0.1182831 | 0.0139909 | 0.3155155 | 0.0780007 | 0.5088869 |

Relative to alpha = 1, alpha = 1.5 produced a weaker relationship with abundance-weighted Faith's PD for the 50 m all-pixel run.
