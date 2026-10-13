# All-Pixel Alpha-Hull Values: 50 m Alpha 1

Date: 2026-10-11

## Purpose

Calculate all-retained-pixel 2D alpha-hull spectral heterogeneity values for the current 50 m quadrats using alpha = 1.
The output table is designed to accumulate future alpha values and spatial scales in the same CSV.

## Scope

- Scale(s): 50m
- Alpha value(s): 1
- PCA basis: vector_normalized_pca
- Metric space: vector-normalized PCA PC1-PC2
- Manual atmospheric/cloud exclusions are retained as explicit skipped rows.

## Summary

| scale | basis | alpha_value | status | total_elapsed_seconds | rows |
| --- | --- | --- | --- | --- | --- |
| 50m | vector_normalized_pca | 1 | ok | 9256.36 | 74 |
| 50m | vector_normalized_pca | 1 | skipped_manual_excluded | NA |  6 |

The completed 50 m alpha = 1 all-pixel calculations took 9,256.36 seconds total for the 74 included quadrats, or about 125.09 seconds per included quadrat on average.

## Correlation With Abundance-Weighted Faith's PD

Correlation was calculated across the 74 included 50 m quadrats with `status = ok`.

| scale | alpha | basis | n | Pearson r | R2 | Pearson p | Spearman r | Spearman p |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| 50m | 1 | vector_normalized_pca | 74 | 0.190 | 0.036 | 0.105 | 0.158 | 0.179 |

## Output Tables

- `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values.csv`
- `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_summary.csv`
- `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_validation.csv`
- `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_phy_afaith_correlations.csv`

## Notes

- `alpha_area_all_pixels` is the all-retained-pixel PC1-PC2 alpha-hull area.
- Duplicate PC1-PC2 coordinates are collapsed before alpha-hull geometry because duplicate points do not change the boundary.
- `previous_sampled_alpha_area` carries the earlier 10,000-point sampled vector-normalized PCA alpha-hull value for comparison.
