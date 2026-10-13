# All-Pixel Alpha-Hull Values: 50 m Alpha 1

Date: 2026-10-11

## Purpose

Calculate all-retained-pixel 2D alpha-hull spectral heterogeneity values for the current 50 m quadrats using alpha = 1.
The output table is designed to accumulate future alpha values and spatial scales in the same CSV.

## Scope

- Scale(s): 10m
- Alpha value(s): 0.5
- PCA basis: vector_normalized_pca
- Metric space: vector-normalized PCA PC1-PC2
- Manual atmospheric/cloud exclusions are retained as explicit skipped rows.

## Summary

| scale | basis | alpha_value | status | total_elapsed_seconds | rows |
| --- | --- | --- | --- | --- | --- |
| 10m | vector_normalized_pca | 0.5 | nonfinite_area | 11.57 |    3 |
| 10m | vector_normalized_pca | 0.5 | ok | 7003.27 | 1741 |
| 10m | vector_normalized_pca | 0.5 | skipped_manual_excluded | NA |  165 |

## Output Tables

- `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values.csv`
- `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_summary.csv`
- `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values_validation.csv`

## Notes

- `alpha_area_all_pixels` is the all-retained-pixel PC1-PC2 alpha-hull area.
- Duplicate PC1-PC2 coordinates are collapsed before alpha-hull geometry because duplicate points do not change the boundary.
- `previous_sampled_alpha_area` carries the earlier 10,000-point sampled vector-normalized PCA alpha-hull value for comparison.
