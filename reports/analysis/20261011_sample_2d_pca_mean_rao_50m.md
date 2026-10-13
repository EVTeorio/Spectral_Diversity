# Sampled 2D PCA Mean Distance and Rao's Q

Date: 2026-10-11

## Purpose

Calculate vector-normalized PCA PC1-PC2 mean Euclidean distance and Rao's Q using the same 10,000-pixel sampling cap used by the sampled alpha-hull workflow.

## Scope

- Scale(s): 50m
- PCA basis: vector_normalized_pca
- Maximum sampled pixels per quadrat: 10000
- Metric space: vector-normalized PCA PC1-PC2
- Rao's Q uses equal pixel weights and squared Euclidean distance, computed as `2 * mean(squared radius)` around the sampled centroid.

## Status Summary

| status | rows | total_elapsed_seconds |
| --- | --- | --- |
| ok | 74 | 294.35 |
| skipped_manual_excluded |  6 |   0.00 |

## Correlations With Abundance-Weighted Faith's PD

| scale | basis | sample_max_points | metric_id | n | pearson_r | r_squared | pearson_p | spearman_r | spearman_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 50m | vector_normalized_pca | 10000 | sampled_2d_mean_distance | 74 | 0.417 | 0.174 | 0.000 | 0.427 | 0.000 |
| 50m | vector_normalized_pca | 10000 | sampled_2d_rao_q | 74 | 0.273 | 0.075 | 0.018 | 0.281 | 0.015 |

## Output Tables

- `reports/tables/sample_2d_pca_mean_rao_values/sample_2d_pca_mean_rao_values.csv`
- `reports/tables/sample_2d_pca_mean_rao_values/sample_2d_pca_mean_rao_values_summary.csv`
- `reports/tables/sample_2d_pca_mean_rao_values/sample_2d_pca_mean_rao_phy_afaith_correlations.csv`
- `reports/tables/sample_2d_pca_mean_rao_values/sample_2d_pca_mean_rao_validation.csv`
