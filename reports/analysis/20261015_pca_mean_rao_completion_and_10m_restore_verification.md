# PCA mean distance and Rao's Q completion; restored 10 m verification

Date: 2026-10-15

## Outcome

The restored 10 m spectra reproduce the saved vector-normalized PC1-PC3 mean
distance and Rao's Q for all 1,744 eligible reference quadrats. Retained-pixel
counts match exactly. All 312 prepared 25 m rasters now have mean distance and
Rao's Q in both original and vector-normalized PCA bases. No 25 m exclusions or
exclusion flags were added, per the user's instruction. No hulls were calculated.
33 m and 100 m were not processed.

## Metric coverage

Counts below are finite metric values, not full-grid inclusion counts.

| Scale | Original mean | Original Rao | Vector-normalized mean | Vector-normalized Rao | Grid cells without prepared rasters |
|---|---:|---:|---:|---:|---:|
| 10 m | 1,744 | 1,744 | 1,744 | 1,744 | 91 |
| 20 m | 436 | 436 | 436 | 436 | 15 |
| 25 m | 312 | 312 | 312 | 312 | 8 |
| 50 m | 74 | 74 | 74 | 74 | 0 |

10 m, 20 m, and 50 m retain their existing reference exclusions (165, 49, and 6
rows). These were not recalculated or expanded. 25 m has no exclusions applied;
the user will make those decisions manually. Existing 10 m, 20 m, and 50 m
component and canonical mean/Rao columns agree, with no missing finite metrics
among their non-excluded spectral-summary rows.

The eight absent 25 m raster IDs are `sub25_223`, `sub25_224`, `sub25_241`,
`sub25_242`, `sub25_243`, `sub25_261`, `sub25_281`, and `sub25_301`. They have
explicit `missing_raster` statuses and NA metrics in the 320-row 25 m outputs.

## Restored 10 m check

- Restored IDs: 1,909, matching the saved spectral-summary footprint exactly.
- Recomputed reference quadrats: 1,744; existing manually excluded rows: 165.
- Raster band names/order were checked against both saved PCA bases.
- All 1,744 mean, Rao, and retained-pixel-count comparisons passed.
- Maximum absolute mean-distance difference before CSV rounding:
  `4.9737991503207013e-14`.
- Maximum absolute Rao difference before CSV rounding:
  `3.751665644813329e-12`.
- Absolute tolerance `1e-10` plus relative tolerance `1e-10 * abs(reference)`
  was used per value. Differences are consistent with floating-point/CSV rounding.

This confirms that the restored rasters reproduce the previous metric workflow;
it is not a byte-for-byte comparison with the deleted files.

## Methods and validation

Reused the current project's masking, cleaning, vector normalization, and
projection functions from `spectral_heterogeneity_all_metrics.R`, while excluding
its automatic workflow invocation. Both saved PCA bases were loaded unchanged.
The shadow mask retains 563 nm values greater than `0.0305476`; invalid and
nonpositive-sum spectra are removed. All retained pixels are projected into
PC1-PC3. Mean distance uses Euclidean distance from the quadrat centroid. Rao's
Q uses equal pixel weights and squared Euclidean distance, calculated as
`2 * mean(squared_radius)`.

Independent pairwise-distance checks on 64 actual 25 m pixel scores in each basis
agree with the centroid Rao formula (maximum absolute difference
`2.8421709430404007e-14`). All 312 available 25 m rasters produced finite metrics
in both bases; no processing failures occurred. All 320 25 m boundary IDs are
represented exactly once in each component and canonical summary. Input raster
sizes and modification times remained unchanged during the run. MD5 hashes
confirm both saved PCA bases and all nine established-scale reference CSVs
remained unchanged. The installed `beepr` completion notification was invoked.

## Outputs

- [Canonical 25 m distance/Rao summary](../../Quad_Values/Spectral_diversitySHPs/spectral_mean_rao_25m_smooth_masked_5nm_summary.csv)
- [Original-PCA 25 m component](../../Quad_Values/25m_PCA_metrics_smooth_masked_5nm_summary.csv)
- [Vector-normalized 25 m component](../../Quad_Values/25m_standardized_PCA_metrics_smooth_masked_5nm_summary.csv)
- [Per-quadrat restored 10 m comparison](../tables/20261015_pca_mean_rao_completion/restored_10m_vector_normalized_comparison.csv)
- [Coverage chart data](../tables/20261015_pca_mean_rao_completion/metric_coverage.csv)
- [Validation checks](../tables/20261015_pca_mean_rao_completion/validation_checks.csv)
- [Independent pairwise Rao check](../tables/20261015_pca_mean_rao_completion/independent_pairwise_rao_check.csv)
- [Computation script](../../scripts/3_Analysis/complete_pca_mean_rao_and_verify_10m.R)

The new canonical 25 m table contains distance/Rao metrics only. No alpha-hull or
spectral-angle values are implied. A joined 25 m biodiversity/environment table
has not been created in this task. Alpha-hull work waits for the user's next step.
