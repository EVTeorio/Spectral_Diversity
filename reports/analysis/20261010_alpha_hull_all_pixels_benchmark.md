# All-Pixel Alpha-Hull Benchmark

Date: 2026-10-10

## Purpose

Benchmark R `alphahull` all-retained-pixel PC1-PC2 alpha-hull areas for 16 selected 10 m quadrats, 8 selected 20 m quadrats, and 4 selected 50 m quadrats.
The benchmark tests multiple alpha parameter values and compares the resulting areas with abundance-weighted Faith's PD (`phy_afaith`).

## Alpha Values

Tested alpha values in current selected subset: 0.1, 0.25, 0.5, 1, 1.5, 2

## Vector-Normalized PCA Correlations With Abundance-Weighted Faith's PD

| alpha | n | pearson_r | r_squared | pearson_p |
| --- | --- | --- | --- | --- |
| 0.10 | 28 | 0.887 | 0.786 | 0.000 |
| 0.25 | 28 | 0.882 | 0.777 | 0.000 |
| 0.50 | 28 | 0.859 | 0.739 | 0.000 |
| 1.00 | 28 | 0.821 | 0.674 | 0.000 |
| 1.50 | 28 | 0.806 | 0.650 | 0.000 |
| 2.00 | 28 | 0.804 | 0.646 | 0.000 |

## Status Counts

| status | rows |
| --- | --- |
| ok | 168 |

## Output Tables

- `reports/tables/alpha_hull_all_pixels_benchmark/selected_quadrats.csv`
- `reports/tables/alpha_hull_all_pixels_benchmark/r_alphahull_all_pixels_results.csv`
- `reports/tables/alpha_hull_all_pixels_benchmark/alpha_hull_phy_afaith_correlations.csv`
- `reports/tables/alpha_hull_all_pixels_benchmark/alpha_hull_runtime_summary.csv`

## Interpretation Notes

- These are benchmark-subset correlations, not final all-quadrat inference.
- Alpha-hull geometry uses all retained valid sunlit pixels projected into PC1-PC2 space; duplicate PC1-PC2 coordinates are collapsed because duplicates do not change the hull boundary.
- Scale-specific correlations are included in the CSV output. The 50 m subset has 4 quadrats, so that scale-specific estimate remains small-n and should be interpreted cautiously.
- Python benchmarking was intentionally deferred.
