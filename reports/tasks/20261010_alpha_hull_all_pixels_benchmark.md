# Task Report: All-Pixel Alpha-Hull Benchmark

Date: 2026-10-10

## Actions

- Added `scripts/3_Analysis/alpha_hull_all_pixels_benchmark.R`.
- Selected representative included quadrats across 10 m, 20 m, and 50 m scales: 10m n=16, 20m n=8, 50m n=4.
- Recomputed all-retained-pixel PC1-PC2 alpha-hull areas with R `alphahull` across multiple alpha values.
- Wrote benchmark results, runtime summaries, correlations, and an analysis report under `reports/`.

## Notes

- The benchmark was run for the primary vector-normalized PCA basis.
- The requested Python comparison was deferred for a later step.
