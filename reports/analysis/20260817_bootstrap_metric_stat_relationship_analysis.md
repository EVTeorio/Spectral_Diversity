# Bootstrap Metric Statistic Relationship Analysis

Date: 2026-08-17

## Purpose

This analysis plots spectral metric values against bootstrap variability statistics: SD, CV, and SE.

## Data Sources

- Final-pipeline SA entropy bootstrap diagnostics use `reports/tables/bootstrap_variation/bootstrap_variation_quadrat_diagnostics.csv` and summarize the 70-iteration SA entropy bootstrap.
- All-metric bootstrap-stat panels use the fixed-4,000-pixel condition from the sample-size sensitivity branch. This is intentionally labeled separately from the final pipeline because the current visible main pipeline does not contain per-iteration PCA metric CSVs.

## Strongest Relationships

| Source | Scale | Metric | Statistic | n | r | R2 | p-value |
| --- | --- | --- | --- | --- | --- | --- | --- |
| final_sa_entropy | 10m | Spectral angle entropy | Bootstrap CV | 1866 | -0.039 | 0.002 | 0.088 |
| final_sa_entropy | 10m | Spectral angle entropy | Bootstrap SD | 1866 | -0.002 | 0.000 | 0.934 |
| final_sa_entropy | 10m | Spectral angle entropy | Bootstrap SE | 1866 | -0.002 | 0.000 | 0.934 |
| final_sa_entropy | 20m | Spectral angle entropy | Bootstrap CV |  477 | -0.151 | 0.023 | 9.55e-04 |
| final_sa_entropy | 20m | Spectral angle entropy | Bootstrap SE |  477 | -0.120 | 0.014 | 0.009 |
| final_sa_entropy | 20m | Spectral angle entropy | Bootstrap SD |  477 | -0.120 | 0.014 | 0.009 |
| final_sa_entropy | 50m | Spectral angle entropy | Bootstrap SD |   79 | 0.171 | 0.029 | 0.131 |
| final_sa_entropy | 50m | Spectral angle entropy | Bootstrap SE |   79 | 0.171 | 0.029 | 0.131 |
| final_sa_entropy | 50m | Spectral angle entropy | Bootstrap CV |   79 | 0.134 | 0.018 | 0.238 |
| sample_size_fixed4000 | 10m | Std PCA Rao's Q | Bootstrap SD |   32 | 0.629 | 0.395 | 1.16e-04 |
| sample_size_fixed4000 | 10m | Std PCA Rao's Q | Bootstrap SE |   32 | 0.629 | 0.395 | 1.16e-04 |
| sample_size_fixed4000 | 10m | Raw PCA Rao's Q | Bootstrap SD |   32 | 0.496 | 0.246 | 0.004 |
| sample_size_fixed4000 | 20m | Std PCA Rao's Q | Bootstrap SE |   16 | 0.970 | 0.940 | 5.90e-10 |
| sample_size_fixed4000 | 20m | Std PCA Rao's Q | Bootstrap SD |   16 | 0.970 | 0.940 | 5.90e-10 |
| sample_size_fixed4000 | 20m | Raw PCA Rao's Q | Bootstrap SE |   16 | 0.940 | 0.884 | 6.31e-08 |
| sample_size_fixed4000 | 50m | Std PCA Rao's Q | Bootstrap SE |    8 | 0.943 | 0.889 | 4.49e-04 |
| sample_size_fixed4000 | 50m | Std PCA Rao's Q | Bootstrap SD |    8 | 0.943 | 0.889 | 4.49e-04 |
| sample_size_fixed4000 | 50m | Std PCA mean distance | Bootstrap SD |    8 | 0.830 | 0.688 | 0.011 |

## Figures

- `reports/figures/bootstrap_metric_relationships/01_final_sa_entropy_vs_bootstrap_stats_10m.png`
- `reports/figures/bootstrap_metric_relationships/01_final_sa_entropy_vs_bootstrap_stats_20m.png`
- `reports/figures/bootstrap_metric_relationships/01_final_sa_entropy_vs_bootstrap_stats_50m.png`
- `reports/figures/bootstrap_metric_relationships/02_fixed4000_metrics_vs_bootstrap_stats_10m.png`
- `reports/figures/bootstrap_metric_relationships/02_fixed4000_metrics_vs_bootstrap_stats_20m.png`
- `reports/figures/bootstrap_metric_relationships/02_fixed4000_metrics_vs_bootstrap_stats_50m.png`

## Output Tables

- `reports/tables/bootstrap_metric_relationships/final_sa_entropy_bootstrap_metric_stats.csv`
- `reports/tables/bootstrap_metric_relationships/fixed4000_sample_size_bootstrap_metric_stats.csv`
- `reports/tables/bootstrap_metric_relationships/bootstrap_metric_stat_relationships.csv`
