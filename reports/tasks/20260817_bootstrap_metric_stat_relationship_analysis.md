# Task Report: Bootstrap Metric Statistic Relationship Analysis

Date: 2026-08-17

## Objective

Create scale-specific scatterplot figures showing bootstrapped spectral metric values against metric SD, metric CV, and metric SE.

## Outputs Created

- `reports/figures/bootstrap_metric_relationships/01_final_sa_entropy_vs_bootstrap_stats_10m.png`
- `reports/figures/bootstrap_metric_relationships/01_final_sa_entropy_vs_bootstrap_stats_20m.png`
- `reports/figures/bootstrap_metric_relationships/01_final_sa_entropy_vs_bootstrap_stats_50m.png`
- `reports/figures/bootstrap_metric_relationships/02_fixed4000_metrics_vs_bootstrap_stats_10m.png`
- `reports/figures/bootstrap_metric_relationships/02_fixed4000_metrics_vs_bootstrap_stats_20m.png`
- `reports/figures/bootstrap_metric_relationships/02_fixed4000_metrics_vs_bootstrap_stats_50m.png`
- `reports/tables/bootstrap_metric_relationships/final_sa_entropy_bootstrap_metric_stats.csv`
- `reports/tables/bootstrap_metric_relationships/fixed4000_sample_size_bootstrap_metric_stats.csv`
- `reports/tables/bootstrap_metric_relationships/bootstrap_metric_stat_relationships.csv`

## Notes

- Final-pipeline bootstrap-stat figures are available for SA entropy.
- All-seven-metric figures use the fixed-4,000-pixel condition from the sample-size branch and should be interpreted as bootstrap stability diagnostics, not final-value uncertainty for PCA metrics.
