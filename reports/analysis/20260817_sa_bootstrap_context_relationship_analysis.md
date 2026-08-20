# SA Bootstrap Context Relationship Analysis

Date: 2026-08-17

## Purpose

This analysis plots the three final SA entropy bootstrap statistics against elevation, species count, crown-equivalent individual count, species composition type, overall brightness, and blue, green, red, and near-infrared brightness.

## Data

- SA bootstrap statistics come from `reports/tables/bootstrap_metric_relationships/final_sa_entropy_bootstrap_metric_stats.csv`.
- Elevation comes from the current `quadrat_analysis_*m.csv` tables.
- Species count and crown-equivalent individual count come from `quadrat_crown_equivalent_individual_totals.csv`.
- Species composition types come from the presence/absence cluster table `species_composition_clusters.csv`.
- Brightness variables come from `quadrat_pixel_brightness_summary.csv`.

## Strongest Context Relationships

| Type | Scale | Statistic | Predictor | n | Effect | R2_or_eta2 | p-value |
| --- | --- | --- | --- | --- | --- | --- | --- |
| continuous | 10m | SA bootstrap CV | NIR brightness | 1866 | 0.269 | 0.072 | 3.40e-32 |
| continuous | 10m | SA bootstrap SE | NIR brightness | 1866 | 0.262 | 0.069 | 1.28e-30 |
| continuous | 10m | SA bootstrap SD | NIR brightness | 1866 | 0.262 | 0.069 | 1.28e-30 |
| continuous | 10m | SA bootstrap CV | Red brightness | 1866 | 0.206 | 0.042 | 2.60e-19 |
| continuous | 10m | SA bootstrap SD | Red brightness | 1866 | 0.202 | 0.041 | 1.28e-18 |
| continuous | 10m | SA bootstrap SE | Red brightness | 1866 | 0.202 | 0.041 | 1.28e-18 |
| continuous | 50m | SA bootstrap CV | Species count |   79 | -0.195 | 0.038 | 0.086 |
| continuous | 50m | SA bootstrap SE | Individual count |   79 | -0.192 | 0.037 | 0.090 |
| continuous | 50m | SA bootstrap SD | Individual count |   79 | -0.192 | 0.037 | 0.090 |
| continuous | 50m | SA bootstrap SD | Species count |   79 | -0.192 | 0.037 | 0.090 |
| composition | 50m | SA bootstrap SE | Species composition type |   79 | ANOVA | 0.027 | 0.360 |
| composition | 50m | SA bootstrap SD | Species composition type |   79 | ANOVA | 0.027 | 0.360 |
| composition | 50m | SA bootstrap CV | Species composition type |   79 | ANOVA | 0.025 | 0.386 |
| composition | 20m | SA bootstrap CV | Species composition type |  477 | ANOVA | 0.016 | 0.187 |
| composition | 20m | SA bootstrap SE | Species composition type |  477 | ANOVA | 0.015 | 0.195 |
| composition | 20m | SA bootstrap SD | Species composition type |  477 | ANOVA | 0.015 | 0.195 |

## Figures

- `reports/figures/bootstrap_metric_relationships/03_sa_bootstrap_stats_vs_continuous_context_10m.png`
- `reports/figures/bootstrap_metric_relationships/03_sa_bootstrap_stats_vs_continuous_context_20m.png`
- `reports/figures/bootstrap_metric_relationships/03_sa_bootstrap_stats_vs_continuous_context_50m.png`
- `reports/figures/bootstrap_metric_relationships/04_sa_bootstrap_stats_vs_species_composition_10m.png`
- `reports/figures/bootstrap_metric_relationships/04_sa_bootstrap_stats_vs_species_composition_20m.png`
- `reports/figures/bootstrap_metric_relationships/04_sa_bootstrap_stats_vs_species_composition_50m.png`

## Tables

- `reports/tables/bootstrap_metric_relationships/sa_bootstrap_context_dataset.csv`
- `reports/tables/bootstrap_metric_relationships/sa_bootstrap_context_relationships.csv`
