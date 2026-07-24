# Multiscale Spectral-Biodiversity Analysis Validation

Last updated: 2026-07-10

## Checks

- Expected pairwise rows: 42
- Observed pairwise rows: 42
- Pairwise row count matches expectation: TRUE
- Missing Pearson r values: 0
- Missing F statistics: 0
- Table CSV count: 3
- Figure PNG count: 3

## Coverage Summary

| Scale | Total quadrats | Primary quadrats | Complete SA entropy | Complete standardized PCA distance | SA all pixels sampled | Std PCA all pixels sampled |
| --- | --- | --- | --- | --- | --- | --- |
| 10m | 2000 | 1656 | 1587 | 1440 | 108 | 1440 |
| 20m |  500 |  414 |  405 |  360 |  12 |  360 |
| 50m |   80 |   80 |   80 |   74 |   1 |   74 |

## Figure Files

- `reports/figures/multiscale_spectral_biodiversity/01_sv_diversity_pairwise_correlation_heatmap.png`
- `reports/figures/multiscale_spectral_biodiversity/02_sa_entropy_diversity_scatterplots.png`
- `reports/figures/multiscale_spectral_biodiversity/03_standardized_pca_distance_diversity_scatterplots.png`

## Table Files

- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_pairwise_correlations.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_top_pairings.csv`

## Notes

- This validation checks that the requested pairwise correlation layer is complete.
- Spatial and environmental sensitivity models remain future second-layer analyses, not part of this direct correlation output.
