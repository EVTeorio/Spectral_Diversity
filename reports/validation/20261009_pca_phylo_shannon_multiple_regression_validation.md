# PCA Phylogenetic-Shannon Multiple Regression Validation

Date: 2026-10-09

## Validation Checks

| check | expected | observed | passed |
| --- | --- | --- | --- |
| required_columns_present | all required fields | all required fields | TRUE |
| full_model_row_count | 18 | 18 | TRUE |
| coefficient_row_count | 90 | 90 | TRUE |
| incremental_r2_row_count | 108 | 108 | TRUE |
| vif_row_count | 72 | 72 | TRUE |
| all_full_models_estimable | 18 | 18 | TRUE |
| finite_full_model_r2 | 18 | 18 | TRUE |
| complete_case_n_positive | 18 | 18 | TRUE |

## Complete-Case Counts

| Scale | Spectral | n | Total_rows |
| --- | --- | --- | --- |
| 10m | Original PCA alpha hull | 1741 | 2000 |
| 10m | Original PCA mean distance | 1744 | 2000 |
| 10m | Original PCA Rao's Q | 1744 | 2000 |
| 10m | Vector-normalized PCA alpha hull | 1744 | 2000 |
| 10m | Vector-normalized PCA mean distance | 1744 | 2000 |
| 10m | Vector-normalized PCA Rao's Q | 1744 | 2000 |
| 20m | Original PCA alpha hull |  436 |  500 |
| 20m | Original PCA mean distance |  436 |  500 |
| 20m | Original PCA Rao's Q |  436 |  500 |
| 20m | Vector-normalized PCA alpha hull |  436 |  500 |
| 20m | Vector-normalized PCA mean distance |  436 |  500 |
| 20m | Vector-normalized PCA Rao's Q |  436 |  500 |
| 50m | Original PCA alpha hull |   74 |   80 |
| 50m | Original PCA mean distance |   74 |   80 |
| 50m | Original PCA Rao's Q |   74 |   80 |
| 50m | Vector-normalized PCA alpha hull |   74 |   80 |
| 50m | Vector-normalized PCA mean distance |   74 |   80 |
| 50m | Vector-normalized PCA Rao's Q |   74 |   80 |

## Validation Outputs

- `reports/tables/pca_phylo_shannon_multiple_regression/pca_phylo_shannon_validation_checks.csv`
- `reports/tables/pca_phylo_shannon_multiple_regression/pca_phylo_shannon_complete_case_counts.csv`

## Notes

- All models use complete cases for the response and the four biodiversity predictors.
- The combined quadrat tables already exclude unavailable PCA values with `NA`, so manually excluded atmospheric/cloud quadrats are absent from complete-case fits for PCA-dependent metrics.
