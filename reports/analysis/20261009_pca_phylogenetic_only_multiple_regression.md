# PCA Phylogenetic-Only Multiple Regression

Date: 2026-10-09

## Purpose

This analysis fits scale-specific multiple regressions for the six PCA-based spectral heterogeneity metrics using only phylogenetic diversity predictors.

## Model Form

`spectral_metric ~ phy_faith + phy_rao + phy_afaith`

## Top Full Models

| Scale | Spectral | n | R2 | Adj_R2 | p |
| --- | --- | --- | --- | --- | --- |
| 50m | Vector-normalized PCA alpha hull |   74 | 0.289 | 0.259 | 2.44e-05 |
| 50m | Vector-normalized PCA mean distance |   74 | 0.222 | 0.189 | 4.99e-04 |
| 50m | Original PCA alpha hull |   74 | 0.177 | 0.141 | 0.003 |
| 20m | Vector-normalized PCA alpha hull |  436 | 0.152 | 0.146 | 2.25e-15 |
| 50m | Original PCA mean distance |   74 | 0.145 | 0.108 | 0.012 |
| 50m | Original PCA Rao's Q |   74 | 0.133 | 0.096 | 0.018 |
| 50m | Vector-normalized PCA Rao's Q |   74 | 0.126 | 0.088 | 0.024 |
| 20m | Vector-normalized PCA mean distance |  436 | 0.096 | 0.090 | 1.61e-09 |
| 20m | Original PCA alpha hull |  436 | 0.048 | 0.041 | 9.63e-05 |
| 10m | Vector-normalized PCA alpha hull | 1744 | 0.038 | 0.037 | 1.05e-14 |
| 20m | Original PCA mean distance |  436 | 0.027 | 0.020 | 0.008 |
| 20m | Vector-normalized PCA Rao's Q |  436 | 0.021 | 0.014 | 0.028 |

## Strongest Single-Predictor Partial Contributions

These rows show the additional R2 from adding one phylogenetic predictor after the other two phylogenetic predictors are already in the model.

| Scale | Spectral | Added | Base_R2 | Full_R2 | Incremental_R2 | p |
| --- | --- | --- | --- | --- | --- | --- |
| 50m | Original PCA alpha hull | Abundance-weighted Faith's PD | 0.118 | 0.177 | 0.059 | 0.028 |
| 50m | Original PCA mean distance | Abundance-weighted Faith's PD | 0.091 | 0.145 | 0.053 | 0.040 |
| 50m | Vector-normalized PCA alpha hull | Abundance-weighted Faith's PD | 0.238 | 0.289 | 0.051 | 0.028 |
| 20m | Vector-normalized PCA alpha hull | Abundance-weighted Faith's PD | 0.116 | 0.152 | 0.036 | 2.15e-05 |
| 50m | Vector-normalized PCA Rao's Q | Phylogenetic Rao's Q | 0.095 | 0.126 | 0.031 | 0.118 |
| 50m | Original PCA Rao's Q | Abundance-weighted Faith's PD | 0.106 | 0.133 | 0.028 | 0.139 |
| 50m | Vector-normalized PCA mean distance | Abundance-weighted Faith's PD | 0.197 | 0.222 | 0.026 | 0.134 |
| 50m | Vector-normalized PCA mean distance | Phylogenetic Rao's Q | 0.197 | 0.222 | 0.025 | 0.134 |
| 20m | Original PCA alpha hull | Abundance-weighted Faith's PD | 0.024 | 0.048 | 0.024 | 0.001 |
| 50m | Vector-normalized PCA alpha hull | Phylogenetic Rao's Q | 0.269 | 0.289 | 0.020 | 0.166 |
| 20m | Vector-normalized PCA mean distance | Abundance-weighted Faith's PD | 0.078 | 0.096 | 0.018 | 0.003 |
| 20m | Original PCA mean distance | Abundance-weighted Faith's PD | 0.014 | 0.027 | 0.013 | 0.018 |

## Paper-Facing Table

- `Documents/Paper/20261009_phylogenetic_only_pca_multiple_regression_table.md`

## Output Tables

- `reports/tables/pca_phylogenetic_only_multiple_regression/pca_phylogenetic_only_full_model_summary.csv`
- `reports/tables/pca_phylogenetic_only_multiple_regression/pca_phylogenetic_only_full_model_coefficients.csv`
- `reports/tables/pca_phylogenetic_only_multiple_regression/pca_phylogenetic_only_partial_r2.csv`
- `reports/tables/pca_phylogenetic_only_multiple_regression/pca_phylogenetic_only_predictor_vif.csv`
- `reports/tables/pca_phylogenetic_only_multiple_regression/pca_phylogenetic_only_top_full_models.csv`

## Interpretation Notes

- These models isolate the phylogenetic predictors by excluding Shannon diversity.
- The top model rankings are nearly identical to the phylogenetic-plus-Shannon layer, confirming that most explanatory power in the four-predictor model came from the phylogenetic block.
- As with the prior OLS layer, p-values should remain cautious because prior spatial diagnostics found residual spatial autocorrelation for several priority relationships.
