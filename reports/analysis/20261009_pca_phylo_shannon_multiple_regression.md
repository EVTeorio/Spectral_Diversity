# PCA Phylogenetic-Shannon Multiple Regression

Date: 2026-10-09

## Purpose

This analysis fits scale-specific multiple regressions for the six PCA-based spectral heterogeneity metrics. Each model uses Faith's PD, phylogenetic Rao's Q, abundance-weighted Faith's PD, and Shannon diversity as biodiversity predictors.

## Model Form

`spectral_metric ~ phy_faith + phy_rao + phy_afaith + sp_shannon`

Models are fit separately for each quadrat scale and each spectral metric. The analysis uses complete cases from the current combined quadrat tables and does not add edge, environmental, brightness, or spatial covariates.

## Top Full Models

| Scale | Spectral | n | R2 | Adj_R2 | p |
| --- | --- | --- | --- | --- | --- |
| 50m | Vector-normalized PCA alpha hull |   74 | 0.292 | 0.251 | 7.27e-05 |
| 50m | Vector-normalized PCA mean distance |   74 | 0.229 | 0.184 | 0.001 |
| 50m | Original PCA alpha hull |   74 | 0.202 | 0.156 | 0.003 |
| 20m | Vector-normalized PCA alpha hull |  436 | 0.174 | 0.166 | 5.07e-17 |
| 50m | Original PCA mean distance |   74 | 0.168 | 0.120 | 0.012 |
| 50m | Original PCA Rao's Q |   74 | 0.158 | 0.109 | 0.017 |
| 50m | Vector-normalized PCA Rao's Q |   74 | 0.133 | 0.083 | 0.040 |
| 20m | Vector-normalized PCA mean distance |  436 | 0.121 | 0.113 | 2.09e-11 |
| 20m | Original PCA alpha hull |  436 | 0.067 | 0.059 | 4.77e-06 |
| 10m | Vector-normalized PCA alpha hull | 1744 | 0.065 | 0.063 | 2.96e-24 |
| 10m | Vector-normalized PCA mean distance | 1744 | 0.044 | 0.042 | 3.03e-16 |
| 20m | Vector-normalized PCA Rao's Q |  436 | 0.043 | 0.035 | 7.17e-04 |

## Strongest Single-Predictor Partial Contributions

These rows show the additional R2 from adding one biodiversity predictor after the other three predictors are already in the model.

| Scale | Spectral | Added | Base_R2 | Full_R2 | Incremental_R2 | p |
| --- | --- | --- | --- | --- | --- | --- |
| 50m | Vector-normalized PCA Rao's Q | Phylogenetic Rao's Q | 0.095 | 0.133 | 0.038 | 0.085 |
| 50m | Vector-normalized PCA mean distance | Phylogenetic Rao's Q | 0.197 | 0.229 | 0.032 | 0.096 |
| 10m | Vector-normalized PCA alpha hull | Shannon diversity | 0.038 | 0.065 | 0.026 | 3.63e-12 |
| 50m | Original PCA alpha hull | Shannon diversity | 0.177 | 0.202 | 0.026 | 0.141 |
| 50m | Original PCA Rao's Q | Phylogenetic Rao's Q | 0.132 | 0.158 | 0.025 | 0.155 |
| 20m | Vector-normalized PCA mean distance | Shannon diversity | 0.096 | 0.121 | 0.025 | 5.17e-04 |
| 50m | Original PCA Rao's Q | Shannon diversity | 0.133 | 0.158 | 0.024 | 0.163 |
| 10m | Vector-normalized PCA mean distance | Shannon diversity | 0.020 | 0.044 | 0.024 | 4.83e-11 |
| 50m | Original PCA mean distance | Shannon diversity | 0.145 | 0.168 | 0.023 | 0.168 |
| 50m | Vector-normalized PCA alpha hull | Phylogenetic Rao's Q | 0.269 | 0.292 | 0.023 | 0.137 |
| 50m | Vector-normalized PCA alpha hull | Abundance-weighted Faith's PD | 0.269 | 0.292 | 0.023 | 0.139 |
| 20m | Vector-normalized PCA Rao's Q | Shannon diversity | 0.021 | 0.043 | 0.023 | 0.002 |

## Biodiversity Block Contributions

The phylogenetic block is Faith's PD, phylogenetic Rao's Q, and abundance-weighted Faith's PD fitted after Shannon diversity. The Shannon row is Shannon diversity fitted after the three phylogenetic predictors.

| Scale | Spectral | Added | Base_R2 | Full_R2 | Incremental_R2 | p |
| --- | --- | --- | --- | --- | --- | --- |
| 10m | Original PCA alpha hull | Phylogenetic block | 0.002 | 0.012 | 0.010 | 4.08e-04 |
| 10m | Original PCA alpha hull | Shannon diversity | 0.005 | 0.012 | 0.007 | 4.46e-04 |
| 10m | Original PCA mean distance | Phylogenetic block | 0.001 | 0.010 | 0.009 | 0.001 |
| 10m | Original PCA mean distance | Shannon diversity | 0.004 | 0.010 | 0.006 | 9.20e-04 |
| 10m | Original PCA Rao's Q | Phylogenetic block | 0.003 | 0.010 | 0.007 | 0.008 |
| 10m | Original PCA Rao's Q | Shannon diversity | 0.001 | 0.010 | 0.009 | 8.07e-05 |
| 10m | Vector-normalized PCA alpha hull | Phylogenetic block | 0.001 | 0.065 | 0.064 | 7.44e-25 |
| 10m | Vector-normalized PCA alpha hull | Shannon diversity | 0.038 | 0.065 | 0.026 | 3.63e-12 |
| 10m | Vector-normalized PCA mean distance | Phylogenetic block | 0.001 | 0.044 | 0.043 | 1.30e-16 |
| 10m | Vector-normalized PCA mean distance | Shannon diversity | 0.020 | 0.044 | 0.024 | 4.83e-11 |
| 10m | Vector-normalized PCA Rao's Q | Phylogenetic block | 0.004 | 0.016 | 0.012 | 9.83e-05 |
| 10m | Vector-normalized PCA Rao's Q | Shannon diversity | 0.002 | 0.016 | 0.013 | 1.18e-06 |
| 20m | Original PCA alpha hull | Phylogenetic block | 0.009 | 0.067 | 0.058 | 9.78e-06 |
| 20m | Original PCA alpha hull | Shannon diversity | 0.048 | 0.067 | 0.019 | 0.003 |
| 20m | Original PCA mean distance | Phylogenetic block | 0.007 | 0.043 | 0.036 | 0.001 |
| 20m | Original PCA mean distance | Shannon diversity | 0.027 | 0.043 | 0.016 | 0.008 |
| 20m | Original PCA Rao's Q | Phylogenetic block | 0.006 | 0.038 | 0.032 | 0.003 |
| 20m | Original PCA Rao's Q | Shannon diversity | 0.017 | 0.038 | 0.022 | 0.002 |
| 20m | Vector-normalized PCA alpha hull | Phylogenetic block | 0.001 | 0.174 | 0.173 | 1.05e-17 |
| 20m | Vector-normalized PCA alpha hull | Shannon diversity | 0.152 | 0.174 | 0.022 | 7.86e-04 |
| 20m | Vector-normalized PCA mean distance | Phylogenetic block | 0.001 | 0.121 | 0.120 | 5.94e-12 |
| 20m | Vector-normalized PCA mean distance | Shannon diversity | 0.096 | 0.121 | 0.025 | 5.17e-04 |
| 20m | Vector-normalized PCA Rao's Q | Phylogenetic block | 0.004 | 0.043 | 0.040 | 5.37e-04 |
| 20m | Vector-normalized PCA Rao's Q | Shannon diversity | 0.021 | 0.043 | 0.023 | 0.002 |
| 50m | Original PCA alpha hull | Phylogenetic block | 0.025 | 0.202 | 0.178 | 0.003 |
| 50m | Original PCA alpha hull | Shannon diversity | 0.177 | 0.202 | 0.026 | 0.141 |
| 50m | Original PCA mean distance | Phylogenetic block | 0.025 | 0.168 | 0.143 | 0.011 |
| 50m | Original PCA mean distance | Shannon diversity | 0.145 | 0.168 | 0.023 | 0.168 |
| 50m | Original PCA Rao's Q | Phylogenetic block | 0.012 | 0.158 | 0.145 | 0.011 |
| 50m | Original PCA Rao's Q | Shannon diversity | 0.133 | 0.158 | 0.024 | 0.163 |
| 50m | Vector-normalized PCA alpha hull | Phylogenetic block | 0.000 | 0.292 | 0.292 | 2.46e-05 |
| 50m | Vector-normalized PCA alpha hull | Shannon diversity | 0.289 | 0.292 | 0.003 | 0.570 |
| 50m | Vector-normalized PCA mean distance | Phylogenetic block | 0.000 | 0.229 | 0.229 | 4.25e-04 |
| 50m | Vector-normalized PCA mean distance | Shannon diversity | 0.222 | 0.229 | 0.007 | 0.445 |
| 50m | Vector-normalized PCA Rao's Q | Phylogenetic block | 0.001 | 0.133 | 0.132 | 0.020 |
| 50m | Vector-normalized PCA Rao's Q | Shannon diversity | 0.126 | 0.133 | 0.007 | 0.446 |

## Output Tables

- `reports/tables/pca_phylo_shannon_multiple_regression/pca_phylo_shannon_full_model_summary.csv`
- `reports/tables/pca_phylo_shannon_multiple_regression/pca_phylo_shannon_full_model_coefficients.csv`
- `reports/tables/pca_phylo_shannon_multiple_regression/pca_phylo_shannon_incremental_r2.csv`
- `reports/tables/pca_phylo_shannon_multiple_regression/pca_phylo_shannon_predictor_vif.csv`
- `reports/tables/pca_phylo_shannon_multiple_regression/pca_phylo_shannon_single_predictor_models.csv`
- `reports/tables/pca_phylo_shannon_multiple_regression/pca_phylo_shannon_top_full_models.csv`
- `reports/tables/pca_phylo_shannon_multiple_regression/pca_phylo_shannon_strongest_partial_predictors.csv`

## Interpretation Notes

- The vector-normalized PCA metrics remain the primary PCA evidence layer; original PCA metrics are retained as brightness-sensitive comparisons.
- These are ordinary least-squares models and should be interpreted cautiously because prior diagnostics found spatial autocorrelation in several spectral metrics, biodiversity metrics, and residuals.
- Predictor collinearity is expected because the biodiversity metrics are not independent. Use the standardized coefficients, incremental R2 rows, and VIF table together rather than interpreting raw coefficients alone.
- The 50 m models have the strongest R2 values but smaller complete-case sample sizes.
