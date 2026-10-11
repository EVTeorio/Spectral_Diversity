# PCA Phylogenetic-Only Multiple Regression ExecPlan

Date: 2026-10-09

## Objective

Create a phylogenetic-only multiple-regression layer for the six PCA-based spectral heterogeneity metrics and produce a paper-facing Markdown table.

## Requested Task

Model the three original PCA metrics and three vector-normalized PCA metrics using only:

- Faith's PD
- Phylogenetic Rao's Q
- Abundance-weighted Faith's PD

Exclude Shannon diversity from this model set.

## Proposed Changes

- Add `scripts/3_Analysis/pca_phylogenetic_only_multiple_regression.R`.
- Write model-output CSVs under `reports/tables/pca_phylogenetic_only_multiple_regression/`.
- Create `Documents/Paper/20261009_phylogenetic_only_pca_multiple_regression_table.md`.
- Add task, analysis, validation, and lightweight test files.

## Validation Plan

- Confirm all 18 scale x spectral-metric models are estimable.
- Confirm coefficient, partial-R2, and VIF row counts match the expected model grid.
- Confirm the paper-facing Markdown table is created.
- Run a base-R output test.

## Risks

- Phylogenetic predictors are correlated, so individual partial effects should be interpreted with VIF and incremental-R2 diagnostics.
- OLS p-values remain cautious because prior diagnostics found spatial autocorrelation in several spectral-biodiversity relationships.
