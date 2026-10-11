# PCA and Spectral Angle Phylogenetic-Only Regression ExecPlan

Date: 2026-10-09

## Objective

Create a curiosity-check model layer that adds spectral angle entropy to the phylogenetic-only regression comparison already built for the six PCA spectral metrics.

## Requested Task

Fit phylogenetic-only OLS models for:

- Spectral angle entropy
- Three original PCA metrics
- Three vector-normalized PCA metrics

using Faith's PD, phylogenetic Rao's Q, and abundance-weighted Faith's PD as predictors.

## Proposed Changes

- Add `scripts/3_Analysis/pca_sa_phylogenetic_only_multiple_regression.R`.
- Write model-output CSVs under `reports/tables/pca_sa_phylogenetic_only_multiple_regression/`.
- Create `Documents/Paper/20261009_phylogenetic_only_pca_sa_multiple_regression_table.md`.
- Add task, analysis, validation, and lightweight test files.

## Validation Plan

- Confirm all 21 scale x spectral-metric models are estimable.
- Confirm coefficient, partial-R2, and VIF row counts match the expected model grid.
- Confirm the paper-facing Markdown table is created.
- Run a base-R output test.

## Risks

- Spectral angle entropy has different provenance and uncertainty behavior than the PCA metrics, so interpretation should distinguish it as an angular non-PCA comparison.
- OLS p-values remain cautious because prior diagnostics found spatial autocorrelation in several spectral-biodiversity relationships.
