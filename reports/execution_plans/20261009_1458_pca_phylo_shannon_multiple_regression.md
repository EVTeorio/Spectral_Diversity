# PCA Phylogenetic-Shannon Multiple Regression ExecPlan

Date: 2026-10-09

## Objective

Add a reproducible multiple-regression analysis for the six PCA-based spectral heterogeneity metrics using phylogenetic diversity metrics plus Shannon diversity as predictors.

## Requested Task

Focus only on relationships between:

- Three original PCA metrics: alpha-hull area, mean Euclidean distance, and spectral Rao's Q.
- Three vector-normalized PCA metrics: alpha-hull area, mean Euclidean distance, and spectral Rao's Q.
- Biodiversity predictors: Faith's PD, phylogenetic Rao's Q, abundance-weighted Faith's PD, and Shannon diversity.

## Files To Review

- `CODEX_AGENT_GUIDELINES.md`
- `reports/project_state.md`
- `reports/data_dictionary.md`
- `reports/combined_quadrat_variable_guide.md`
- `scripts/3_Analysis/final_research_direction_analysis.R`
- Existing paper-facing PCA spectral-biodiversity tables under `Documents/Paper/`

## Proposed Changes

- Add `scripts/3_Analysis/pca_phylo_shannon_multiple_regression.R`.
- Write model outputs under `reports/tables/pca_phylo_shannon_multiple_regression/`.
- Write analysis, task, and validation notes under `reports/analysis/`, `reports/tasks/`, and `reports/validation/`.
- Update project-state and directory-map documentation.

## Expected Modifications

- Fit scale-specific full models of each spectral metric against all four biodiversity predictors.
- Fit reduced models for Shannon-only, phylogenetic-only, and each individual biodiversity predictor.
- Report coefficients, standardized coefficients, partial p-values, model R2, adjusted R2, AIC, BIC, nested-model incremental R2, VIF-style collinearity diagnostics, residual diagnostics, and top model summaries.

## Validation Plan

- Confirm required columns exist in all three combined quadrat tables.
- Confirm model rows cover all 3 scales x 6 spectral metrics.
- Confirm full-model coefficient rows cover all 3 scales x 6 spectral metrics x 4 predictors.
- Confirm complete-case counts match model inputs and key numeric outputs are finite where models are estimable.
- Run the R script end to end with R 4.2.3.

## Risks

- Predictors are correlated, so coefficient signs and partial p-values may differ from pairwise correlations.
- 50 m models have small sample sizes after PCA exclusions.
- Residual spatial autocorrelation is not solved by this OLS layer; results should remain cautious until a spatial sensitivity layer is selected.
