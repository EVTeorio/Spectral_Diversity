# ExecPlan: Multiscale Spectral-Biodiversity Analysis Reports

## Objective

Create a scientifically defensible downstream analysis that evaluates relationships between hyperspectral spectral heterogeneity, species/phylogenetic diversity, environmental covariates, and spatial scale using the current combined quadrat analysis tables.

## Requested Task

- Lean on `scripts/3_Analysis/LLM.R`, `scripts/3_Analysis/Analysis_PDF.R`, and `RESEARCH_OBJECTIVES.md`.
- Analyze current root-level combined quadrat tables for 10 m, 20 m, and 50 m scales.
- Generate user-friendly PDF reports summarizing findings.
- Save PDFs under `Documents/PDFs/`.
- Save every figure used in the PDFs under a figures folder, consistent with prior workflows.

## Files To Review

- `RESEARCH_OBJECTIVES.md`
- `reports/project_state.md`
- `reports/data_dictionary.md`
- `reports/combined_quadrat_variable_guide.md`
- `reports/analysis/20260618_bootstrap_variation_analysis.md`
- `scripts/3_Analysis/LLM.R`
- `scripts/3_Analysis/Analysis_PDF.R`
- `scripts/3_Analysis/combine_quadrat_analysis_tables.R`

## Relevant Inputs

- `quadrat_analysis_10m.csv`
- `quadrat_analysis_20m.csv`
- `quadrat_analysis_50m.csv`
- Optional spectral bootstrap QC companion files under `Quad_Values/*_SA_entropy_smooth_masked_5nm_summary.csv`

## Proposed Changes

- Add a reproducible R workflow under `scripts/3_Analysis/`.
- Create an analysis output folder under `reports/figures/`.
- Create analysis tables under `reports/tables/`.
- Create PDF reports under `Documents/PDFs/`.
- Create task and validation reports documenting outputs and checks.
- Update `reports/project_state.md`, `reports/data_dictionary.md`, and `reports/directory_map.md` if output inventories change materially.

## Analysis Strategy

- Treat `spec_sa` as the primary spectral heterogeneity response.
- Treat PCA-space metrics (`spec_pca_mean`, `spec_rao_q`, `spec_alpha`) as secondary response/sensitivity metrics.
- Evaluate primary biodiversity predictors from the research objectives: `phy_faith`, `phy_afaith`, and `sp_shannon`.
- Evaluate environmental covariates: `env_elev`, `env_tri5`, `env_tri11`, and `env_tri21`.
- Use per-scale complete-case datasets for each response so known upstream spectral missingness is explicit rather than imputed.
- Fit interpretable model sets per scale:
  - null model;
  - single biodiversity models;
  - environmental model;
  - biodiversity plus environment models;
  - interaction with elevation for key biodiversity predictors where sample size supports it.
- Use standardized predictors for coefficient comparability.
- Report coefficients, standard errors, confidence intervals, p-values, R2/adjusted R2, AIC, and model comparisons.
- Use k-nearest-neighbor Moran's I diagnostics on model residuals to assess remaining spatial autocorrelation.
- Include scale comparisons and visual summaries.

## Expected Outputs

- `Documents/PDFs/spectral_biodiversity_multiscale_findings.pdf`
- `Documents/PDFs/spectral_biodiversity_model_appendix.pdf`
- Figure PNGs under `reports/figures/multiscale_spectral_biodiversity/`
- Analysis tables under `reports/tables/multiscale_spectral_biodiversity/`
- `reports/tasks/20260624_multiscale_spectral_biodiversity_analysis.md`
- `reports/validation/20260624_multiscale_spectral_biodiversity_analysis_validation.md`

## Validation Plan

- Confirm all three combined CSVs are readable.
- Confirm expected variables are present.
- Confirm PDF files are created and non-empty.
- Confirm all figures referenced in the PDFs are saved as standalone PNG files.
- Render PDFs to PNG pages with Poppler where available and inspect representative pages for legibility, alignment, and clipping.
- Run the completion beep after successful workflow execution.

## Risks

- The 50 m scale has only 80 quadrats, so multivariable and interaction models should be interpreted cautiously.
- Known spectral missingness means sample sizes differ by response metric and scale.
- Spatial autocorrelation may remain in residuals; if so, OLS model results are still useful descriptively, but inference should be interpreted with caution.
- Strong collinearity among biodiversity metrics can make multivariable biodiversity models unstable; single-predictor and grouped model comparisons will be emphasized.
