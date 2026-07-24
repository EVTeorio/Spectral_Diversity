# ExecPlan: Direct Multiscale SV-Diversity Pairwise Analysis

## Objective

Create a scientifically defensible first-layer analysis that directly evaluates relationships between spectral variation and phylogenetic/species diversity across spatial scales using the current combined quadrat analysis tables.

## Requested Task

- Root the analysis in the project question: what is the relationship between spectral variation and phylogenetic/species diversity?
- Use current root-level combined quadrat tables for 10 m, 20 m, and 50 m scales.
- Treat `spec_spca_mean` and `spec_sa` as the two primary spectral variation measures.
- Independently pair each spectral variation measure with every phylogenetic and species diversity measure.
- Report Pearson `r`, `R2`, simple-regression F statistic, F-test p-value, slope, intercept, and Spearman rank correlation for each pair.
- Write a Markdown analysis report, CSV tables, and figures.

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
- Create a Markdown analysis report under `reports/analysis/`.
- Create task and validation reports documenting outputs and checks.
- Update `reports/project_state.md`, `reports/data_dictionary.md`, and `reports/directory_map.md` if output inventories change materially.

## Analysis Strategy

- Primary spectral variation responses: `spec_spca_mean` and `spec_sa`.
- Diversity predictors: `phy_faith`, `phy_rao`, `phy_afaith`, `sp_rich`, `sp_shannon`, `sp_simpson`, and `sp_even`.
- For each scale, spectral variation measure, and diversity measure, fit `SV_measure ~ diversity_measure`.
- Use complete cases per pair so known upstream spectral missingness is explicit rather than imputed.
- Keep documented 10 m and 20 m edge-quadrat exclusions; retain all 50 m quadrats unless a separate 50 m edge rule is documented.
- Summarize the strongest pairings by absolute Pearson `r` for each scale and SV measure.
- Treat environmental and spatial models as later sensitivity layers rather than the active first-layer analysis.

## Expected Outputs

- `reports/analysis/20260710_sv_diversity_pairwise_correlations.md`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_pairwise_correlations.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_top_pairings.csv`
- `reports/figures/multiscale_spectral_biodiversity/01_sv_diversity_pairwise_correlation_heatmap.png`
- `reports/figures/multiscale_spectral_biodiversity/02_sa_entropy_diversity_scatterplots.png`
- `reports/figures/multiscale_spectral_biodiversity/03_standardized_pca_distance_diversity_scatterplots.png`
- `reports/tasks/20260624_multiscale_spectral_biodiversity_analysis.md`
- `reports/validation/20260624_multiscale_spectral_biodiversity_analysis_validation.md`

## Validation Plan

- Confirm all three combined CSVs are readable.
- Confirm expected variables are present.
- Confirm the pairwise table has 42 rows: 3 scales x 2 SV measures x 7 diversity measures.
- Confirm no requested pair is missing Pearson `r` or F statistic values.
- Confirm all requested output CSVs and PNGs are created and non-empty.

## Risks

- The 50 m scale has only 80 quadrats, and standardized PCA distance has fewer complete 50 m rows because manual PCA exclusions are preserved.
- Known spectral missingness means sample sizes differ by response metric and scale.
- Spatial autocorrelation may remain, so simple-regression p-values should be treated as descriptive until a second-layer spatial sensitivity analysis is added.
- Diversity metrics are correlated with one another, so this first layer emphasizes independent pairwise relationships rather than causal separation among predictors.
