# ExecPlan: Species And Phylogenetic Diversity Correlation Prep

## Objective

Prepare the repository documentation for a focused correlation analysis comparing species-diversity measures with phylogenetic-diversity measures across 10 m, 20 m, and 50 m quadrats.

## Requested Task

Read the project Markdown context and update any documents necessary before beginning the analysis of correlations between the species and phylogenetic measures.

## Files Reviewed

- `RESEARCH_OBJECTIVES.md`
- `CODEX_AGENT_GUIDELINES.md`
- `reports/project_state.md`
- `reports/data_dictionary.md`
- `reports/directory_map.md`
- `reports/combined_quadrat_variable_guide.md`
- `reports/analysis/20260710_sv_diversity_pairwise_correlations.md`
- `reports/analysis/20260725_edge_bootstrap_sensitivity_sv_diversity.md`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset_column_guide.md`
- `Documents/Paper/20260730_paper_workspace_index.md`
- `Documents/Paper/20260803_markdown_context_digest.md`

## Relevant Inputs

- `quadrat_analysis_10m.csv`
- `quadrat_analysis_20m.csv`
- `quadrat_analysis_50m.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`

## Proposed Analysis Direction

Use the current analysis-ready quadrat tables to quantify direct relationships among:

- Species measures: `sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even`
- Phylogenetic measures: `phy_faith`, `phy_rao`, `phy_afaith`

Run the analysis within each spatial scale, preserve the current primary-analysis and edge-flag logic, and report Pearson and Spearman correlations with sample sizes and missingness. This analysis should be interpreted as a biodiversity-metric concordance layer, distinct from the already completed spectral-variation versus biodiversity correlation layer.

## Proposed Outputs

- Analysis report under `reports/analysis/`
- Correlation tables under `reports/tables/`
- Figures under `reports/figures/`
- Task and validation notes under `reports/tasks/` and `reports/validation/`

## Validation Plan

- Confirm input row counts and required columns for each scale.
- Confirm no unexpected duplicate `scale` + `quad_id` records.
- Confirm correlation tables use complete cases per species-phylogenetic pair.
- Confirm edge and primary-analysis filters match the current SV-diversity workflow.
- Confirm generated tables and figures exist and are non-empty.

## Risks

- Species and phylogenetic diversity measures are not independent because they are calculated from the same crown-overlap species matrix.
- Correlations may be inflated by shared richness or abundance structure.
- 10 m and 20 m edge-quadrat handling should be clearly separated from all-quadrat descriptive summaries.
- This preparatory update does not run the correlation analysis itself.
