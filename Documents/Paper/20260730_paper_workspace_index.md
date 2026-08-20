# Paper Workspace Index

Created: 2026-07-30

This folder is the working location for manuscript-related outputs. New manuscript files should be saved in `Documents/Paper/` and include the creation date in the filename using `YYYYMMDD_...`.

## Current Files

- `20260730_working_outline.md`: first working manuscript outline for publication planning.
- `20260730_paper_workspace_index.md`: index and source map for the paper workspace.
- `20260803_writing_style_notes.md`: R-generated notes on writing style and outline framing from accessible prior Word/context documents.
- `20260803_markdown_context_digest.md`: R-generated digest of repository Markdown files read before creating the outline.
- `20260803_detailed_manuscript_outline.docx`: detailed Word outline focused on abundance-weighted Faith's PD and Shannon diversity against focal spectral heterogeneity metrics across scale.
- `20260803_remote_sensing_mdpi_author_requirements.md`: target-journal requirements note for MDPI Remote Sensing, including length, APC, section format, submission, references, and fit notes.
- `20260805_results_pca_mean_distance_afaith.docx`: Results-section draft focused on standardized PCA mean distance versus abundance-weighted Faith's PD across 10 m, 20 m, and 50 m scales.
- `20260805_results_spectral_rao_q_afaith.docx`: Results-section draft focused on spectral Rao's Q versus abundance-weighted Faith's PD across 10 m, 20 m, and 50 m scales.
- `20260805_discussion_pca_distance_vs_rao_q_scale.docx`: Discussion-section draft explaining why standardized PCA mean distance may outperform spectral Rao's Q, why correlations strengthen with scale, and how evolutionary/physiological constraints may limit the spectral-phylogenetic slope.
- `20260805_results_spectral_rao_q_afaith_centroid_clarified.docx`: clarified Rao's Q Results draft noting that the script computes equal-weight pairwise squared Euclidean Rao's Q via the equivalent centroid formula `2 * mean(squared_radius)`.
- `20260805_discussion_pca_distance_vs_rao_q_scale_centroid_clarified.docx`: clarified Discussion draft noting that both PCA mean distance and spectral Rao's Q use the quadrat centroid computationally, while Rao's Q represents equal-weight pairwise squared Euclidean dissimilarity.
- `SVH_Results.docx`: user-created Word results document for manual final edits. Agents may read this file for context but should not edit, overwrite, rename, or regenerate it unless the user explicitly changes this instruction.
- `20260819_final_research_direction.md`: focused final paper direction that prioritizes standardized PCA alpha hull, standardized PCA mean distance, and spectral angle entropy against phylogenetic diversity metrics, with species-diversity contrasts and sensitivity outputs assigned to main or supplemental roles.
- `20260819_final_research_narrative.md`: revised active manuscript narrative using the working title "Assessing the Spectral Variation Hypothesis Across Spatial Scales Using PCA-based Canopy Spectral Heterogeneity" and centering scale dependence, PCA metric behavior, biodiversity metric behavior, transformations, and correlation mechanisms.

## Naming Convention

- Use `YYYYMMDD_short_description.ext`.
- For Word outputs in this paper workflow, use `C:/Program Files/R/R-4.2.3/bin/Rscript.exe` unless a later instruction changes the tooling.
- Examples:
  - `20260730_working_outline.md`
  - `20260730_target_journal_notes.md`
  - `20260730_results_narrative_draft.md`
  - `20260730_figure_inventory.md`

## Current Manuscript Direction

Working manuscript focus: multiscale relationships between drone hyperspectral spectral heterogeneity and field-derived phylogenetic/species diversity in the Paint Rock Forest Dynamics Plot.

Primary spectral variation measures:

- `spec_spca_mean`: standardized PCA mean Euclidean distance after vector-normalizing spectra.
- `spec_sa`: spectral angle entropy.

Primary biodiversity emphasis:

- Phylogenetic Rao's Q
- Abundance-weighted Faith's PD
- Faith's PD

Supporting or contrast biodiversity measures:

- Species richness
- Shannon diversity
- Simpson diversity
- Species evenness

## Analysis-Ready Inputs

- `quadrat_analysis_10m.csv`
- `quadrat_analysis_20m.csv`
- `quadrat_analysis_50m.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_pairwise_correlations.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_top_pairings.csv`
- `reports/tables/edge_bootstrap_sensitivity/analysis_dataset_with_bootstrap_fields.csv`

For the species-diversity versus phylogenetic-diversity correlation pass, all available biodiversity quadrats are used because these metrics are not affected by the spectral atmospheric-distortion exclusions. `sv_diversity_analysis_dataset.csv` remains useful because it carries `edge_flag` for edge-highlight figures, but `primary_analysis` is not used to filter the biodiversity-only concordance tables.

## Key Analysis Reports

- `reports/analysis/20260710_sv_diversity_pairwise_correlations.md`
- `reports/analysis/20260725_edge_bootstrap_sensitivity_sv_diversity.md`
- `reports/analysis/20260817_species_phylogenetic_correlation_analysis.md`
- `reports/analysis/20260707_pca_loading_spectral_regions.md`
- `reports/analysis/20260707_standardized_pca_metric_sample_size_effects.md`
- `reports/analysis/20260618_bootstrap_variation_analysis.md`

## Methods Documentation Sources

- `RESEARCH_OBJECTIVES.md`
- `reports/project_state.md`
- `reports/data_dictionary.md`
- `reports/combined_quadrat_variable_guide.md`
- `reports/directory_map.md`
- `reports/validation/20260709_standardized_pca_diagnostics_and_sample_size_validation.md`
- `reports/validation/20260624_multiscale_spectral_biodiversity_analysis_validation.md`
- `reports/validation/20260624_combined_quadrat_analysis_tables_validation.md`

## Current Figure Sources

- `reports/figures/multiscale_spectral_biodiversity/`
- `reports/figures/species_phylogenetic_correlation/`
- `reports/figures/edge_bootstrap_sensitivity/`
- `reports/figures/pca_loading_spectral_regions/`
- `reports/figures/sample_size_effects/`
- `reports/figures/bootstrap_variation/`
- `reports/figures/spectral_heterogeneity/`
- `reports/figures/methods_tree_counts/`

## Current Table Sources

- `reports/tables/multiscale_spectral_biodiversity/`
- `reports/tables/species_phylogenetic_correlation/`
- `reports/tables/edge_bootstrap_sensitivity/`
- `reports/tables/pca_loading_spectral_regions/`
- `reports/tables/sample_size_effects/`
- `reports/tables/bootstrap_variation/`
- `reports/tables/methods_tree_counts/`

## 2026-08-11 Methods Draft Additions

- Created `Documents/Paper/20260811_methods_draft_spectral_phylogenetic_diversity.docx` as a methods-section draft for the current spectral-phylogenetic diversity manuscript.
- Methods draft was generated with R using `tools/create_methods_docx_and_tree_counts.R`.
- Style/context references read for the draft: `Documents/Methods-SampleT.docx`, `Documents/Spectral Variation Paper.docx`, and the read-only paper context in `Documents/Paper/SVH_Paper.docx`.
- The draft covers site context, quadrat construction, crown-overlap abundance estimation, species and phylogenetic diversity calculation, spectral preprocessing, spectral angle entropy, standardized PCA mean Euclidean distance, spectral Rao's Q, quadrat-level statistics, exclusion rules, and scale-specific correlation analysis.
- Tree-count support outputs were generated from the same R workflow. The all-plot tree pool contains 4,007 eligible tree crowns across 50 species after applying the plant-diversity filters used for phylogenetic diversity.
- Scale-specific included-tree tables count unique crowns intersecting focal complete analysis quadrats for `spec_spca_mean` versus abundance-weighted Faith's PD. These inclusion sets differ by quadrat scale because 10 m and 20 m use edge-excluded primary quadrats while 50 m does not use that edge rule.

## Immediate Next Step

Use `20260819_final_research_narrative.md` and root-level `Checklist.md` to update the manuscript outline, Results drafts, figure inventory, and table inventory. The active direction is now a scale-centered PCA-based Spectral Variation Hypothesis paper that asks how spectral metrics correlate, why they correlate with biodiversity metrics, and how scale changes those correlations. Complete transformation and spatial-autocorrelation checks before treating p-values as final inferential evidence.
