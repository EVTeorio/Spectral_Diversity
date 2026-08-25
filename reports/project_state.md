# Project State

Last updated: 2026-08-19

## Current Objective

Build a reproducible, well-documented research workflow and publication manuscript framework for evaluating relationships between hyperspectral spectral heterogeneity and biodiversity/phylogenetic diversity metrics in the Paint Rock Forest Dynamics Plot.

Immediate analysis objective: quantify concordance among species-diversity and phylogenetic-diversity measures before extending interpretation of the spectral-diversity relationships.

## Active Research Questions

- Does spectral heterogeneity increase with Faith's phylogenetic diversity?
- Do relationships differ across 10 m, 20 m, and 50 m spatial grains?
- Do environmental covariates such as elevation, slope, aspect, or topographic exposure explain additional variation in spectral heterogeneity?
- Which spectral heterogeneity metric is most useful for the biodiversity relationship: spectral angle entropy, band entropy, PCA-space metrics, or another metric?
- How strongly do species-diversity measures (`sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even`) correlate with phylogenetic-diversity measures (`phy_faith`, `phy_rao`, `phy_afaith`) within each quadrat scale?
- Do species-phylogenetic correlations change across 10 m, 20 m, and 50 m grains, and do they help explain why phylogenetic measures currently track spectral variation more strongly than species-diversity measures?

## Completed Work

- Confirmed `CODEX_AGENT_GUIDELINES.md` and `RESEARCH_OBJECTIVES.md` as governing/context documents.
- Updated documentation to reflect the user's independent testing confirmation that partitioned quadrat spectra in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m` are the spectral inputs to use moving forward.
- Updated agent guidance to require review of existing `scripts/` R code style and loop patterns before writing or restructuring R scripts.
- Applied Savitzky-Golay smoothing to confirmed partitioned quadrat spectra at 10 m, 20 m, and 50 m scales.
- Created current smoothed spectra output folders:
  - `Quad_Spectra/10m_smooth` with 1,909 rasters
  - `Quad_Spectra/20m_smooth` with 485 rasters
  - `Quad_Spectra/50m_smooth` with 80 rasters
- Updated `scripts/1_Data Processing/smoothing.R` to process all confirmed scales, create `_smooth` output folders, skip completed outputs, log progress, and use scale-aware `terra::app()` core counts.
- Resampled the current smoothed spectra to 5 nm spacing.
- Created current smoothed 5 nm spectra output folders:
  - `Quad_Spectra/10m_smooth_5nm` with 1,909 rasters
  - `Quad_Spectra/20m_smooth_5nm` with 485 rasters
  - `Quad_Spectra/50m_smooth_5nm` with 80 rasters
- Updated `scripts/1_Data Processing/resampling.R` to process all current smoothed scales, create `_smooth_5nm` output folders, skip completed outputs, and log progress.
- Reviewed all current Markdown project-state, context, task, execution-plan, and validation notes before beginning the spectral heterogeneity workflow update.
- Confirmed that the active pickup point for the next analysis stage is the current smoothed and 5 nm resampled quadrat spectra in `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, and `Quad_Spectra/50m_smooth_5nm`.
- Refactored `scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R` so it can calculate per-quadrat spectral angle entropy for 10 m, 20 m, and 50 m current `_smooth_5nm` spectra.
- Added a guarded all-pixel entropy path for small enough masked quadrats and automatic fallback to `boot_mean` from bootstrap subsamples when all-pixel pair counts are too large.
- Added sampled-pair entropy inside large bootstrap replicates to make the workflow computationally feasible while preserving bootstrap subsampling of pixels.
- Preserved and documented the current shadow masking threshold in the spectral heterogeneity workflow: `563 nm`, threshold `0.0305476`, retaining sunlit pixels greater than the threshold under the current direction setting.
- Added `scripts/2_Indices Creation/Spectral_diversity/run_sa_entropy_scale.R` to run the spectral heterogeneity workflow one scale at a time.
- Added `tests/test_sa_entropy_bootstrapping.R` with base R helper checks for identifier parsing, spectral normalization, pairwise angle generation, sampled-pair angle generation, and entropy calculation.
- Added `tests/smoke_sa_entropy_bootstrapping.R` for a one-raster-per-scale smoke test with reduced bootstrap settings.
- Confirmed R package availability through the R 4.2 user library when Codex is allowed to access `C:/Users/PaintRock/AppData/Local/R/win-library/4.2`; `terra` 1.7.71 and `beepr` are available there.
- Ran the spectral heterogeneity workflow for all three scales from current `_smooth_5nm` spectra.
- Created current spectral heterogeneity outputs:
  - `Quad_Values/10m_SA_entropy_smooth_masked_5nm_summary.csv`
  - `Quad_Values/20m_SA_entropy_smooth_masked_5nm_summary.csv`
  - `Quad_Values/50m_SA_entropy_smooth_masked_5nm_summary.csv`
  - matching long and wide bootstrap CSVs for 10 m, 20 m, and 50 m
  - `Quad_Values/Spectral_diversitySHPs/SA_entropy_10m_smooth_masked_5nm.shp`
  - `Quad_Values/Spectral_diversitySHPs/SA_entropy_20m_smooth_masked_5nm.shp`
  - `Quad_Values/Spectral_diversitySHPs/SA_entropy_50m_smooth_masked_5nm.shp`
- Added `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R` to create all-scale spectral heterogeneity outputs with four primary measures:
  - existing spectral angle entropy
  - PCA-space mean Euclidean distance from the quadrat centroid
  - alpha-hull area in global PC1-PC2 space
  - Rao's Q in global PC1-PC3 space using equal pixel weights and squared Euclidean distance
- Rebuilt the global PCA basis to avoid nested multiscale resampling of the same footprint. The current PCA basis samples only current 10 m smoothed 5 nm rasters after the 563 nm shadow mask is applied.
- Current PCA-dependent outputs supersede and replace the earlier 2026-06-22 PCA-dependent outputs and the later multiscale-nested PCA basis. Those older PCA values should be disregarded for downstream interpretation.
- Created the raw global PCA basis from 1,909 current 10 m rasters with up to 450 retained illuminated pixels per raster, producing 854,767 sampled pixels and 121 spectral bands. The sample count is below 1,909 x 450 because some 10 m rasters have fewer than 450 valid illuminated pixels after masking.
- Created a second PCA basis with the same 10 m sampling logic after vector-normalizing each retained spectrum; standardized PCA outputs use the `standardized_PCA...` naming convention.
- Created PCA variance diagnostic outputs:
  - `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`
  - `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm_variance_explained.csv`
  - `reports/figures/spectral_heterogeneity/global_pca_smooth_masked_5nm_variance_explained.png`
  - `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds`
  - `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm_variance_explained.csv`
  - `reports/figures/spectral_heterogeneity/standardized_PCA_global_pca_smooth_masked_5nm_variance_explained.png`
- Raw global PCA variance explained: PC1 = 66.23%, PC2 = 20.72%, PC3 = 3.76%, cumulative PC1-PC3 = 90.71%.
- Standardized PCA variance explained: PC1 = 45.46%, PC2 = 22.54%, PC3 = 6.80%, cumulative PC1-PC3 = 74.80%.
- Created all-scale combined spectral heterogeneity CSVs and shapefiles:
  - `Quad_Values/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
  - `Quad_Values/20m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
  - `Quad_Values/50m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
  - `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_10m_smooth_masked_5nm_summary.csv`
  - `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_20m_smooth_masked_5nm_summary.csv`
  - `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_50m_smooth_masked_5nm_summary.csv`
  - `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_10m_smooth_masked_5nm.shp`
  - `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_20m_smooth_masked_5nm.shp`
  - `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_50m_smooth_masked_5nm.shp`
- Created PCA-only metric CSVs:
  - `Quad_Values/10m_PCA_metrics_smooth_masked_5nm_summary.csv`
  - `Quad_Values/20m_PCA_metrics_smooth_masked_5nm_summary.csv`
  - `Quad_Values/50m_PCA_metrics_smooth_masked_5nm_summary.csv`
- Created standardized PCA-only metric CSVs:
  - `Quad_Values/10m_standardized_PCA_metrics_smooth_masked_5nm_summary.csv`
  - `Quad_Values/20m_standardized_PCA_metrics_smooth_masked_5nm_summary.csv`
  - `Quad_Values/50m_standardized_PCA_metrics_smooth_masked_5nm_summary.csv`
- Added `reports/validation/20260622_spectral_heterogeneity_all_metrics_validation.md`.
- Created a bootstrap variation QC analysis comparing within-quadrat bootstrap variation to between-quadrat spectral heterogeneity variation.
- Added `scripts/3_Analysis/bootstrap_variation_analysis.R`.
- Created `reports/analysis/20260618_bootstrap_variation_analysis.md` with seven figures and three diagnostic tables under `reports/figures/bootstrap_variation/` and `reports/tables/bootstrap_variation/`.
- Added `reports/figures/bootstrap_variation/spectral_entropy_histograms_by_scale.png` to visualize spectral heterogeneity distribution shape and skewness by scale.
- Bootstrap variation analysis conclusion: `spectral_entropy` / `boot_mean` values are usable as primary heterogeneity means, but `boot_sd`, `boot_cv`, bootstrap-mean CI fields, and `method` should be carried forward as quality-control fields and high-CV or wide-CI quadrats should be used in sensitivity checks.
- Added `scripts/3_Analysis/sa_entropy_sample_size_effects.R` to compare SA entropy means across fixed and percentage-based bootstrap pixel sample sizes for 32 selected 10 m quadrats, 16 selected 20 m quadrats, and 8 selected 50 m quadrats. The script retains the original six pilot quadrats; adds fixed 2,000 and 3,000 pixel rules for 10 m and 20 m, fixed 6,000 and 8,000 pixel rules for 50 m, and percent-based 1%, 2%, and 3% rules; handles full-pixel conditions deterministically so 100% retained-pixel rows have zero artificial bootstrap variation; and writes outputs to `reports/analysis/20260703_sa_entropy_sample_size_effects.md`, `reports/tables/sample_size_effects/sa_entropy/`, and `reports/figures/sample_size_effects/sa_entropy/`.
- Added `scripts/3_Analysis/spectral_metric_sample_size_effects.R` to extend the same sample-size sensitivity design to PCA mean distance, spectral Rao's Q, and alpha-hull area for both regular PCA and vector-normalized standardized PCA. The script reuses the finalized SA entropy quadrat/sample design, calculates full-pixel conditions deterministically, writes regular metric-specific tables under `reports/tables/sample_size_effects/{pca_mean_distance,spectral_rao_q,alpha_hull_area}/`, writes standardized metric-specific tables under `reports/tables/sample_size_effects/standardized_PCA_*`, writes matching figure folders, draws 95% CI error bars on the mean-by-sample-size figures, and summarizes the outputs in `reports/analysis/20260704_pca_metric_sample_size_effects.md` and `reports/analysis/20260707_standardized_pca_metric_sample_size_effects.md`.
- Reviewed current Markdown context before the next task and refreshed `reports/directory_map.md` to match the live working tree. As of 2026-07-06, the live derived-output directory is `Quad_Values/`.
- Updated scripts, reports, execution plans, validation notes, task reports, sample-size design CSVs, governance docs, and `combined_quadrat_variable_guide.docx` so `Quad_Values/` is the canonical derived-output directory. Added `reports/tasks/20260706_quad_values_path_update.md` and `reports/validation/20260706_quad_values_path_update_validation.md`.
- Added `scripts/3_Analysis/pc1_mean_reflectance_correlation_50m.R` to test the relationship between per-pixel mean reflectance and PC1, PC2, and PC3 for both regular PCA and vector-normalized standardized PCA. The analysis samples the same 250 valid sunlit pixels from each current 50 m quadrat for each PCA basis, producing 20,000 sampled pixels per basis from 80 quadrats, and writes combined and basis-specific tables, figures, an analysis report, a task note, and validation notes under `reports/`.
- The 50 m mean-reflectance diagnostic now uses standardized PCA as the key brightness-reduced interpretation. Under standardized PCA, PC1 has weak pixel-level association with mean reflectance (r = -0.1178, R2 = 0.0139) and a modest quadrat-mean association (r = -0.3482, R2 = 0.1213); PC2 has a stronger pixel-level association (r = -0.3312, R2 = 0.1097) but almost no quadrat-mean association (r = 0.0598, R2 = 0.0036). Regular PCA results are retained as supporting comparison and show the expected brightness-dominated PC1 relationship.
- Added `scripts/3_Analysis/pca_loading_spectral_region_analysis.R` to identify where PC1 and PC2 loadings are not uniform across wavelength for both regular and standardized PCA. The key finding is stated from standardized PCA: standardized PC1 has its strongest nonuniform loading region at 843-903 nm, standardized PC2 has concentrated structure at 753-813 nm, and the brightness-reduced standardized PC2 region to carry forward is 753-778 nm. Regular PCA loading regions are retained only as supporting comparison.
- Refactored the plant-diversity workflow into `scripts/2_Indices Creation/Plant_diversity/plant_diversity_all_scales.R`.
- Updated the older plant-diversity entry scripts (`sp_weighted_matrix.R`, `species_diversity.R`, and `phylogenetic_diversity.R`) to call the all-scale workflow.
- Created all-scale plant diversity outputs under `Quad_Values/Diversity_SHPs/`:
  - `plant_diversity_10m.csv` and `plant_diversity_10m.shp`
  - `plant_diversity_20m.csv` and `plant_diversity_20m.shp`
  - `plant_diversity_50m.csv` and `plant_diversity_50m.shp`
- Confirmed plant-diversity `quad_id` values align with current spectral diversity IDs: 10 m uses `sub_id` values such as `0_a`, 20 m uses `Name`, and 50 m uses `Name` values such as `sub50_1`.
- Added `scripts/2_Indices Creation/Enviro_Variables/environmental_variables_all_scales.R` to calculate all-scale environmental covariates from `PRFPD_DTM_leafOff.tiff`.
- Created all-scale environmental outputs under `Quad_Values/Enviro_SHPs/`:
  - `enviro_variables_10m.csv` and `enviro_variables_10m.shp`
  - `enviro_variables_20m.csv` and `enviro_variables_20m.shp`
  - `enviro_variables_50m.csv` and `enviro_variables_50m.shp`
- Environmental outputs include mean DTM elevation (`elev_mean`) and mean Riley topographic roughness index at 5x5, 11x11, and 21x21 windows (`tri5_mean`, `tri11_mean`, `tri21_mean`).
- Confirmed environmental output row counts match source quadrat shapefiles: 2,000 rows for 10 m, 500 rows for 20 m, and 80 rows for 50 m.
- Confirmed no missing values in the required environmental metrics across all three scales.
- Added `tests/test_environmental_variables_all_scales.R`.
- Added `reports/tasks/20260623_environmental_variables_all_scales.md`.
- Added `reports/validation/20260623_environmental_variables_all_scales_validation.md`.
- Added `scripts/3_Analysis/combine_quadrat_analysis_tables.R` to join current spectral heterogeneity, species composition, species diversity, phylogenetic diversity, environmental/topographic values, and quadrat center coordinates into one analysis-ready root-level CSV per scale.
- Created current combined quadrat analysis tables:
  - `quadrat_analysis_10m.csv` with 2,000 rows and 32 columns
  - `quadrat_analysis_20m.csv` with 500 rows and 32 columns
  - `quadrat_analysis_50m.csv` with 80 rows and 32 columns
- Updated the combined tables to exclude per-species composition columns while retaining species diversity summaries, phylogenetic diversity summaries, spectral heterogeneity metrics, environmental/topographic metrics, and center coordinates.
- Added shortened-column documentation in `reports/combined_quadrat_variable_guide.md` and a user-friendly Word guide in `combined_quadrat_variable_guide.docx`.
- Added `scripts/3_Analysis/create_combined_variable_guide_docx.R` to regenerate the Word guide.
- Added `reports/tasks/20260624_combined_quadrat_analysis_tables.md`.
- Added `reports/validation/20260624_combined_quadrat_analysis_tables_validation.md`.
- Validation found no duplicate quadrat IDs, no missing center coordinates, no missing environmental elevation values, no per-species composition columns, and no pixel-count/pair-count/bootstrap/method/geometry metadata columns in the combined CSVs.
- Combined-table spectral missingness reflects upstream raster availability and manual PCA exclusions: missing `spec_sa` values are 97 at 10 m, 15 at 20 m, and 0 at 50 m; missing `spec_pca_mean` values are 256 at 10 m, 64 at 20 m, and 6 at 50 m.
- Revamped `scripts/3_Analysis/multiscale_spectral_biodiversity_analysis.R` so the active first-layer analysis directly answers: what is the relationship between spectral variation and phylogenetic/species diversity?
- The active multiscale analysis now uses two primary spectral variation responses: standardized PCA mean Euclidean distance (`spec_spca_mean`) and spectral angle entropy (`spec_sa`).
- Each spectral variation measure is independently paired with each diversity measure: `phy_faith`, `phy_rao`, `phy_afaith`, `sp_rich`, `sp_shannon`, `sp_simpson`, and `sp_even`.
- Each pairwise relationship is fitted as a simple linear model within scale and reported with Pearson `r`, `R2`, F statistic, F-test p-value, slope, intercept, and Spearman rank correlation.
- The current direct SV-diversity outputs are `reports/analysis/20260710_sv_diversity_pairwise_correlations.md`, `reports/tables/multiscale_spectral_biodiversity/sv_diversity_pairwise_correlations.csv`, `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`, and `reports/tables/multiscale_spectral_biodiversity/sv_diversity_top_pairings.csv`.
- `sv_diversity_analysis_dataset.csv` now carries upstream sampling metadata for both primary SV measures: `sa_n_pixels`, `sa_method`, and `sa_all_pixels_sampled` for SA entropy; and `spca_n_pixels`, `spca_metric_method`, and `spca_euclidean_all_pixels_sampled` for standardized PCA Euclidean distance. `sa_all_pixels_sampled` is TRUE when all retained pixels were sampled for SA entropy and FALSE when the 5,000-pixel cap was used. `spca_euclidean_all_pixels_sampled` is TRUE when standardized PCA Euclidean distance used all retained pixels.
- Current strongest pairings by absolute Pearson `r` are phylogenetic rather than species-diversity measures: 10 m standardized PCA distance with `phy_rao` (r = 0.121, R2 = 0.015), 10 m SA entropy with `phy_rao` (r = 0.106, R2 = 0.011), 20 m standardized PCA distance with `phy_afaith` (r = 0.291, R2 = 0.085), 20 m SA entropy with `phy_rao` (r = 0.180, R2 = 0.033), 50 m standardized PCA distance with `phy_rao` (r = 0.444, R2 = 0.197), and 50 m SA entropy with `phy_rao` (r = 0.340, R2 = 0.115).
- The earlier environment-adjusted AIC model-ranking workflow and PDF reports are superseded as the active first-layer analysis. Environmental and spatial models should be revisited only after the direct pairwise SV-diversity relationships are interpreted.
- Added edge-quadrat, environmental, and SA entropy bootstrap sensitivity outputs in `reports/analysis/20260725_edge_bootstrap_sensitivity_sv_diversity.md`, `reports/figures/edge_bootstrap_sensitivity/`, and `reports/tables/edge_bootstrap_sensitivity/`. These outputs show that edge removal generally weakens all-quadrat correlations slightly at 10 m and 20 m, elevation is negatively associated with spectral variation, environmental screening models retain some incremental phylogenetic-diversity signal, and low-CV SA entropy subsets retain or strengthen phylogenetic correlations while high-CV subsets show weak or absent relationships.
- Documented 10 m and 20 m edge quadrats remain excluded from primary inference; 50 m uses all quadrats because no separate 50 m edge rule is documented.
- Added `reports/tasks/20260624_multiscale_spectral_biodiversity_analysis.md`.
- Added `reports/validation/20260624_multiscale_spectral_biodiversity_analysis_validation.md`.
- Created the manuscript workspace under `Documents/Paper/`. All manuscript-related outputs should be saved there and should include the creation date in the filename using `YYYYMMDD_...`.
- Added `Documents/Paper/20260730_working_outline.md` as the first working manuscript outline. It frames the paper around scale-dependent positive relationships between spectral heterogeneity and phylogenetic diversity, with stronger support for phylogenetic metrics than species-diversity metrics.
- Added `Documents/Paper/20260730_paper_workspace_index.md` to document the paper folder naming convention, current manuscript direction, analysis-ready inputs, source reports, and likely next writing steps.
- Added `Documents/Paper/20260803_writing_style_notes.md` using an R-based workflow to record writing-style guidance from accessible prior Word/context documents. `Documents/Spectral Variation Paper.docx` was locked by another process during this pass, so the note records that limitation; `Documents/SVH Outline.docx` contributed 66 extracted paragraphs.
- Added `Documents/Paper/20260803_markdown_context_digest.md` using R to document the repository Markdown files read before creating the outline.
- Added `Documents/Paper/20260803_detailed_manuscript_outline.docx` as a detailed Word outline focused on correlations between spectral heterogeneity and two focal biodiversity measures: abundance-weighted Faith's PD and Shannon species diversity. Environmental factors are intentionally excluded from this outline stage.
- Added `tools/create_paper_outline_docx.R` as the R workflow used to extract accessible Word text, read Markdown context, write paper notes, and generate the outline DOCX.
- Added `Documents/Paper/20260803_remote_sensing_mdpi_author_requirements.md` for the MDPI Remote Sensing target-journal requirements, including manuscript length, APC, section format, submission notes, references, and fit notes.
- Added `Documents/Paper/20260805_results_pca_mean_distance_afaith.docx` as a Results-section draft focused on standardized PCA mean Euclidean distance versus abundance-weighted Faith's PD across 10 m, 20 m, and 50 m scales. The draft includes figure and table placeholders and avoids environmental-factor interpretation.
- Added `tools/create_pca_afaith_results_docx.R` as the R workflow used to extract the accessible style-reference documents, read latest Markdown/project-result context, and generate the focused Results DOCX.
- Recorded `Documents/Paper/SVH_Results.docx` as a user-created Word results document for manual final edits. Agents may read it for context but should not edit, overwrite, rename, or regenerate it unless explicitly instructed by the user.
- Added `Documents/Paper/20260805_results_spectral_rao_q_afaith.docx` as a separate Results-section draft focused on spectral Rao's Q versus abundance-weighted Faith's PD across 10 m, 20 m, and 50 m scales. The draft uses primary-analysis rows only, includes figure and table placeholders, and avoids environmental-factor interpretation.
- Added `tools/create_rao_q_afaith_results_docx.R` as the R workflow used to read current manuscript/project context and generate the focused spectral Rao's Q Results DOCX.
- Added `Documents/Paper/20260805_discussion_pca_distance_vs_rao_q_scale.docx` as a Discussion-section draft expanding the metric-comparison and scale-dependence interpretation for standardized PCA mean Euclidean distance, spectral Rao's Q, and abundance-weighted Faith's PD.
- Added `tools/create_metric_discussion_docx.R` as the R workflow used to read current manuscript/project context and generate the focused discussion DOCX.
- Added clarified copies `Documents/Paper/20260805_results_spectral_rao_q_afaith_centroid_clarified.docx` and `Documents/Paper/20260805_discussion_pca_distance_vs_rao_q_scale_centroid_clarified.docx` after verifying the actual Rao's Q implementation in `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R`. The clarified wording states that spectral Rao's Q is conceptually equal-weight pairwise squared Euclidean Rao's Q but is computed with the equivalent centroid formula `2 * mean(squared_radius)`.
- Created `Documents/Paper/20260811_methods_draft_spectral_phylogenetic_diversity.docx` as the current methods draft for the manuscript. The draft was generated with R from `tools/create_methods_docx_and_tree_counts.R`, using `Documents/Methods-SampleT.docx`, `Documents/Spectral Variation Paper.docx`, and project markdown context.
- Created tree-count support outputs for methods under `reports/figures/methods_tree_counts/` and `reports/tables/methods_tree_counts/`. The all-plot eligible phylogenetic-diversity tree pool contains 4,007 tree crowns across 50 species. Included-tree summaries are scale-specific and based on focal complete `spec_spca_mean` versus abundance-weighted Faith's PD analysis quadrats.
- Reviewed current Markdown context on 2026-08-17 before beginning the species-diversity versus phylogenetic-diversity correlation analysis.
- Added `reports/execution_plans/20260817_1019_species_phylogenetic_correlation_prep.md`, `reports/tasks/20260817_species_phylogenetic_correlation_prep.md`, and `reports/validation/20260817_species_phylogenetic_correlation_prep_validation.md` to document the analysis setup.
- Added `scripts/3_Analysis/species_phylogenetic_correlation_analysis.R` to compare four species-diversity metrics with three phylogenetic-diversity metrics across all available 10 m, 20 m, and 50 m biodiversity quadrats.
- Created species-phylogenetic concordance outputs under `reports/analysis/20260817_species_phylogenetic_correlation_analysis.md`, `reports/tables/species_phylogenetic_correlation/`, and `reports/figures/species_phylogenetic_correlation/`.
- Species-phylogenetic concordance results now use all available biodiversity quadrats, not edge-filtered primary-analysis quadrats. Faith's PD converges most strongly with species richness at all scales: 10 m r = 0.710, 20 m r = 0.727, and 50 m r = 0.761.
- Phylogenetic Rao's Q aligns most strongly with Simpson diversity at 10 m and 20 m, but only moderately with species richness at 50 m; abundance-weighted Faith's PD weakens strongly with scale and is nearly decoupled from richness at 50 m.
- Incremental spectral-variation models use all complete spectral-diversity cases without edge filtering and show phylogenetic metrics add explanatory signal beyond paired species metrics, with the largest gains from phylogenetic Rao's Q at 50 m for standardized PCA mean distance and SA entropy.
- Added species presence/absence composition clusters and residual-divergence summaries for all 2,580 biodiversity quadrats to identify quadrats with similar species lists but unusually high or low phylogenetic diversity relative to species-diversity expectations.
- Added elevation-gradient all-diversity metric scatterplot figures for 10 m, 20 m, and 50 m, plus `reports/tables/species_phylogenetic_correlation/elevation_adjusted_all_diversity_metric_models.csv`. Each panel reports the incremental R2 and p-value for adding mean elevation after the x-axis diversity metric.
- Named the presence/absence composition types using diagnostic species codes and added species-code heatmaps for 10 m, 20 m, and 50 m under `reports/figures/species_phylogenetic_correlation/`, with companion tables `species_composition_type_names.csv` and `species_composition_type_presence_frequencies.csv`.
- Added `09_all_diversity_metric_pairwise_scatter_composition_cluster_*m.png` figures to show the same all-metric scatterplot panels as the `05_...` figures, with quadrats colored by named species composition type. Each panel legend includes mean silhouette width for each type calculated from that scatterplot's two displayed diversity metrics, with values archived in `species_composition_scatterplot_silhouette.csv`.
- Added `10_all_diversity_metric_pairwise_scatter_individual_ramp_*m.png` figures to show the same all-metric scatterplot panels with quadrats colored by crown-equivalent individuals. The first panel of each scale figure includes the color-ramp legend, and `quadrat_crown_equivalent_individual_totals.csv` stores the quadrat-level values.
- Added `11_all_diversity_metric_pairwise_scatter_species_count_ramp_*m.png` figures to show the same all-metric scatterplot panels with quadrats colored by the number of species present. The first panel of each scale figure includes the color-ramp legend, and the plotted values are stored as `present_species_count`.
- Added standardized moderation models testing whether species count or crown-equivalent individuals change all-metric relationship slopes. The full model table is `species_individual_moderated_metric_relationships.csv`, with BH-significant interaction terms in `significant_species_individual_moderated_metric_relationships.csv`.
- Added `scripts/3_Analysis/spectral_heterogeneity_relationship_analysis.R` to compare all seven focal spectral heterogeneity measures across scales: spectral angle entropy, raw PCA alpha-hull area, raw PCA spectral Rao's Q, raw PCA mean Euclidean distance, standardized PCA alpha-hull area, standardized PCA spectral Rao's Q, and standardized PCA mean Euclidean distance.
- Created spectral heterogeneity relationship outputs under `reports/analysis/20260817_spectral_heterogeneity_relationship_analysis.md`, `reports/tables/spectral_heterogeneity_relationships/`, and `reports/figures/spectral_heterogeneity_relationships/`. The figure sets include integrated all-seven-metric scatterplots plus elevation, species-composition, crown-equivalent-individual, species-count, mean 563 nm pixel-brightness, region-specific brightness highlighted versions, and direct spectral-metric versus regional-illumination scatterplots for each quadrat scale. Pixel brightness is summarized per quadrat from retained sunlit pixels in `quadrat_pixel_brightness_summary.csv`, including 563 nm, blue 450-495 nm, green 500-570 nm, red 620-750 nm, and near-infrared 750-998 nm. Direct illumination relationships are archived in `spectral_metric_regional_illumination_correlations.csv`.
- Added `scripts/3_Analysis/spectral_biodiversity_all_metric_scatter_analysis.R` to create one 7-row by 3-column scatterplot figure for each focal spectral heterogeneity measure, with biodiversity metrics as rows and 10 m, 20 m, and 50 m scales as columns. Each panel reports complete-case `n`, Pearson `r`, simple-regression `R2`, and p-value.
- Created all-metric spectral-biodiversity scatterplot outputs under `reports/analysis/20260819_spectral_biodiversity_all_metric_scatter_analysis.md`, `reports/tables/spectral_biodiversity_all_metrics/`, and `reports/figures/spectral_biodiversity_all_metrics/`. The strongest absolute relationships are concentrated at 50 m and involve phylogenetic metrics, led by standardized PCA alpha hull with abundance-weighted Faith's PD (r = 0.497, R2 = 0.247), standardized PCA alpha hull with phylogenetic Rao's Q (r = 0.487, R2 = 0.238), and standardized PCA mean distance with phylogenetic Rao's Q (r = 0.444, R2 = 0.197).
- Added `scripts/3_Analysis/bootstrap_metric_stat_relationship_analysis.R` to create scale-specific scatterplot figures comparing bootstrapped metric values against bootstrap SD, CV, and SE.
- Created bootstrap metric-stat relationship outputs under `reports/analysis/20260817_bootstrap_metric_stat_relationship_analysis.md`, `reports/tables/bootstrap_metric_relationships/`, and `reports/figures/bootstrap_metric_relationships/`. The final-pipeline figure set plots spectral angle entropy against 70-iteration bootstrap SD/CV/SE. A separate all-seven-metric figure set uses the fixed-4,000-pixel condition from the sample-size sensitivity branch and is labeled as stability diagnostics rather than final PCA metric uncertainty.
- Added `scripts/3_Analysis/sa_bootstrap_context_relationship_analysis.R` to plot the final SA entropy bootstrap SD, CV, and SE against elevation, species count, crown-equivalent individual count, species composition type, overall brightness, and blue, green, red, and near-infrared brightness.
- Created SA bootstrap context outputs under `reports/analysis/20260817_sa_bootstrap_context_relationship_analysis.md`, `reports/tables/bootstrap_metric_relationships/`, and `reports/figures/bootstrap_metric_relationships/`. The strongest continuous relationships are 10 m brightness relationships, especially NIR brightness with SA bootstrap CV (r = 0.269, R2 = 0.072) and NIR brightness with SA bootstrap SD/SE (r = 0.262, R2 = 0.069). Species composition effects are small: the largest eta-squared is about 0.027 at 50 m and is not significant, while 10 m composition CV is statistically detectable but tiny (eta-squared = 0.006).
- Created `Documents/Paper/20260819_final_research_direction.md` to focus the final manuscript direction around standardized PCA alpha hull, standardized PCA mean distance, and spectral angle entropy as primary spectral metrics; phylogenetic Rao's Q and abundance-weighted Faith's PD as primary phylogenetic metrics; and Shannon diversity as the main species-diversity contrast. Added supporting execution, task, and validation notes under `reports/`.
- Created `Documents/Paper/20260819_final_research_narrative.md` and root-level `Checklist.md` after reading `20260815_evan_next_steps.docx` and the read-only current draft `Documents/Paper/SVH_Paper.docx`. The active manuscript direction now uses the working title "Assessing the Spectral Variation Hypothesis Across Spatial Scales Using PCA-based Canopy Spectral Heterogeneity" and centers scale dependence, PCA metric behavior, biodiversity metric behavior, transformations, and correlation mechanisms.
- Added `scripts/3_Analysis/final_research_direction_analysis.R` to complete the current `Checklist.md` analysis-needed items that can be addressed from existing 10 m, 20 m, and 50 m outputs. The script writes `reports/analysis/20260819_final_research_direction_analysis.md`, nine CSV tables under `reports/tables/final_research_direction/`, and four figures under `reports/figures/final_research_direction/`.
- The final research direction analysis produced 273 all spectral-biodiversity relationship rows, 234 spectral metric pairwise rows, 63 biodiversity metric pairwise rows, 20,412 transformation rows, 252 best-transform summary rows, 468 metric-driver rows, 42 PCA computation-method rows, 45 Moran's I diagnostic rows, and 120 variogram rows.
- Current final research direction results confirm that the strongest PCA spectral-biodiversity relationships remain concentrated at 50 m and involve phylogenetic metrics. The top relationships are standardized PCA alpha hull with abundance-weighted Faith's PD (r = 0.497, R2 = 0.247), standardized PCA alpha hull with phylogenetic Rao's Q (r = 0.487, R2 = 0.238), and standardized PCA mean distance with phylogenetic Rao's Q (r = 0.444, R2 = 0.197).
- Transformation diagnostics show the largest Pearson-r gains for standardized PCA Rao's Q relationships, especially at 10 m and 20 m, but these transformations should be interpreted as diagnostics rather than automatic replacements for untransformed results.
- Spatial diagnostics found significant Moran's I in priority variables and many OLS residuals, so ordinary least-squares p-values should remain cautious until a final spatial modeling or spatial sensitivity layer is selected.
- PCA computation-method checks show final PCA mean-distance and Rao-style metrics are marked `all_pixels` for included quadrats, while alpha-hull and 3D hull metrics use all pixels for many 10 m quadrats and sampled-pixel/fallback methods for larger 20 m and 50 m quadrats. Fixed 4,000-pixel outputs remain sample-size sensitivity results rather than the final main PCA metric values.
- Updated `reports/directory_map.md` to reflect the `Documents/Paper/` manuscript-output convention and the current edge/bootstrap sensitivity output folders.
- Created required governance directories under `reports/` and `logs/`.
- Created an execution plan for repository documentation and cleanup.
- Created baseline `reports/directory_map.md`.
- Created baseline `reports/data_dictionary.md`.
- Removed obsolete `.Rhistory` session files from the repository.
- Updated `CODEX_AGENT_GUIDELINES.md` with R/Rtools access guidance, subsampled unit testing expectations, and current data assumptions.
- Removed authorized legacy directories named `old`, `Outdated`, and `Currently not relevant`.
- Created `tests/` for future `testthat` unit tests.
- Updated R package library use to the standard R 4.2 paths:
  - `C:/Users/PaintRock/AppData/Local/R/win-library/4.2`
  - `C:/Program Files/R/R-4.2.3/library`
- Discarded raster checking/reproduction artifacts at user request; those temporary scripts, reports, and outputs should not guide future analysis.

## Pending Work

- Review current R scripts for package dependencies, inputs, outputs, assumptions, and reproducibility risks.
- Decide how to handle missing 10 m and 20 m quadrats that do not have current raster summaries.
- Carry bootstrap quality-control fields into downstream tables, especially `boot_sd`, `boot_cv`, and `method`.
- Run downstream model sensitivity checks with high-CV quadrats flagged or excluded.
- When modifying R workflows, first review nearby scripts in `scripts/` to match the project's existing loop structure and file-processing pattern.
- Standardize script and directory names that contain spaces or spelling errors after checking dependencies.
- Add formal `testthat` tests once core functions are modularized.
- Expand shapefile and raster data dictionary entries using R geospatial packages.
- After interpreting the direct pairwise SV-diversity layer, add spatial and environmental sensitivity checks as a second analysis layer.
- Join bootstrap quality-control fields (`boot_sd`, `boot_cv`, and `method`) into a companion modeling table and rerun sensitivity checks that flag or exclude high-CV quadrats.
- Decide which species-phylogenetic concordance results belong in the main manuscript versus supplement.
- Use the species-composition residual-divergence table to inspect individual quadrats or clusters where phylogenetic diversity diverges from species-diversity expectations.
- Use `Documents/Paper/20260819_final_research_narrative.md`, `Checklist.md`, and `reports/analysis/20260819_final_research_direction_analysis.md` to update the manuscript outline, Results drafts, figure inventory, and table inventory.
- Verify remaining methods details flagged from `SVH_Paper.docx`, especially tree filtering, species counts, CRS language, crown-radius calculation, and abundance-weighted Faith's PD.
- Interpret whether metric magnitude is confounded with bootstrap uncertainty, especially for the fixed-4,000 sample-size branch where larger PCA-space dispersion can mechanically increase SD and SE.

## Known Issues

- R/Rscript was not available on the command path during this documentation pass.
- Several script and directory names contain spaces or misspellings, which can complicate reproducible command-line execution.
- Some current CSV files contain unnamed first columns.
- `Quad_Values/20m_spectral_sp.csv` has 64 missing values for each primary spectral metric.
- `Rscript` is available at `C:/Program Files/R/R-4.2.3/bin/Rscript.exe`, matching the required R 4.2.3 compatibility target, though it is not available as `Rscript` on the active shell path.
- Codex can access the R 4.2 user package library under `C:/Users/PaintRock/AppData/Local/R/win-library/4.2` when commands are run with permission outside the workspace sandbox; sandboxed Rscript may still report only the base R library.

## Technical Debt

- Some scripts appear exploratory or duplicated, including `modeling_test.R` and `scratch paper.R`.
- Current workflows likely use hard-coded paths and should be audited.
- The temporary project-local `r_libs/` package library was removed from the working directory.

## Important Assumptions

- Legacy directories named `old`, `Outdated`, and `Currently not relevant` have been removed from the active repository.
- As of 2026-07-06, current derived spectral, diversity, environmental, and topographic outputs are visible under `Quad_Values/`.
- As of 2026-07-30, manuscript-related outputs belong under `Documents/Paper/` and should include their creation date in the filename.
- `Documents/Paper/SVH_Results.docx` is a user-created manual-editing document. Treat it as read-only context for agent work unless the user explicitly asks the agent to edit that specific file.
- As of 2026-08-11, `Documents/Paper/20260811_methods_draft_spectral_phylogenetic_diversity.docx` is the current generated methods draft. Its supporting tree-count charts are in `reports/figures/methods_tree_counts/`, and its supporting species-count CSV tables are in `reports/tables/methods_tree_counts/`.
- The methods tree-count tables use the headers `Genus`, `Species`, and `Number of individuals` for the individual all-plot and included-scale CSVs. The combined summary table adds scale-specific included-count columns for comparison.
- The primary integrated 20 m analysis table is `Quad_Values/20m_spectral_sp.csv`.
- The current all-scale plant-diversity outputs are `Quad_Values/Diversity_SHPs/plant_diversity_10m.csv`, `plant_diversity_20m.csv`, and `plant_diversity_50m.csv`, with matching shapefiles.
- The current all-scale environmental outputs are `Quad_Values/Enviro_SHPs/enviro_variables_10m.csv`, `enviro_variables_20m.csv`, and `enviro_variables_50m.csv`, with matching shapefiles.
- The current combined analysis-ready tables are `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`.
- The combined analysis tables intentionally exclude per-species composition columns, pixel-count, pair-count, bootstrap replicate, processing method, manual exclusion, and geometry metadata columns.
- Quadrat center coordinate columns in the combined analysis tables are `center_x` and `center_y`, derived from plant-diversity shapefile polygon centroids in NAD83 / UTM zone 16N.
- Raw spectral inputs are in `HSI_NA_trimmed/`.
- The user has independently tested and confirmed the partitioned quadrat spectra as the inputs to use moving forward.
- Confirmed partitioned quadrat spectra are in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m`.
- Current smoothed spectra derived from those confirmed inputs are in `Quad_Spectra/10m_smooth`, `Quad_Spectra/20m_smooth`, and `Quad_Spectra/50m_smooth`.
- Current smoothed 5 nm spectra are in `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, and `Quad_Spectra/50m_smooth_5nm`.
- The current spectral heterogeneity workflow should consume the current `_smooth_5nm` spectra, not older `_resampled_5nm` or `_smoothed_5nm` folders.
- Spectral angle entropy values should be produced per quadrat at each resolution; 10 m uses `sub_id` IDs such as `0_a`, 20 m uses numeric `Name` IDs, and 50 m uses `Name` IDs such as `sub50_1`.
- Final-pipeline per-iteration bootstrap outputs currently visible for metric-stat uncertainty are spectral angle entropy outputs only. All-seven-metric bootstrap SD/CV/SE figures use the fixed-4,000-pixel sample-size sensitivity branch and should not be interpreted as final-value uncertainty for PCA-based metrics unless main-pipeline per-iteration PCA metric outputs are later located or regenerated.
- Plant-diversity outputs now use the same `quad_id` convention as spectral diversity outputs, so downstream joins should use `scale` plus `quad_id`.
- Environmental outputs now use the same `quad_id` convention as plant-diversity and spectral diversity outputs, so downstream joins should use `scale` plus `quad_id`.
- Environmental outputs were derived from `PRFPD_DTM_leafOff.tiff`; quadrat polygons were transformed to the DTM CRS before raster extraction.
- Riley topographic roughness index outputs were calculated from 5x5, 11x11, and 21x21 DTM moving windows and then averaged within each quadrat.
- Plant-diversity species matrix values represent summed crown-overlap proportions by species and quadrat, preserving the existing crown-buffer/intersection logic.
- Exact all-pixel spectral angle entropy is only reasonable for small masked quadrats because pairwise angle counts scale quadratically with pixel count; the workflow should use bootstrap `boot_mean` as the primary fallback value for larger quadrats.
- Current spectral heterogeneity outputs used 70 bootstrap replicates, 5,000 sampled pixels per bootstrap, and 10,000 sampled pixel-pairs per large bootstrap replicate.
- Validation found 10 m summary rows for 1,909 rasters, with 1,866 bootstrap-mean values, 37 exact all-pixel values, and 6 insufficient-pixel outputs.
- Validation found 20 m summary rows for 485 rasters, with 477 bootstrap-mean values and 8 exact all-pixel values.
- Validation found 50 m summary rows for 80 rasters, with 79 bootstrap-mean values and 1 exact all-pixel value.
- The 10 m output shapefile has 97 missing `spec_ent` values: 91 shapefile quadrats without matching current raster outputs and 6 insufficient-pixel raster summaries.
- The 20 m output shapefile has 15 missing `spec_ent` values from shapefile quadrats without matching current raster outputs.
- The 50 m output shapefile has no missing `spec_ent` values.
- The current PCA basis samples the entire available 10 m raster footprint after shadow masking and does not sample from 20 m or 50 m rasters. Downstream PCA metric values still apply the manual atmospheric/cloud exclusion policy for metric reporting. Current raster-output exclusions are 165 of 1,909 at 10 m, 49 of 485 at 20 m, and 6 of 80 at 50 m.
- The current 50 m manual atmospheric/cloud exclusions are `sub50_80`, `sub50_79`, `sub50_71`, `sub50_70`, `sub50_62`, and `sub50_53`.
- The combined 10 m spectral heterogeneity CSV has 1,909 rows and 47 columns. It has 6 missing SA entropy values from insufficient-pixel SA entropy summaries, and 165 missing raw and standardized PCA Euclidean, Rao Q, and 3D PCA hull values from manual exclusions. Raw and standardized alpha-hull values have 168 missing values: 165 manual exclusions plus 3 alpha-hull failures after deterministic sampled fallback.
- The combined 20 m spectral heterogeneity CSV has 485 rows and 47 columns. It has no missing SA entropy values and 49 missing raw and standardized PCA-dependent values from manual exclusions. The combined 20 m shapefile has 15 missing SA entropy values from shapefile quadrats without matching current raster outputs, and 64 missing PCA-dependent values from those unmatched polygons plus 49 manual exclusions.
- The combined 50 m spectral heterogeneity CSV and shapefile have 80 rows/features. They have no missing SA entropy values and 6 missing raw and standardized PCA-dependent values from manual exclusions.
- The supplemental three-axis hull output is a convex hull volume/area in global PC1-PC3 space, not a true 3D alpha hull, because `alphashape3d` is not currently installed.
- Bootstrap QC found low median CV at all scales, about 0.4% to 0.5%, but a non-trivial high-CV tail.
- Using a 5% bootstrap-CV flag, 372 of 1,866 10 m bootstrap-mean quadrats, 60 of 477 20 m quadrats, and 11 of 79 50 m quadrats would be flagged.
- Using a 10% bootstrap-CV flag, 127 of 1,866 10 m bootstrap-mean quadrats and 12 of 477 20 m quadrats would be flagged; no 50 m quadrats exceed 10%.
- Spectral heterogeneity distributions show mild right skew at 10 m, near-symmetry at 20 m, and stronger right skew at 50 m.
- `Quad_Spectra/10m_test`, `Quad_Spectra/20m_test`, and `Quad_Spectra/50m_test` are testing/validation artifacts unless the user explicitly promotes them.
- Quadrat boundary files are in `Quad_Scale_SHPs/`.
- Spectral products that produce spectral heterogeneity values are in `Quad_Values/Spectral_diversitySHPs/`.
- Unnamed first CSV columns are likely index columns.

## Next Recommended Actions

1. Use `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv` as the current analysis-ready tables for downstream biodiversity/spectral/environment modeling.
2. Continue manuscript preparation in `Documents/Paper/`, using dated filenames for all outline, draft, figure-inventory, table-inventory, and narrative files.
3. Decide whether 10 m and 20 m missing spectral quadrats should remain `NA`, be excluded, or be tracked in a missing-quadrat report.
4. Decide whether bootstrap quality-control fields such as `boot_sd`, `boot_cv`, and `method` should be added to a separate QC companion table rather than the value-only combined tables.
5. Use the July 25 edge/bootstrap sensitivity outputs to decide which robustness results belong in the main manuscript versus supplement.
6. Use `Documents/Paper/20260819_final_research_narrative.md` as the current paper-facing narrative document for final manuscript outputs.
7. Use root-level `Checklist.md` as the current gap list for missing analyses, writing, methods/results updates, and draft inaccuracies.
8. Decide whether the spatial diagnostics are sufficient as a sensitivity layer or whether to extend them into spatial GLS/GAM or spatial eigenvector models before final inferential language.
9. Use `C:/Program Files/R/R-4.2.3/bin/Rscript.exe` for R execution unless the shell path is updated later.
