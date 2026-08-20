# Directory Map

Last updated: 2026-08-17

Working directory:
`C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity`

Current live derived-output directory:
`Quad_Values/`

`Quad_Values/` is the canonical directory for derived biodiversity, spectral diversity, environmental, topographic, and index outputs.

Current manuscript-output directory:
`Documents/Paper/`

All manuscript-related outputs should be saved under `Documents/Paper/` with the creation date in the filename using `YYYYMMDD_...`.

## Top-Level Inventory

| Directory | Purpose | File Count | Approx. Size |
|---|---|---:|---:|
| `.agents/` | Agent-side metadata and coordination files | 0 | <0.001 GB |
| `Documents/` | Manuscript drafts, paper workspace files, reports, presentations, workshop materials, and generated PDF reports | 34 | 0.022 GB |
| `Functions/` | Legacy/support R functions, package archive, notebooks, and helper modules | 20 | 0.033 GB |
| `HSI_NA_trimmed/` | Current trimmed raw hyperspectral imagery in ENVI-style raster format | 75 | 57.022 GB |
| `logs/` | Project logs required by governance standards | 10 | <0.001 GB |
| `Quad_Scale_SHPs/` | Quadrat boundary shapefiles and KMLs for 10 m, 20 m, and 50 m analysis scales | 14 | 0.002 GB |
| `Quad_Spectra/` | Quadrat-level spectral raster extracts and derived products, including current `_smooth` and `_smooth_5nm` outputs | 37,547 | 391.497 GB |
| `Quad_Values/` | Current derived biodiversity, spectral diversity, environmental, topographic, and index outputs | 136 | 0.722 GB |
| `reports/` | Governance reports, project inventories, task reports, plans, validation notes, analysis reports, figures, and tables | 481 | 0.178 GB |
| `scripts/` | R workflow scripts for preprocessing, index creation, analysis, utilities, and visuals | 67 | 0.001 GB |
| `Shadow_vs_sunlit_SHP/` | Shadow/sunlit training or classification shapefiles | 12 | <0.001 GB |
| `Test/` | Ad hoc test/projection script area | 1 | <0.001 GB |
| `tests/` | Formal unit tests, smoke checks, and lightweight reproduced-raster validation artifacts | 4 | 0.565 GB |
| `VegIndex_NA_trimmed/` | Vegetation-index raster outputs aligned with trimmed imagery | 75 | 16.485 GB |
| `Visuals/` | Current exported figures and visual raster products | 5 | 0.001 GB |

## Repository Structure

```text
Spectral_Diversity/
    CODEX_AGENT_GUIDELINES.md
    RESEARCH_OBJECTIVES.md
    .gitignore
    .Rhistory
    markdown_script.R
    global_pca.rds
    PR_tree_DL.csv
    51sp_taxanomy.csv
    species_summary.csv
    spec_shannaon_DIV_with Geometry.csv
    51sp_codeTree.pdf
    51sp_codeTree.png
    PlantTrait_Phylogene_Chart.drawio
    hsi_illumination_class.tif.aux.xml
    PRFPD_DTM_leafOff.tiff
    quadrat_analysis_10m.csv
    quadrat_analysis_20m.csv
    quadrat_analysis_50m.csv
    combined_quadrat_variable_guide.docx

    Documents/
        manuscript drafts, PDF model reports, ARD presentation, workshop document
        Paper/
            20260730_paper_workspace_index.md
            20260730_working_outline.md
            20260803_detailed_manuscript_outline.docx
            20260803_markdown_context_digest.md
            20260803_remote_sensing_mdpi_author_requirements.md
            20260803_writing_style_notes.md
            20260805_discussion_pca_distance_vs_rao_q_scale.docx
            20260805_discussion_pca_distance_vs_rao_q_scale_centroid_clarified.docx
            20260805_results_pca_mean_distance_afaith.docx
            20260805_results_spectral_rao_q_afaith.docx
            20260805_results_spectral_rao_q_afaith_centroid_clarified.docx
            SVH_Results.docx
        PDFs/
            spectral_biodiversity_multiscale_findings.pdf
            spectral_biodiversity_model_appendix.pdf

    Functions/
        1_LCE_derivs.R
        2_LCE_veg_index.R
        3_LCE_main_function.R
        dataframe_operations.R
        lecospectR.R
        lecospectR_1.0.tar.gz
        model_support.R
        pfts.R
        pipeline.R
        probably_not_used.R
        raster_operations.R
        site_specific_processing.R
        spectral_operations.R
        training_utilities.R
        type_conversion.R
        utilities.R
        validation.R
        visualization.R
        Segmentation/
        test-visalization.ipynb

    HSI_NA_trimmed/
        trimmed hyperspectral raster cubes and paired metadata

    Quad_Values/
        20m_spectral_sp.csv
        20m_SA_smooth_masked.csv
        20m_SA_smooth_masked_7_11.csv
        20m_SA_entrop_boot_results.csv
        20m_SA_entrop_boot100_results.csv
        10m_SA_entropy_smooth_masked_5nm_summary.csv
        20m_SA_entropy_smooth_masked_5nm_summary.csv
        50m_SA_entropy_smooth_masked_5nm_summary.csv
        10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv
        20m_spectral_heterogeneity_smooth_masked_5nm_summary.csv
        50m_spectral_heterogeneity_smooth_masked_5nm_summary.csv
        *_PCA_metrics_smooth_masked_5nm_summary.csv
        *_SA_entropy_smooth_masked_5nm_boot_long.csv
        *_SA_entropy_smooth_masked_5nm_boot_wide.csv
        Diversity_Rasters/
        Diversity_SHPs/
            plant_diversity_10m.*
            plant_diversity_20m.*
            plant_diversity_50m.*
        Enviro_SHPs/
            enviro_variables_10m.*
            enviro_variables_20m.*
            enviro_variables_50m.*
        Other_variables/
            Slope_Aspect.csv
            Slope_Aspect_rasters/
        Spectral_diversitySHPs/
            global_pca_smooth_masked_5nm.rds
            global_pca_smooth_masked_5nm_variance_explained.csv
            standardized_PCA_global_pca_smooth_masked_5nm.rds
            standardized_PCA_global_pca_smooth_masked_5nm_variance_explained.csv
            spectral_heterogeneity_10m_smooth_masked_5nm.*
            spectral_heterogeneity_20m_smooth_masked_5nm.*
            spectral_heterogeneity_50m_smooth_masked_5nm.*
        Spectral_diversitySHPs_first/

    Quad_Scale_SHPs/
        PR_10m.*
        PR_20m.*
        PR_50m.*
        KMLs/

    Quad_Spectra/
        10m/
        10m_smooth/
        10m_smooth_5nm/
        10m_smoothed_5nm_first/
        10m_test/
        10m_VegIndex/
        20m/
        20m_RGB/
        20m_smooth/
        20m_smooth_5nm/
        20m_smoothed_5nm_first/
        20m_test/
        20m_VegIndex/
        50m/
        50m_RGB/
        50m_RGB_smooth/
        50m_smooth/
        50m_smooth_5nm/
        50m_smoothed_5nm_first/
        50m_test/
        50m_VegIndex/

    scripts/
        1_Data Processing/
        2_Indices Creation/
            Enviro_Variables/
            Plant_diversity/
            Spectral_diversity/
        3_Analysis/
        auxilary/
        visuals/

    Shadow_vs_sunlit_SHP/
        Shadow_vs_sunlit.*
        Shadow_vs_sunlit_Rev.*

    Test/
        UTM projected.R

    tests/
        test_sa_entropy_bootstrapping.R
        smoke_sa_entropy_bootstrapping.R
        test_environmental_variables_all_scales.R
        reproduced_quad_spectra/

    VegIndex_NA_trimmed/
        vegetation-index raster outputs

    Visuals/
        exported figures and map products

    reports/
        analysis/
        architecture/
        combined_quadrat_variable_guide.md
        execution_plans/
        figures/
            bootstrap_metric_relationships/
            bootstrap_variation/
            edge_bootstrap_sensitivity/
            multiscale_spectral_biodiversity/
            pca_loading_spectral_regions/
            pc1_mean_reflectance_correlation/
            sample_size_effects/
                alpha_hull_area/
                pca_mean_distance/
                sa_entropy/
                spectral_rao_q/
            sa_entropy_sample_size_effects/
            spectral_biodiversity_all_metrics/
            species_phylogenetic_correlation/
            spectral_heterogeneity/
        history/
        tables/
            bootstrap_metric_relationships/
            bootstrap_variation/
            edge_bootstrap_sensitivity/
            multiscale_spectral_biodiversity/
            pca_loading_spectral_regions/
            pc1_mean_reflectance_correlation/
            sample_size_effects/
                alpha_hull_area/
                pca_mean_distance/
                sa_entropy/
                spectral_rao_q/
            sa_entropy_sample_size_effects/
            spectral_biodiversity_all_metrics/
            species_phylogenetic_correlation/
        tasks/
        validation/

    logs/
        resampling_*.log
        smoothing_*.log
        sa_entropy_bootstrapping_progress.log
        spectral_heterogeneity_all_metrics_*.log
```

## Current Script Inventory

### `scripts/1_Data Processing/`

- `aspect_slope.R`
- `HSI_NA_trim.R`
- `HSI_quad_crop_refined.R`
- `resampling.R`
- `smoothing.R`
- `Sunlit_vs_Shadow_AUC.R`
- `testing_reproducibilty.R`

### `scripts/2_Indices Creation/`

- `adding_other_variables.R`
- `adding_to_CSV.R`
- `Indices_CSV.R`
- `topo_Index.R`
- `Enviro_Variables/elevation.R`
- `Enviro_Variables/environmental_variables_all_scales.R`
- `Plant_diversity/Missing_quads.R`
- `Plant_diversity/plant_diversity_all_scales.R`
- `Plant_diversity/phylogenetic_diversity.R`
- `Plant_diversity/sp_weighted_matrix.R`
- `Plant_diversity/species_diversity.R`
- `Spectral_diversity/HSI_global_PCA.R`
- `Spectral_diversity/SA_entropy_bootstrapping.R`
- `Spectral_diversity/run_sa_entropy_scale.R`
- `Spectral_diversity/spectral_heterogeneity_all_metrics.R`
- `Spectral_diversity/spectral_angle_entropy.R`
- `Spectral_diversity/spectral_angle_entropy_movingwindow.R`
- `Spectral_diversity/spectral_angle_entropy_raster.R`

### `scripts/3_Analysis/`

- `Analysis_1.R`
- `Analysis_PDF.R`
- `bootstrap_metric_stat_relationship_analysis.R`
- `bootstrap_variation_analysis.R`
- `Bootstrap_Analysis.R`
- `combine_quadrat_analysis_tables.R`
- `create_combined_variable_guide_docx.R`
- `edge_and_bootstrap_sensitivity_sv_diversity.R`
- `final_research_direction_analysis.R`
- `global_vs_local.R`
- `linearity_test.R`
- `LLM.R`
- `multiscale_spectral_biodiversity_analysis.R`
- `modeling.R`
- `modeling_test.R`
- `ouliers.R`
- `random Forest.R`
- `scratch paper.R`
- `sa_entropy_sample_size_effects.R`
- `sa_bootstrap_context_relationship_analysis.R`
- `species_phylogenetic_correlation_analysis.R`
- `Spatial_random_Forest.R`
- `spectral_biodiversity_all_metric_scatter_analysis.R`
- `spectral_heterogeneity_relationship_analysis.R`
- `spectral_metric_sample_size_effects.R`
- `within_between_variation.R`

### `scripts/auxilary/`

- `kml_shp_conversion.R`
- `missing tiles.R`
- `RGB_Conversion.R`
- `sp_summary.R`
- `Stats.R`
- `vegindex summerization effects.R`

### `scripts/visuals/`

- `edge_quads.R`
- `interactrive_crown_radius.R`
- `shadow_mask.R`
- `smooth_spectra_compare.R`
- `Spatial Correlation.R`
- `spatial_corr_plot.R`
- `Spectral_Signatures.R`

### Remaining Review Candidate Scripts

- `scripts/3_Analysis/scratch paper.R`
- `scripts/3_Analysis/modeling_test.R`

## Current Function Inventory

### `Functions/`

- `1_LCE_derivs.R`
- `2_LCE_veg_index.R`
- `3_LCE_main_function.R`
- `dataframe_operations.R`
- `lecospectR.R`
- `lecospectR_1.0.tar.gz`
- `model_support.R`
- `pfts.R`
- `pipeline.R`
- `probably_not_used.R`
- `raster_operations.R`
- `site_specific_processing.R`
- `spectral_operations.R`
- `training_utilities.R`
- `type_conversion.R`
- `utilities.R`
- `validation.R`
- `visualization.R`
- `Segmentation/segmentation_test.ipynb/segmentation_test.ipynb`
- `test-visalization.ipynb`

## File-Type Summary

| Extension | Count | Notes |
|---|---:|---|
| no extension | 16,656 | ENVI raster data files and extensionless raster products |
| `.xml` | 12,438 | Raster auxiliary metadata, mainly spectral products |
| `.hdr` | 12,348 | ENVI raster headers |
| `.tif` | 579 | Raster imagery and derived products |
| `.png` | 363 | Exported figures, including bootstrap metric-stat, bootstrap context, multiscale, all-metric spectral-biodiversity, species-phylogenetic, spectral-relationship, composition-type, and sample-size diagnostics |
| `.R` | 93 | R workflow scripts, including root-level, function, test, tool, and ad hoc scripts |
| `.csv` | 139 | Tabular tree, taxonomy, spectral, biodiversity, environmental, bootstrap, validation summary, and combined analysis data |
| `.md` | 89 | Project context, governance, reports, tasks, plans, validation notes, variable guides, and dated paper workspace notes |
| `.shp/.dbf/.prj/.shx` | 26 each | Shapefile component sets |
| `.envi` | 25 | Vegetation-index raster products |
| `.pdf` | 15 | Reports and phylogeny output |
| `.sample` | 14 | Sample/sidecar files associated with raster products |
| `.docx` | 16 | Manuscript/workshop documents, paper outline/results drafts, and variable guides |
| `.log` | 10 | Processing stdout, stderr, and progress logs under `logs/` |
| `.cpg` | 8 | Shapefile code-page sidecar files |
| `.Rhistory` | 6 | R session-history files are present and should be treated as non-analytical artifacts |
| `.rds` | 2 | Serialized R objects, including PCA model outputs |
| `.ipynb` | 2 | Function-development or segmentation notebooks |

## Cleanup Notes

- Removed earlier `.Rhistory` session files on 2026-06-15; as of 2026-06-23, `.Rhistory` files are again present in the live tree and should not be treated as analytical inputs.
- Removed authorized `old/`, `Outdated/`, and `Currently not relevant/` directories on 2026-06-15.
- The cleanup removed approximately 107.8 GB from `Quad_Spectra/old/`.
- Current live derived outputs are in `Quad_Values/`.
- Current spectral heterogeneity products are in `Quad_Values/Spectral_diversitySHPs/`.
- `Quad_Values/Spectral_diversitySHPs_first/` and `Quad_Spectra/*_smoothed_5nm_first/` are retained first-pass/earlier output folders unless explicitly promoted.
- The user has independently tested and confirmed the partitioned quadrat spectra to use moving forward.
- Confirmed partitioned quadrat spectra are in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m`.
- Current smoothed spectra generated from those confirmed inputs are in `Quad_Spectra/10m_smooth`, `Quad_Spectra/20m_smooth`, and `Quad_Spectra/50m_smooth`.
- Current smoothed 5 nm spectra generated from the smoothed outputs are in `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, and `Quad_Spectra/50m_smooth_5nm`.
- The current spectral angle entropy workflow should pick up from the `_smooth_5nm` folders and write outputs under `Quad_Values/`.
- The current PCA-dependent spectral heterogeneity workflow is `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R`. The PCA basis samples only current 10 m smoothed 5 nm rasters after shadow masking, requesting up to 450 illuminated pixels per raster. Downstream raw and standardized PCA metric outputs still preserve documented atmospheric/cloud exclusions as missing metric values.
- Current per-scale spectral heterogeneity outputs have been generated for 10 m, 20 m, and 50 m using the current `_smooth_5nm` spectra.
- Current per-scale plant diversity outputs are visible for 10 m, 20 m, and 50 m under `Quad_Values/Diversity_SHPs/`, with `quad_id` values aligned to the current spectral heterogeneity outputs.
- Current per-scale environmental outputs are visible for 10 m, 20 m, and 50 m under `Quad_Values/Enviro_SHPs/`, with `quad_id` values aligned to current plant-diversity and spectral outputs.
- Current root-level combined analysis tables are `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`; their shortened-column guides are `reports/combined_quadrat_variable_guide.md` and `combined_quadrat_variable_guide.docx`.
- Current multiscale spectral-biodiversity analysis is the direct SV-diversity pairwise correlation workflow in `scripts/3_Analysis/multiscale_spectral_biodiversity_analysis.R`. The active report is `reports/analysis/20260710_sv_diversity_pairwise_correlations.md`, with companion figures and tables under `reports/figures/multiscale_spectral_biodiversity/` and `reports/tables/multiscale_spectral_biodiversity/`. Older PDF reports in `Documents/PDFs/` are superseded historical context.
- Current all-metric spectral-biodiversity scatterplot analysis is in `scripts/3_Analysis/spectral_biodiversity_all_metric_scatter_analysis.R`. The active report is `reports/analysis/20260819_spectral_biodiversity_all_metric_scatter_analysis.md`, with companion tables under `reports/tables/spectral_biodiversity_all_metrics/` and 13 figures under `reports/figures/spectral_biodiversity_all_metrics/`. Each figure holds one spectral heterogeneity measure, including spectral angle entropy plus raw and standardized PCA alpha, Rao's Q, mean distance, convex hull, 3D hull volume, and 3D hull area metrics, with the seven biodiversity metrics as rows and 10 m, 20 m, and 50 m as columns.
- Current final research direction analysis is in `scripts/3_Analysis/final_research_direction_analysis.R`. The active report is `reports/analysis/20260819_final_research_direction_analysis.md`, with companion tables under `reports/tables/final_research_direction/` and figures under `reports/figures/final_research_direction/`. This manuscript-prep analysis compiles PCA-centered spectral-biodiversity relationships, spectral metric pairwise relationships, biodiversity metric pairwise relationships, transformation diagnostics, metric driver screens, PCA computation-method checks, Moran's I spatial diagnostics, and variogram summaries.
- Current species-diversity versus phylogenetic-diversity concordance analysis is in `scripts/3_Analysis/species_phylogenetic_correlation_analysis.R`. The active report is `reports/analysis/20260817_species_phylogenetic_correlation_analysis.md`, with companion figures and tables under `reports/figures/species_phylogenetic_correlation/` and `reports/tables/species_phylogenetic_correlation/`. Scale-specific all-diversity metric scatterplot figures are named `05_all_diversity_metric_pairwise_scatter_10m.png`, `05_all_diversity_metric_pairwise_scatter_20m.png`, and `05_all_diversity_metric_pairwise_scatter_50m.png`. Edge-highlight versions for 10 m and 20 m are named `06_all_diversity_metric_pairwise_scatter_edge_highlight_10m.png` and `06_all_diversity_metric_pairwise_scatter_edge_highlight_20m.png`.
- Elevation-gradient all-diversity metric scatterplot figures are named `07_all_diversity_metric_pairwise_scatter_elevation_gradient_10m.png`, `07_all_diversity_metric_pairwise_scatter_elevation_gradient_20m.png`, and `07_all_diversity_metric_pairwise_scatter_elevation_gradient_50m.png`; their companion model table is `reports/tables/species_phylogenetic_correlation/elevation_adjusted_all_diversity_metric_models.csv`.
- Species composition types are now clustered from species presence/absence rather than abundance and named with diagnostic species codes. The composition-type heatmaps are `08_species_presence_composition_type_heatmap_10m.png`, `08_species_presence_composition_type_heatmap_20m.png`, and `08_species_presence_composition_type_heatmap_50m.png`; companion tables are `species_composition_type_names.csv`, `species_composition_type_presence_frequencies.csv`, and `species_composition_cluster_divergence_summary_presence_based.csv`.
- Composition-cluster-highlighted all-diversity scatterplot figures are named `09_all_diversity_metric_pairwise_scatter_composition_cluster_10m.png`, `09_all_diversity_metric_pairwise_scatter_composition_cluster_20m.png`, and `09_all_diversity_metric_pairwise_scatter_composition_cluster_50m.png`. Each panel legend shows each composition type with mean silhouette width calculated from that panel's two plotted diversity metrics; the companion table is `species_composition_scatterplot_silhouette.csv`.
- Crown-equivalent-individual ramp all-diversity scatterplot figures are named `10_all_diversity_metric_pairwise_scatter_individual_ramp_10m.png`, `10_all_diversity_metric_pairwise_scatter_individual_ramp_20m.png`, and `10_all_diversity_metric_pairwise_scatter_individual_ramp_50m.png`. Their first panels include the color-ramp legend, and the companion table is `quadrat_crown_equivalent_individual_totals.csv`.
- Species-count ramp all-diversity scatterplot figures are named `11_all_diversity_metric_pairwise_scatter_species_count_ramp_10m.png`, `11_all_diversity_metric_pairwise_scatter_species_count_ramp_20m.png`, and `11_all_diversity_metric_pairwise_scatter_species_count_ramp_50m.png`. Their first panels include the color-ramp legend, and the plotted values are stored as `present_species_count` in `quadrat_crown_equivalent_individual_totals.csv`.
- Species-count and crown-equivalent-individual moderation results are in `species_individual_moderated_metric_relationships.csv`, with BH-significant interaction terms summarized in `significant_species_individual_moderated_metric_relationships.csv`.
- Current spectral heterogeneity metric relationship analysis is in `scripts/3_Analysis/spectral_heterogeneity_relationship_analysis.R`. The active report is `reports/analysis/20260817_spectral_heterogeneity_relationship_analysis.md`, with companion tables under `reports/tables/spectral_heterogeneity_relationships/` and 33 figures under `reports/figures/spectral_heterogeneity_relationships/`. Each pairwise figure compares all combinations among 13 focal spectral heterogeneity measures: spectral angle entropy plus raw and standardized PCA alpha-hull area, spectral Rao's Q, mean Euclidean distance, convex hull, 3D hull volume, and 3D hull area metrics. The figure families are plain pairwise spectral metric scatterplots, elevation-gradient scatterplots, species-composition-highlighted scatterplots, crown-equivalent-individual ramp scatterplots, species-count ramp scatterplots, mean 563 nm pixel-brightness ramp scatterplots, region-specific brightness ramp scatterplots for blue 450-495 nm, green 500-570 nm, red 620-750 nm, and near-infrared 750-998 nm, and direct spectral-metric versus regional-illumination scatterplots.
- Current bootstrap metric-stat relationship analysis is in `scripts/3_Analysis/bootstrap_metric_stat_relationship_analysis.R`. The active report is `reports/analysis/20260817_bootstrap_metric_stat_relationship_analysis.md`, with companion tables under `reports/tables/bootstrap_metric_relationships/` and six scale-specific figures under `reports/figures/bootstrap_metric_relationships/`. The final-pipeline figure family plots spectral angle entropy against 70-iteration bootstrap SD, CV, and SE; the all-seven-metric figure family uses the fixed-4,000-pixel sample-size sensitivity branch and should be interpreted as metric stability diagnostics rather than final PCA uncertainty.
- Current SA bootstrap context relationship analysis is in `scripts/3_Analysis/sa_bootstrap_context_relationship_analysis.R`. The active report is `reports/analysis/20260817_sa_bootstrap_context_relationship_analysis.md`, with companion tables under `reports/tables/bootstrap_metric_relationships/` and six scale-specific figures under `reports/figures/bootstrap_metric_relationships/`. The continuous-context figures are named `03_sa_bootstrap_stats_vs_continuous_context_*m.png`; species-composition figures are named `04_sa_bootstrap_stats_vs_species_composition_*m.png`.
- Current paper workspace is `Documents/Paper/`. Manuscript-related outputs should be created there and should include the creation date in the filename. Current paper-prep files include `Documents/Paper/20260730_working_outline.md`, `Documents/Paper/20260730_paper_workspace_index.md`, `Documents/Paper/20260803_writing_style_notes.md`, `Documents/Paper/20260803_markdown_context_digest.md`, `Documents/Paper/20260803_remote_sensing_mdpi_author_requirements.md`, `Documents/Paper/20260803_detailed_manuscript_outline.docx`, `Documents/Paper/20260805_results_pca_mean_distance_afaith.docx`, `Documents/Paper/20260805_results_spectral_rao_q_afaith.docx`, `Documents/Paper/20260805_results_spectral_rao_q_afaith_centroid_clarified.docx`, `Documents/Paper/20260805_discussion_pca_distance_vs_rao_q_scale.docx`, `Documents/Paper/20260805_discussion_pca_distance_vs_rao_q_scale_centroid_clarified.docx`, and `Documents/Paper/SVH_Results.docx`.
- `Documents/Paper/SVH_Results.docx` is a user-created Word results document for manual final edits. Agents may read it for context but should not edit, overwrite, rename, or regenerate it unless explicitly instructed by the user.
- `Documents/Paper/20260811_methods_draft_spectral_phylogenetic_diversity.docx` is the current generated methods draft. It was built with `tools/create_methods_docx_and_tree_counts.R` and should be treated as a dated draft that can be regenerated or revised if the methods scope changes.
- `reports/figures/methods_tree_counts/` contains dated PNG charts of total eligible tree crowns by species for the full plot and the included 10 m, 20 m, and 50 m focal analysis quadrats.
- `reports/tables/methods_tree_counts/` contains the corresponding dated CSV tables, including individual all-plot and included-scale tables with `Genus`, `Species`, and `Number of individuals` headers plus a combined scale-comparison summary.
- Current edge-quadrat and bootstrap sensitivity outputs are under `reports/analysis/20260725_edge_bootstrap_sensitivity_sv_diversity.md`, `reports/figures/edge_bootstrap_sensitivity/`, and `reports/tables/edge_bootstrap_sensitivity/`. These are manuscript-relevant robustness/sensitivity sources, but the files are currently untracked in git.
- The species-diversity versus phylogenetic-diversity correlation pass was prepared in `reports/execution_plans/20260817_1019_species_phylogenetic_correlation_prep.md`, with prep task and validation notes in `reports/tasks/20260817_species_phylogenetic_correlation_prep.md` and `reports/validation/20260817_species_phylogenetic_correlation_prep_validation.md`.
- Current sample-size sensitivity outputs are under `reports/analysis/`, `reports/figures/sample_size_effects/`, and `reports/tables/sample_size_effects/`, including SA entropy, PCA mean distance, spectral Rao's Q, and alpha-hull area.
- Current 50 m PC1-PC3 mean-reflectance correlation outputs are under `reports/analysis/`, `reports/figures/pc1_mean_reflectance_correlation/`, and `reports/tables/pc1_mean_reflectance_correlation/`.
- Current PCA loading spectral-region outputs are under `reports/analysis/`, `reports/figures/pca_loading_spectral_regions/`, and `reports/tables/pca_loading_spectral_regions/`.
- `Quad_Spectra/10m_test`, `Quad_Spectra/20m_test`, and `Quad_Spectra/50m_test` are testing/validation artifacts and should not be treated as primary analytical inputs without explicit user direction.
