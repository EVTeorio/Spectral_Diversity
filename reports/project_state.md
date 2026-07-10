# Project State

Last updated: 2026-07-09

## Current Objective

Build a reproducible and well-documented research workflow for evaluating relationships between hyperspectral spectral heterogeneity and biodiversity/phylogenetic diversity metrics in the Paint Rock Forest Dynamics Plot.

## Active Research Questions

- Does spectral heterogeneity increase with Faith's phylogenetic diversity?
- Do relationships differ across 10 m, 20 m, and 50 m spatial grains?
- Do environmental covariates such as elevation, slope, aspect, or topographic exposure explain additional variation in spectral heterogeneity?
- Which spectral heterogeneity metric is most useful for the biodiversity relationship: spectral angle entropy, band entropy, PCA-space metrics, or another metric?

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
- Added `scripts/3_Analysis/multiscale_spectral_biodiversity_analysis.R` to evaluate spectral-biodiversity relationships across 10 m, 20 m, and 50 m quadrat scales using the current combined analysis tables.
- Created user-facing PDF reports under `Documents/PDFs/`:
  - `spectral_biodiversity_multiscale_findings.pdf`
  - `spectral_biodiversity_model_appendix.pdf`
- Created analysis figures under `reports/figures/multiscale_spectral_biodiversity/` and analysis tables under `reports/tables/multiscale_spectral_biodiversity/`.
- The multiscale analysis treats `spec_sa` as the primary spectral heterogeneity response and uses `spec_pca_mean`, `spec_rao_q`, `spec_alpha`, `spec_spca_mean`, `spec_spca_rao`, and `spec_spca_alpha` as secondary response metrics. The script has been updated for standardized PCA fields, but biodiversity/environment relationship analyses were intentionally not rerun during the PCA-only diagnostic update.
- The multiscale analysis uses `phy_faith`, `phy_afaith`, and `sp_shannon` as primary biodiversity predictors and `env_elev` plus `env_tri11` as environmental controls.
- Documented 10 m and 20 m edge quadrats were excluded from primary inference in the multiscale analysis; 50 m used all quadrats because no separate 50 m edge rule is documented.
- Best-supported primary-response candidate models by AIC used abundance-weighted Faith's PD: biodiversity plus environment at 10 m and 20 m, and abundance-weighted Faith's PD by elevation interaction at 50 m. Adjusted R2 values were approximately 0.038, 0.094, and 0.242, respectively.
- Residual Moran's I remained significant in primary spectral angle entropy models at all three scales, so OLS p-values should be interpreted cautiously until spatial model sensitivity checks are added.
- Added `reports/tasks/20260624_multiscale_spectral_biodiversity_analysis.md`.
- Added `reports/validation/20260624_multiscale_spectral_biodiversity_analysis_validation.md`.
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
- Add spatial model sensitivity checks for the multiscale spectral-biodiversity analysis because residual Moran's I remained significant after the current OLS candidate models.
- Join bootstrap quality-control fields (`boot_sd`, `boot_cv`, and `method`) into a companion modeling table and rerun sensitivity checks that flag or exclude high-CV quadrats.

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
2. Decide whether 10 m and 20 m missing spectral quadrats should remain `NA`, be excluded, or be tracked in a missing-quadrat report.
3. Decide whether bootstrap quality-control fields such as `boot_sd`, `boot_cv`, and `method` should be added to a separate QC companion table rather than the value-only combined tables.
4. Run downstream model sensitivity checks with high-CV quadrats flagged or excluded, using 5% and 10% bootstrap-CV thresholds.
5. Compare the new current 20 m spectral heterogeneity output against older 20 m spectral entropy products before replacing downstream analysis columns.
6. Use `C:/Program Files/R/R-4.2.3/bin/Rscript.exe` for R execution unless the shell path is updated later.
7. Extend the multiscale PDF analysis with spatial GLS/GAM or spatial eigenvector sensitivity models before treating coefficient p-values as final manuscript inference.
