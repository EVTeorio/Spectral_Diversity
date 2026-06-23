# Project State

Last updated: 2026-06-22

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
  - `Indices_SHPs/10m_SA_entropy_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/20m_SA_entropy_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/50m_SA_entropy_smooth_masked_5nm_summary.csv`
  - matching long and wide bootstrap CSVs for 10 m, 20 m, and 50 m
  - `Indices_SHPs/Spectral_diversitySHPs/SA_entropy_10m_smooth_masked_5nm.shp`
  - `Indices_SHPs/Spectral_diversitySHPs/SA_entropy_20m_smooth_masked_5nm.shp`
  - `Indices_SHPs/Spectral_diversitySHPs/SA_entropy_50m_smooth_masked_5nm.shp`
- Added `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R` to create all-scale spectral heterogeneity outputs with four primary measures:
  - existing spectral angle entropy
  - PCA-space mean Euclidean distance from the quadrat centroid
  - alpha-hull area in global PC1-PC2 space
  - Rao's Q in global PC1-PC3 space using equal pixel weights and squared Euclidean distance
- Rebuilt the global PCA basis after applying manual atmospheric/cloud exclusions from `RESEARCH_OBJECTIVES.md`. The current PCA basis excludes affected quadrats from PCA sampling and from PCA-dependent metric calculation.
- Current PCA-dependent outputs supersede and replace the earlier 2026-06-22 PCA-dependent outputs created before manual exclusions were applied; the earlier pre-exclusion PCA values should be disregarded.
- Created a global PCA basis from a uniform per-eligible-quadrat sample across current 10 m, 20 m, and 50 m smoothed 5 nm spectra, with a requested target of up to 500 retained sunlit pixels per eligible quadrat raster and a configured maximum PCA sample of 800,000 rows. The effective uniform sample was 354 pixels per eligible raster, producing 797,916 sampled pixels from 2,254 eligible rasters.
- Created PCA variance diagnostic outputs:
  - `Indices_SHPs/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`
  - `Indices_SHPs/Spectral_diversitySHPs/global_pca_smooth_masked_5nm_variance_explained.csv`
  - `reports/figures/spectral_heterogeneity/global_pca_smooth_masked_5nm_variance_explained.png`
- Global PCA variance explained after manual exclusions: PC1 = 66.98%, PC2 = 19.93%, PC3 = 3.80%, cumulative PC1-PC3 = 90.71%.
- Created all-scale combined spectral heterogeneity CSVs and shapefiles:
  - `Indices_SHPs/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/20m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/50m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_10m_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_20m_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_50m_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_10m_smooth_masked_5nm.shp`
  - `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_20m_smooth_masked_5nm.shp`
  - `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_50m_smooth_masked_5nm.shp`
- Created PCA-only metric CSVs:
  - `Indices_SHPs/10m_PCA_metrics_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/20m_PCA_metrics_smooth_masked_5nm_summary.csv`
  - `Indices_SHPs/50m_PCA_metrics_smooth_masked_5nm_summary.csv`
- Added `reports/validation/20260622_spectral_heterogeneity_all_metrics_validation.md`.
- Created a bootstrap variation QC analysis comparing within-quadrat bootstrap variation to between-quadrat spectral heterogeneity variation.
- Added `scripts/3_Analysis/bootstrap_variation_analysis.R`.
- Created `reports/analysis/20260618_bootstrap_variation_analysis.md` with seven figures and three diagnostic tables under `reports/figures/bootstrap_variation/` and `reports/tables/bootstrap_variation/`.
- Added `reports/figures/bootstrap_variation/spectral_entropy_histograms_by_scale.png` to visualize spectral heterogeneity distribution shape and skewness by scale.
- Bootstrap variation analysis conclusion: `spectral_entropy` / `boot_mean` values are usable as primary heterogeneity estimates, but `boot_sd`, `boot_cv`, and `method` should be carried forward as quality-control fields and high-CV quadrats should be used in sensitivity checks.
- Refactored the plant-diversity workflow into `scripts/2_Indices Creation/Plant_diversity/plant_diversity_all_scales.R`.
- Updated the older plant-diversity entry scripts (`sp_weighted_matrix.R`, `species_diversity.R`, and `phylogenetic_diversity.R`) to call the all-scale workflow.
- Created all-scale plant diversity outputs under `Indices_SHPs/Diversity_SHPs/`:
  - `plant_diversity_10m.csv` and `plant_diversity_10m.shp`
  - `plant_diversity_20m.csv` and `plant_diversity_20m.shp`
  - `plant_diversity_50m.csv` and `plant_diversity_50m.shp`
- Confirmed plant-diversity `quad_id` values align with current spectral diversity IDs: 10 m uses `sub_id` values such as `0_a`, 20 m uses `Name`, and 50 m uses `Name` values such as `sub50_1`.
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
- Join the new per-scale spectral heterogeneity values into downstream biodiversity/species/environment analysis tables using `quad_id`.
- Decide how to handle missing 10 m and 20 m quadrats that do not have current raster summaries.
- Carry bootstrap quality-control fields into downstream tables, especially `boot_sd`, `boot_cv`, and `method`.
- Run downstream model sensitivity checks with high-CV quadrats flagged or excluded.
- Join the new all-scale combined spectral heterogeneity outputs to plant-diversity outputs by `scale` plus `quad_id`.
- When modifying R workflows, first review nearby scripts in `scripts/` to match the project's existing loop structure and file-processing pattern.
- Standardize script and directory names that contain spaces or spelling errors after checking dependencies.
- Add formal `testthat` tests once core functions are modularized.
- Expand shapefile and raster data dictionary entries using R geospatial packages.

## Known Issues

- R/Rscript was not available on the command path during this documentation pass.
- Several script and directory names contain spaces or misspellings, which can complicate reproducible command-line execution.
- Some current CSV files contain unnamed first columns.
- `Indices_SHPs/20m_spectral_sp.csv` has 64 missing values for each primary spectral metric.
- `Rscript` is available at `C:/Program Files/R/R-4.2.3/bin/Rscript.exe`, matching the required R 4.2.3 compatibility target, though it is not available as `Rscript` on the active shell path.
- Codex can access the R 4.2 user package library under `C:/Users/PaintRock/AppData/Local/R/win-library/4.2` when commands are run with permission outside the workspace sandbox; sandboxed Rscript may still report only the base R library.

## Technical Debt

- Some scripts appear exploratory or duplicated, including `modeling_test.R` and `scratch paper.R`.
- Current workflows likely use hard-coded paths and should be audited.
- The temporary project-local `r_libs/` package library was removed from the working directory.

## Important Assumptions

- Legacy directories named `old`, `Outdated`, and `Currently not relevant` have been removed from the active repository.
- The primary integrated 20 m analysis table is `Indices_SHPs/20m_spectral_sp.csv`.
- The current all-scale plant-diversity outputs are `Indices_SHPs/Diversity_SHPs/plant_diversity_10m.csv`, `plant_diversity_20m.csv`, and `plant_diversity_50m.csv`, with matching shapefiles.
- Raw spectral inputs are in `HSI_NA_trimmed/`.
- The user has independently tested and confirmed the partitioned quadrat spectra as the inputs to use moving forward.
- Confirmed partitioned quadrat spectra are in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m`.
- Current smoothed spectra derived from those confirmed inputs are in `Quad_Spectra/10m_smooth`, `Quad_Spectra/20m_smooth`, and `Quad_Spectra/50m_smooth`.
- Current smoothed 5 nm spectra are in `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, and `Quad_Spectra/50m_smooth_5nm`.
- The current spectral heterogeneity workflow should consume the current `_smooth_5nm` spectra, not older `_resampled_5nm` or `_smoothed_5nm` folders.
- Spectral angle entropy values should be produced per quadrat at each resolution; 10 m uses `sub_id` IDs such as `0_a`, 20 m uses numeric `Name` IDs, and 50 m uses `Name` IDs such as `sub50_1`.
- Plant-diversity outputs now use the same `quad_id` convention as spectral diversity outputs, so downstream joins should use `scale` plus `quad_id`.
- Plant-diversity species matrix values represent summed crown-overlap proportions by species and quadrat, preserving the existing crown-buffer/intersection logic.
- Exact all-pixel spectral angle entropy is only reasonable for small masked quadrats because pairwise angle counts scale quadratically with pixel count; the workflow should use bootstrap `boot_mean` as the primary fallback value for larger quadrats.
- Current spectral heterogeneity outputs used 70 bootstrap replicates, 5,000 sampled pixels per bootstrap, and 10,000 sampled pixel-pairs per large bootstrap replicate.
- Validation found 10 m summary rows for 1,909 rasters, with 1,866 bootstrap-mean estimates, 37 exact all-pixel estimates, and 6 insufficient-pixel outputs.
- Validation found 20 m summary rows for 485 rasters, with 477 bootstrap-mean estimates and 8 exact all-pixel estimates.
- Validation found 50 m summary rows for 80 rasters, with 79 bootstrap-mean estimates and 1 exact all-pixel estimate.
- The 10 m output shapefile has 97 missing `spec_ent` values: 91 shapefile quadrats without matching current raster outputs and 6 insufficient-pixel raster summaries.
- The 20 m output shapefile has 15 missing `spec_ent` values from shapefile quadrats without matching current raster outputs.
- The 50 m output shapefile has no missing `spec_ent` values.
- The current PCA-dependent spectral heterogeneity workflow applies the manual atmospheric/cloud exclusion policy `manual_atmospheric_cloud_exclusions_20260622`. Current raster-output exclusions are 165 of 1,909 at 10 m, 49 of 485 at 20 m, and 6 of 80 at 50 m.
- The current 50 m manual atmospheric/cloud exclusions are `sub50_80`, `sub50_79`, `sub50_71`, `sub50_70`, `sub50_62`, and `sub50_53`.
- The combined 10 m spectral heterogeneity CSV has 1,909 rows. It has 6 missing SA entropy values from insufficient-pixel SA entropy summaries, and 165 missing PCA Euclidean, Rao Q, and 3D PCA hull values from manual exclusions. It has 168 missing alpha-hull values: 165 manual exclusions plus 3 alpha-hull failures after deterministic sampled fallback.
- The combined 20 m spectral heterogeneity CSV has 485 rows. It has no missing SA entropy values and 49 missing PCA-dependent values from manual exclusions. The combined 20 m shapefile has 15 missing SA entropy values from shapefile quadrats without matching current raster outputs, and 64 missing PCA-dependent values from those unmatched polygons plus 49 manual exclusions.
- The combined 50 m spectral heterogeneity CSV and shapefile have 80 rows/features. They have no missing SA entropy values and 6 missing PCA-dependent values from manual exclusions.
- The supplemental three-axis hull output is a convex hull volume/area in global PC1-PC3 space, not a true 3D alpha hull, because `alphashape3d` is not currently installed.
- Bootstrap QC found low median CV at all scales, about 0.4% to 0.5%, but a non-trivial high-CV tail.
- Using a 5% bootstrap-CV flag, 372 of 1,866 10 m bootstrap-estimated quadrats, 60 of 477 20 m quadrats, and 11 of 79 50 m quadrats would be flagged.
- Using a 10% bootstrap-CV flag, 127 of 1,866 10 m bootstrap-estimated quadrats and 12 of 477 20 m quadrats would be flagged; no 50 m quadrats exceed 10%.
- Spectral heterogeneity distributions show mild right skew at 10 m, near-symmetry at 20 m, and stronger right skew at 50 m.
- `Quad_Spectra/10m_test`, `Quad_Spectra/20m_test`, and `Quad_Spectra/50m_test` are testing/validation artifacts unless the user explicitly promotes them.
- Quadrat boundary files are in `Quad_Scale_SHPs/`.
- Spectral products that produce spectral heterogeneity values are in `Indices_SHPs/Spectral_diversitySHPs/`.
- Unnamed first CSV columns are likely index columns.

## Next Recommended Actions

1. Join the new `*_SA_entropy_smooth_masked_5nm_summary.csv` values and bootstrap QC fields into the new all-scale `plant_diversity_*m.csv` tables by `quad_id`.
2. Decide whether 10 m and 20 m missing quadrats should remain `NA`, be excluded, or be tracked in a missing-quadrat report.
3. Run downstream model sensitivity checks with high-CV quadrats flagged or excluded, using 5% and 10% bootstrap-CV thresholds.
4. Compare the new current 20 m spectral heterogeneity output against older 20 m spectral entropy products before replacing downstream analysis columns.
5. Use `C:/Program Files/R/R-4.2.3/bin/Rscript.exe` for R execution unless the shell path is updated later.
