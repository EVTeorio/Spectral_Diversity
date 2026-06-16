# Project State

Last updated: 2026-06-15

## Current Objective

Build a reproducible and well-documented research workflow for evaluating relationships between hyperspectral spectral heterogeneity and biodiversity/phylogenetic diversity metrics in the Paint Rock Forest Dynamics Plot.

## Active Research Questions

- Does spectral heterogeneity increase with Faith's phylogenetic diversity?
- Do relationships differ across 10 m, 20 m, and 50 m spatial grains?
- Do environmental covariates such as elevation, slope, aspect, or topographic exposure explain additional variation in spectral heterogeneity?
- Which spectral heterogeneity metric is most useful for the biodiversity relationship: spectral angle entropy, band entropy, PCA-space metrics, or another metric?

## Completed Work

- Confirmed `CODEX_AGENT_GUIDELINES.md` and `RESEARCH_OBJECTIVES.md` as governing/context documents.
- Created required governance directories under `reports/` and `logs/`.
- Created an execution plan for repository documentation and cleanup.
- Created baseline `reports/directory_map.md`.
- Created baseline `reports/data_dictionary.md`.
- Removed obsolete `.Rhistory` session files from the repository.
- Updated `CODEX_AGENT_GUIDELINES.md` with R/Rtools access guidance, subsampled unit testing expectations, and current data assumptions.
- Removed authorized legacy directories named `old`, `Outdated`, and `Currently not relevant`.
- Created `tests/` for future `testthat` unit tests.
- Added `tests/test_quad_crop_reproducibility.R` to compare a small quadrat subset against the core crop logic from `HSI_quad_crop_refined.R`.
- Updated R package library use to the standard R 4.2 paths:
  - `C:/Users/PaintRock/AppData/Local/R/win-library/4.2`
  - `C:/Program Files/R/R-4.2.3/library`

## Pending Work

- Review current R scripts for package dependencies, inputs, outputs, assumptions, and reproducibility risks.
- Standardize script and directory names that contain spaces or spelling errors after checking dependencies.
- Add formal `testthat` tests once core functions are modularized.
- Expand shapefile and raster data dictionary entries using R geospatial packages.
- Revisit the quadrat crop reproducibility test after source HSI rasters are hydrated locally from OneDrive.

## Known Issues

- R/Rscript was not available on the command path during this documentation pass.
- Several script and directory names contain spaces or misspellings, which can complicate reproducible command-line execution.
- Some current CSV files contain unnamed first columns.
- `Indices_SHPs/20m_spectral_sp.csv` has 64 missing values for each primary spectral metric.
- `Rscript` is available at `C:/Program Files/R/R-4.2.3/bin/Rscript.exe`, matching the required R 4.2.3 compatibility target, though it is not available as `Rscript` on the active shell path.
- Quadrat crop reproducibility testing could not be completed because most `HSI_NA_trimmed/` source rasters were unreadable/cloud-only OneDrive files during the test attempt. Only `raw_0_rd_rf_or` and `raw_11159_rd_rf_or` loaded successfully.
- The final crop reproducibility test run was manually interrupted before completion.

## Technical Debt

- Some scripts appear exploratory or duplicated, including `modeling_test.R` and `scratch paper.R`.
- Current workflows likely use hard-coded paths and should be audited.
- The project now has a `tests/` directory and an initial crop reproducibility test, but that test has not passed end-to-end yet.
- The temporary project-local `r_libs/` package library was removed from the working directory.
- The user R library has a broken/locked `rlang` installation: `dplyr` and `testthat` fail to load because the `rlang` DLL version does not match the package metadata. This likely requires closing active R/RStudio sessions and reinstalling `rlang`.

## Important Assumptions

- Legacy directories named `old`, `Outdated`, and `Currently not relevant` have been removed from the active repository.
- The primary integrated 20 m analysis table is `Indices_SHPs/20m_spectral_sp.csv`.
- Raw spectral inputs are in `HSI_NA_trimmed/`.
- Quadrat spectral products are in `Quad_Spectra/`.
- Quadrat boundary files are in `Quad_Scale_SHPs/`.
- Spectral products that produce spectral heterogeneity values are in `Indices_SHPs/Spectral_diversitySHPs/`.
- Unnamed first CSV columns are likely index columns.

## Next Recommended Actions

1. Review `scripts/1_Data Processing/HSI_quad_crop_refined.R` and all scripts that feed `Indices_SHPs/20m_spectral_sp.csv`.
2. Use `C:/Program Files/R/R-4.2.3/bin/Rscript.exe` for R execution unless the shell path is updated later.
3. Rename or replace files with spaces and misspellings after dependency references are checked.
4. Begin tests for diversity metric joins and missing-quadrat handling using small representative subsamples.
5. Expand the data dictionary with shapefile attribute counts and raster metadata once R geospatial tooling is available.
6. Close active R/RStudio sessions if needed, repair `rlang` in `C:/Users/PaintRock/AppData/Local/R/win-library/4.2`, hydrate the `HSI_NA_trimmed/` source rasters locally, then rerun `tests/test_quad_crop_reproducibility.R`.
