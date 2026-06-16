# Task Report: Quadrat Crop Reproducibility Attempt

Date: 2026-06-15

## Objective

Test whether existing spectra in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m` align with outputs produced by the crop logic in `scripts/1_Data Processing/HSI_quad_crop_refined.R`.

## Actions Performed

- Reviewed `HSI_quad_crop_refined.R` and confirmed the script's current checked-in parameters point to `VegIndex_NA_trimmed` and `Quad_Spectra/10m_VegIndex`, while its core crop/mask/resample/nearest-tile-selection logic is applicable to base spectra.
- Initially installed required R packages into temporary project-local `r_libs/`:
  - `terra`
  - `sf`
  - `dplyr`
  - `testthat`
- Added `tests/test_quad_crop_reproducibility.R`.
- Attempted to run the test with `C:/Program Files/R/R-4.2.3/bin/Rscript.exe`.

## Result

The test did not complete successfully.

Initial failures showed that most `HSI_NA_trimmed/` source rasters could not be opened by `terra`, likely because they were OneDrive cloud-only placeholders at the time of testing. Only these source rasters loaded successfully:

- `HSI_NA_trimmed/raw_0_rd_rf_or`
- `HSI_NA_trimmed/raw_11159_rd_rf_or`

The test was adjusted to skip unreadable source rasters and choose candidate quadrats that intersect readable rasters, but the final run was manually interrupted before completion.

## Files Created or Modified

- `tests/test_quad_crop_reproducibility.R`
- `reports/tasks/20260615_quad_crop_reproducibility_attempt.md`
- `reports/validation/20260615_quad_crop_reproducibility_validation.md`
- `reports/project_state.md`

## Library Location Update

The test was updated to use the standard R 4.2 library paths:

- `C:/Users/PaintRock/AppData/Local/R/win-library/4.2`
- `C:/Program Files/R/R-4.2.3/library`

The temporary project-local `r_libs/` directory was removed from the working directory.

## No Reprocessing Performed

`HSI_quad_crop_refined.R` was not rerun for any scale because the source rasters needed for a reliable validation were not fully readable locally.
