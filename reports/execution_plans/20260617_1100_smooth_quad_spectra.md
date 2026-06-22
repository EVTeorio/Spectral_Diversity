# ExecPlan: Smooth Confirmed Quadrat Spectra

Date: 2026-06-17

## Objective

Apply the Savitzky-Golay smoothing workflow to the confirmed partitioned quadrat spectra at 10 m, 20 m, and 50 m scales.

## Requested Task

Use `scripts/1_Data Processing/smoothing.R` on each confirmed quad-scale directory and create output folders with the `_smooth` suffix.

## Files To Review

- `CODEX_AGENT_GUIDELINES.md`
- `scripts/1_Data Processing/smoothing.R`
- `scripts/1_Data Processing/HSI_quad_crop_refined.R`
- `reports/project_state.md`
- `reports/data_dictionary.md`
- `reports/directory_map.md`

## Relevant Inputs

- `Quad_Spectra/10m`
- `Quad_Spectra/20m`
- `Quad_Spectra/50m`

## Proposed Outputs

- `Quad_Spectra/10m_smooth`
- `Quad_Spectra/20m_smooth`
- `Quad_Spectra/50m_smooth`

## Proposed Changes

- Update `smoothing.R` to loop over the three confirmed scale folders.
- Preserve the current Savitzky-Golay parameters unless execution reveals they are invalid for the input band structure.
- Skip ENVI/GDAL sidecar files.
- Use `overwrite = FALSE` so interrupted runs can resume without rewriting completed rasters.
- Create output directories when absent.

## Validation Plan

- Confirm output folders are created.
- Confirm output raster counts match the input raster counts for each scale.
- Confirm at least one output raster per scale can be opened by `terra::rast()`.
- Update task and validation reports after processing.

## Risks

- Processing may take a long time, especially for 50 m rasters.
- OneDrive cloud-only files may need hydration before R can read them.
- `terra`, `signal`, or `beepr` may be unavailable or have package/library issues.
