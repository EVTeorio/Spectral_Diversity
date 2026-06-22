# ExecPlan: Resample Smoothed Quadrat Spectra to 5 nm

Date: 2026-06-18

## Objective

Run the smoothed quadrat spectra through the resampling workflow and create new `_5nm` output folders for each scale.

## Requested Task

Use `scripts/1_Data Processing/resampling.R` on the smoothed spectra outputs from the previous processing step.

## Inputs

- `Quad_Spectra/10m_smooth`
- `Quad_Spectra/20m_smooth`
- `Quad_Spectra/50m_smooth`

## Outputs

- `Quad_Spectra/10m_smooth_5nm`
- `Quad_Spectra/20m_smooth_5nm`
- `Quad_Spectra/50m_smooth_5nm`

## Proposed Changes

- Update `resampling.R` to use the current working directory and standard R 4.2 user library path.
- Loop over the three smoothed scale folders.
- Resample spectral bands to 5 nm spacing using target bands `seq(398, 999, 5)`.
- Create output directories automatically.
- Skip completed outputs when both the ENVI data file and `.hdr` are present.
- Log progress to `logs/resampling_progress.log`.

## Validation Plan

- Confirm output raster counts match smoothed input raster counts.
- Confirm one sample raster per scale opens with `terra::rast()`.
- Confirm each sample output has the expected 121 bands.
- Update project documentation, task report, and validation report.

## Risks

- 50 m rasters are large and may require slower, cautious processing.
- OneDrive file sync or hydration can slow writes.
- Resampling changes the spectral band count from 328 to 121, so downstream scripts should use the `_5nm` folders deliberately.
