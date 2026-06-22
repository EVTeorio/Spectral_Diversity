# Task Report: Resample Smoothed Quadrat Spectra

Date: 2026-06-18

## Objective

Run the smoothed spectra through the resampling workflow and create new output folders with `_5nm` added to each smoothed scale folder name.

## Inputs

- `Quad_Spectra/10m_smooth`
- `Quad_Spectra/20m_smooth`
- `Quad_Spectra/50m_smooth`

## Outputs

- `Quad_Spectra/10m_smooth_5nm`
- `Quad_Spectra/20m_smooth_5nm`
- `Quad_Spectra/50m_smooth_5nm`

## Script Updated

- `scripts/1_Data Processing/resampling.R`

## Processing Notes

- Added the standard R 4.2 user library path before loading packages.
- Updated the script to loop over `10m_smooth`, `20m_smooth`, and `50m_smooth`.
- Created output folders automatically.
- Skipped existing completed outputs so interrupted runs can resume.
- Added progress logging to `logs/resampling_progress.log`.
- Resampled to target bands from 398 nm through 998 nm in 5 nm steps.
- Used scale-aware `terra::app()` core counts:
  - 10 m smoothed: 8 cores
  - 20 m smoothed: 4 cores
  - 50 m smoothed: 2 cores

## Results

| Scale | Output Folder | Raster Count |
|---|---|---:|
| 10 m | `Quad_Spectra/10m_smooth_5nm` | 1,909 |
| 20 m | `Quad_Spectra/20m_smooth_5nm` | 485 |
| 50 m | `Quad_Spectra/50m_smooth_5nm` | 80 |

## Notes

- No source spectra or smoothed spectra were modified.
- The current implementation uses band-wise interpolation in `terra::app()` to avoid loading full rasters into R data frames.
