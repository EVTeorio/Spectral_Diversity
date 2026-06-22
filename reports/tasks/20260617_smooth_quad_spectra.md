# Task Report: Smooth Quadrat Spectra

Date: 2026-06-17

## Objective

Apply the smoothing workflow to each confirmed partitioned quadrat spectra scale and write new output folders with the `_smooth` suffix.

## Inputs

- `Quad_Spectra/10m`
- `Quad_Spectra/20m`
- `Quad_Spectra/50m`

## Outputs

- `Quad_Spectra/10m_smooth`
- `Quad_Spectra/20m_smooth`
- `Quad_Spectra/50m_smooth`

## Script Updated

- `scripts/1_Data Processing/smoothing.R`

## Processing Notes

- Added the standard R 4.2 user library path before loading packages.
- Preserved Savitzky-Golay smoothing parameters:
  - polynomial order: 3
  - window size: 15
- Updated the script to loop over `10m`, `20m`, and `50m`.
- Created output folders automatically.
- Skipped existing completed outputs so interrupted runs could resume.
- Added progress logging to `logs/smoothing_progress.log`.
- Used scale-aware `terra::app()` core counts:
  - 10 m: 8 cores
  - 20 m: 4 cores
  - 50 m: 2 cores

## Results

| Scale | Output Folder | Raster Count |
|---|---|---:|
| 10 m | `Quad_Spectra/10m_smooth` | 1,909 |
| 20 m | `Quad_Spectra/20m_smooth` | 485 |
| 50 m | `Quad_Spectra/50m_smooth` | 80 |

## Notes

- The initial single-worker run was stopped and resumed with scale-aware parallel processing.
- No source spectra were modified.
