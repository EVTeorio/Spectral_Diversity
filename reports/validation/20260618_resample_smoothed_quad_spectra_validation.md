# Validation Report: Resample Smoothed Quadrat Spectra

Date: 2026-06-18

## Validation Scope

Validate the `_smooth_5nm` raster outputs generated from current smoothed quadrat spectra.

## Count Validation

| Scale | Input Raster Count | Output Raster Count | Result |
|---|---:|---:|---|
| 10 m | 1,909 | 1,909 | Pass |
| 20 m | 485 | 485 | Pass |
| 50 m | 80 | 80 | Pass |

## Raster Open Validation

Sample outputs were opened with `terra::rast()` using R 4.2.3 and `terra` 1.7.71.

| Output Folder | Sample Raster | Bands | Rows | Columns | First Band | Last Band | Result |
|---|---|---:|---:|---:|---|---|---|
| `Quad_Spectra/10m_smooth_5nm` | `0_a` | 121 | 136 | 136 | `398 nm` | `998 nm` | Pass |
| `Quad_Spectra/20m_smooth_5nm` | `0` | 121 | 272 | 272 | `398 nm` | `998 nm` | Pass |
| `Quad_Spectra/50m_smooth_5nm` | `sub50_1` | 121 | 680 | 680 | `398 nm` | `998 nm` | Pass |

## Result

Resampling completed successfully for all requested smoothed quad-scale directories. Output counts match input counts, and representative rasters from each output folder are readable with the expected 121-band 5 nm structure.
