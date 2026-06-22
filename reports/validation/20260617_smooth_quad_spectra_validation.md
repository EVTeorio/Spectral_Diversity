# Validation Report: Smooth Quadrat Spectra

Date: 2026-06-17

## Validation Scope

Validate the `_smooth` raster outputs generated from the confirmed partitioned quadrat spectra.

## Count Validation

| Scale | Input Raster Count | Output Raster Count | Result |
|---|---:|---:|---|
| 10 m | 1,909 | 1,909 | Pass |
| 20 m | 485 | 485 | Pass |
| 50 m | 80 | 80 | Pass |

## Raster Open Validation

Sample outputs were opened with `terra::rast()` using R 4.2.3 and `terra` 1.7.71.

| Output Folder | Sample Raster | Bands | Rows | Columns | Result |
|---|---|---:|---:|---:|---|
| `Quad_Spectra/10m_smooth` | `0_a` | 328 | 136 | 136 | Pass |
| `Quad_Spectra/20m_smooth` | `0` | 328 | 272 | 272 | Pass |
| `Quad_Spectra/50m_smooth` | `sub50_1` | 328 | 680 | 680 | Pass |

## Result

Smoothing completed successfully for all requested quad-scale directories. Output counts match input counts, and representative rasters from each output folder are readable.
