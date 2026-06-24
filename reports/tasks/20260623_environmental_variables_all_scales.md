# Task Report: Environmental Variables at All Scales

## Objective

Create mean elevation and mean Riley topographic roughness index variables for 10 m, 20 m, and 50 m quadrats from `PRFPD_DTM_leafOff.tiff`.

## Script Added

- `scripts/2_Indices Creation/Enviro_Variables/environmental_variables_all_scales.R`

## Key Changes

- Added an all-scale environmental workflow using `terra`, `sf`, `dplyr`, and `readr`.
- Transforms quadrat shapefiles from WGS84 to the DTM CRS before raster extraction.
- Calculates mean elevation per quadrat.
- Calculates Riley TRI rasters from the DTM using 5x5, 11x11, and 21x21 windows.
- Uses compact output field names that are shapefile-safe:
  - `elev_mean`
  - `tri5_mean`
  - `tri11_mean`
  - `tri21_mean`

## Outputs Created

- `Indices_SHPs/Enviro_SHPs/enviro_variables_10m.csv`
- `Indices_SHPs/Enviro_SHPs/enviro_variables_10m.shp`
- `Indices_SHPs/Enviro_SHPs/enviro_variables_20m.csv`
- `Indices_SHPs/Enviro_SHPs/enviro_variables_20m.shp`
- `Indices_SHPs/Enviro_SHPs/enviro_variables_50m.csv`
- `Indices_SHPs/Enviro_SHPs/enviro_variables_50m.shp`

## Tests Added

- `tests/test_environmental_variables_all_scales.R`

The tests cover Riley TRI helper calculations, missing-value handling, invalid window-size rejection, raster output consistency, and rolling-window sum correctness.

## Notes

- The original `scripts/2_Indices Creation/Enviro_Variables/elevation.R` was left unchanged as a historical single-scale reference.
- The first custom-focal TRI implementation was too slow for the all-scale workflow. It was replaced with an in-memory rolling-sum implementation using the identity behind Riley TRI.
