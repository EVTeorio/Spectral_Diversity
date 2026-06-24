# Validation Report: Environmental Variables at All Scales

## Validation Scope

Validate the all-scale environmental workflow that creates mean elevation and mean Riley TRI values from `PRFPD_DTM_leafOff.tiff`.

## Checks Performed

- Confirmed required R packages are available in the R 4.2 user library:
  - `terra`
  - `sf`
  - `dplyr`
  - `readr`
  - `testthat`
- Inspected the DTM raster:
  - 694 rows by 700 columns
  - 1 m resolution
  - elevation range approximately 192.15 m to 303.91 m
  - CRS: NAD83(CSRS) / UTM zone 16N with CGVD2013 height
- Confirmed quadrat shapefile row counts:
  - 10 m: 2,000 features
  - 20 m: 500 features
  - 50 m: 80 features
- Ran `tests/test_environmental_variables_all_scales.R`.
- Ran `scripts/2_Indices Creation/Enviro_Variables/environmental_variables_all_scales.R`.
- Read generated CSVs and shapefiles back into R.
- Confirmed all output row counts match the source quadrat shapefiles.
- Confirmed no missing values in `elev_mean`, `tri5_mean`, `tri11_mean`, or `tri21_mean`.

## Row Counts

| Scale | CSV Rows | Shapefile Features |
|---|---:|---:|
| 10 m | 2,000 | 2,000 |
| 20 m | 500 | 500 |
| 50 m | 80 | 80 |

## Output Value Ranges

| Scale | Metric | Minimum | Maximum |
|---|---|---:|---:|
| 10 m | `elev_mean` | 200.0 | 291.9 |
| 10 m | `tri5_mean` | 0.2799 | 6.5791 |
| 10 m | `tri11_mean` | 1.102 | 27.506 |
| 10 m | `tri21_mean` | 3.807 | 78.215 |
| 20 m | `elev_mean` | 201.9 | 289.3 |
| 20 m | `tri5_mean` | 0.3191 | 4.2006 |
| 20 m | `tri11_mean` | 1.174 | 18.669 |
| 20 m | `tri21_mean` | 5.088 | 61.070 |
| 50 m | `elev_mean` | 206.1 | 282.6 |
| 50 m | `tri5_mean` | 0.9169 | 2.9599 |
| 50 m | `tri11_mean` | 3.934 | 13.441 |
| 50 m | `tri21_mean` | 13.19 | 47.79 |

## Result

Validation passed. The current environmental variable outputs are ready for downstream joins using `scale` plus `quad_id`.
