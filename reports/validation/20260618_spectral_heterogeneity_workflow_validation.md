# Validation Report: Spectral Heterogeneity Workflow Update

Date: 2026-06-18

## Validation Scope

Validate the refactored spectral angle entropy bootstrap workflow and completed 10 m, 20 m, and 50 m spectral heterogeneity outputs from current smoothed 5 nm spectra.

## Package Availability

The user confirmed from RStudio:

- `terra 1.7.71`
- `C:/Users/PaintRock/AppData/Local/R/win-library/4.2/terra`

Codex also confirmed package visibility when Rscript was allowed to read the same user-library path:

- `terra`: available
- `beepr`: available

The earlier package issue was due to sandboxed Rscript visibility, not a missing R package.

## Workflow Settings

- Inputs: `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, `Quad_Spectra/50m_smooth_5nm`
- Shadow mask band: `563 nm`
- Shadow mask threshold: `0.0305476`
- Retained pixels: reflectance greater than the threshold under the current `DIRECTION <- "<"` setting
- Bootstrap replicates: 70
- Pixel subsample per bootstrap: 5,000
- Pair sample per bootstrap when exact bootstrap pairs are too large: 10,000
- Exact all-pixel entropy threshold: 2,000,000 possible pairs

## Checks Performed

| Check | Result | Notes |
|---|---|---|
| Markdown context reviewed | Pass | Reviewed governing, objective, state, data dictionary, directory map, task, execution-plan, and validation Markdown files. |
| Current input stage confirmed | Pass | Workflow starts from `Quad_Spectra/*_smooth_5nm`. |
| Shadow mask documented | Pass | Script carries forward `563 nm`, threshold `0.0305476`, and the retained sunlit condition used by the current direction setting. |
| Script parse | Pass | `SA_entropy_bootstrapping.R` and `run_sa_entropy_scale.R` parse successfully under R 4.2.3. |
| Helper tests | Pass | `tests/test_sa_entropy_bootstrapping.R` passed using base R assertions. |
| Smoke test | Pass | `tests/smoke_sa_entropy_bootstrapping.R` processed one raster per scale with reduced bootstrap settings. |
| Full raster execution | Pass | 50 m, 20 m, and 10 m workflows completed and wrote CSV plus shapefile outputs. |
| Bootstrap file dimensions | Pass | Long-format bootstrap rows equal `summary_rows * 70`; wide-format outputs have 70 bootstrap columns. |
| Completion beep | Pass | `library(beepr); beep()` completed using the R 4.2 user library. |

## Output Summary

| Scale | Summary Rows | Method Counts | Spectral Entropy Missing | Shapefile Features | Shapefile Missing `spec_ent` |
|---|---:|---|---:|---:|---:|
| 10 m | 1,909 | 1,866 `bootstrap_mean`; 37 `exact_all_pixels`; 6 `insufficient_pixels` | 6 | 2,000 | 97 |
| 20 m | 485 | 477 `bootstrap_mean`; 8 `exact_all_pixels` | 0 | 500 | 15 |
| 50 m | 80 | 79 `bootstrap_mean`; 1 `exact_all_pixels` | 0 | 80 | 0 |

## Bootstrap Output Dimensions

| Scale | Long Rows | Wide Rows | Bootstrap Columns |
|---|---:|---:|---:|
| 10 m | 133,630 | 1,909 | 70 |
| 20 m | 33,950 | 485 | 70 |
| 50 m | 5,600 | 80 | 70 |

## Output Files

- `Quad_Values/10m_SA_entropy_smooth_masked_5nm_summary.csv`
- `Quad_Values/10m_SA_entropy_smooth_masked_5nm_boot_long.csv`
- `Quad_Values/10m_SA_entropy_smooth_masked_5nm_boot_wide.csv`
- `Quad_Values/20m_SA_entropy_smooth_masked_5nm_summary.csv`
- `Quad_Values/20m_SA_entropy_smooth_masked_5nm_boot_long.csv`
- `Quad_Values/20m_SA_entropy_smooth_masked_5nm_boot_wide.csv`
- `Quad_Values/50m_SA_entropy_smooth_masked_5nm_summary.csv`
- `Quad_Values/50m_SA_entropy_smooth_masked_5nm_boot_long.csv`
- `Quad_Values/50m_SA_entropy_smooth_masked_5nm_boot_wide.csv`
- `Quad_Values/Spectral_diversitySHPs/SA_entropy_10m_smooth_masked_5nm.shp`
- `Quad_Values/Spectral_diversitySHPs/SA_entropy_20m_smooth_masked_5nm.shp`
- `Quad_Values/Spectral_diversitySHPs/SA_entropy_50m_smooth_masked_5nm.shp`

## Notes

- The all-pixel path was explored automatically and used only where possible under the configured pair-count threshold.
- Most quadrats used `bootstrap_mean`, which is expected because representative masked quadrats had millions to billions of possible all-pixel spectral angle pairs.
- Missing 10 m shapefile values include 91 shapefile quadrats without matching current raster outputs and 6 raster outputs with insufficient valid sunlit pixels.
- Missing 20 m shapefile values are 15 shapefile quadrats without matching current raster outputs.

## Result

The spectral heterogeneity workflow completed successfully for 10 m, 20 m, and 50 m current smoothed 5 nm spectra.
