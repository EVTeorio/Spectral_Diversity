# Validation Report: Quadrat Crop Reproducibility Attempt

Date: 2026-06-15

## Validation Checks

| Check | Result | Notes |
|---|---|---|
| R 4.2.3 available | Passed | `C:/Program Files/R/R-4.2.3/bin/Rscript.exe` was used. |
| Required R package library path updated | Passed | Test configuration now uses `C:/Users/PaintRock/AppData/Local/R/win-library/4.2` and `C:/Program Files/R/R-4.2.3/library`. |
| Temporary project-local package library removed | Passed | `r_libs/` was removed from the working directory. |
| Required R packages load from user library | Failed | `terra` and `sf` load, but `dplyr` and `testthat` fail because `rlang` is not properly installed in the user library. |
| Test file created | Passed | `tests/test_quad_crop_reproducibility.R` exists. |
| Source raster readability | Failed | Most `HSI_NA_trimmed/` rasters could not be opened by `terra`, likely due to OneDrive cloud-only status. |
| Crop reproducibility test completed | Failed/interrupted | The final test run was manually interrupted before completion. |
| Crop script rerun | Not performed | No scale was regenerated because validation could not reliably determine alignment. |

## Readable Source Rasters During Attempt

- `HSI_NA_trimmed/raw_0_rd_rf_or`
- `HSI_NA_trimmed/raw_11159_rd_rf_or`

## Blocker

The crop reproducibility test requires locally readable source HSI rasters. Until the relevant `HSI_NA_trimmed/` files are hydrated from OneDrive, the comparison cannot reliably confirm whether all existing 10 m, 20 m, and 50 m quadrat spectra match the crop script output.

A second environment blocker is present in the user R package library: `dplyr` and `testthat` fail to load because the installed `rlang` DLL version does not match the package version metadata. Active R/RStudio sessions may be locking the DLL.

## Recommended Next Step

Hydrate the `HSI_NA_trimmed/` directory locally, then rerun:

```powershell
& 'C:\Program Files\R\R-4.2.3\bin\Rscript.exe' -e ".libPaths(c('C:/Users/PaintRock/AppData/Local/R/win-library/4.2','C:/Program Files/R/R-4.2.3/library')); testthat::test_file('tests/test_quad_crop_reproducibility.R', reporter='summary')"
```
