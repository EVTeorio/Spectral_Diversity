# Task Report: Spectral Heterogeneity Workflow Update

Date: 2026-06-18

## Objective

Calculate per-quadrat spectral heterogeneity values for each resolution from the current smoothed and 5 nm resampled spectra.

## Inputs

- `Quad_Spectra/10m_smooth_5nm`
- `Quad_Spectra/20m_smooth_5nm`
- `Quad_Spectra/50m_smooth_5nm`

## Script Updated

- `scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R`
- `scripts/2_Indices Creation/Spectral_diversity/run_sa_entropy_scale.R`

## Key Changes

- Updated the workflow to process 10 m, 20 m, and 50 m spectra instead of only the older 20 m folder.
- Replaced the older `20m_smoothed_5nm` input with the current `*_smooth_5nm` inputs.
- Preserved the shadow masking parameters:
  - band: `563 nm`
  - threshold: `0.0305476`
  - current retained sunlit condition: reflectance greater than the threshold when `DIRECTION <- "<"`.
- Added scale-specific shapefile joins:
  - 10 m joins raster IDs such as `0_a` to `PR_10m.shp` field `sub_id`
  - 20 m joins numeric raster IDs to `PR_20m.shp` field `Name`
  - 50 m joins IDs such as `sub50_1` to `PR_50m.shp` field `Name`
- Added a guarded exact all-pixel entropy path when pairwise angle counts are below `MAX_EXACT_PAIRS`.
- Added automatic fallback to bootstrap mean entropy when all-pixel pair counts are too large.
- Reworked pairwise angle calculation to use preallocated vectors and block matrix multiplication rather than repeatedly growing angle vectors.
- Added sampled-pair entropy within large bootstrap replicates to avoid materializing all pairwise angles among 5,000 sampled pixels.
- Removed unnecessary tidyverse and `sf` dependencies from the workflow; full raster execution now requires `terra`, with optional completion sound from `beepr`.

## Run Settings

- Bootstrap replicates: 70
- Pixel subsample per bootstrap: 5,000
- Pair sample per large bootstrap replicate: 10,000
- Exact all-pixel entropy threshold: 2,000,000 possible pairs
- Exact bootstrap-pair threshold: 250,000 possible pairs

## Outputs Created

- `Quad_Values/10m_SA_entropy_smooth_masked_5nm_summary.csv`
- `Quad_Values/20m_SA_entropy_smooth_masked_5nm_summary.csv`
- `Quad_Values/50m_SA_entropy_smooth_masked_5nm_summary.csv`
- matching long and wide bootstrap CSVs for each scale
- `Quad_Values/Spectral_diversitySHPs/SA_entropy_10m_smooth_masked_5nm.shp`
- `Quad_Values/Spectral_diversitySHPs/SA_entropy_20m_smooth_masked_5nm.shp`
- `Quad_Values/Spectral_diversitySHPs/SA_entropy_50m_smooth_masked_5nm.shp`

## Results

| Scale | Rasters Processed | Summary Rows | Primary Method Counts |
|---|---:|---:|---|
| 10 m | 1,909 | 1,909 | 1,866 bootstrap mean; 37 exact all-pixel; 6 insufficient pixels |
| 20 m | 485 | 485 | 477 bootstrap mean; 8 exact all-pixel |
| 50 m | 80 | 80 | 79 bootstrap mean; 1 exact all-pixel |

## Notes

- Do not describe the current SA entropy values as if all pairwise spectral angle combinations were calculated for every quadrat. The `n_pairs` field records the possible number of all-retained-pixel pairs, but `method = "bootstrap_mean"` rows report the mean of 70 bootstrap iterations, not exhaustive all-pair entropy.
- Exact all-pixel entropy remains quadratic in pixel count. The script explores that path only for small enough masked quadrats and otherwise uses `boot_mean` as the primary spectral heterogeneity value.
- Representative masked quadrats had millions to billions of possible all-pixel spectral angle pairs, so fallback to `boot_mean` was expected for most quadrats.
- The earlier package visibility issue was sandbox-specific. `terra` and `beepr` were found in `C:/Users/PaintRock/AppData/Local/R/win-library/4.2` when Codex was allowed to access that R user library path.
