# Task Report: 10 m Footprint PCA and Standardized PCA Rebuild

Date: 2026-07-07

## Request

Replace the previous PCA basis that sampled nested 10 m, 20 m, and 50 m quadrat rasters with a PCA basis sampled from the full available 10 m raster footprint only. Build a second PCA using the same sampling logic after vector-normalizing spectra, name those outputs with `standardized_PCA`, and regenerate downstream PCA-dependent spectral heterogeneity outputs.

## Work Completed

- Updated `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R` so PCA sampling uses only `Quad_Spectra/10m_smooth_5nm`.
- Confirmed the 563 nm shadow mask is applied before PCA sampling and before vector normalization.
- Requested 450 illuminated pixels per 10 m raster.
- Rebuilt the raw PCA object from 854,767 sampled spectra across 1,909 current 10 m rasters.
- Built a separate vector-normalized standardized PCA object from the same sampled spectra.
- Regenerated raw PCA metric CSVs, standardized PCA metric CSVs, combined spectral heterogeneity CSVs, and combined spectral heterogeneity shapefile outputs.
- Added standardized PCA metric fields to the root-level combined quadrat analysis tables.
- Regenerated `reports/combined_quadrat_variable_guide.md` and `combined_quadrat_variable_guide.docx`.
- Reran the 50 m PC mean-reflectance diagnostic and the PC1/PC2 loading spectral-region diagnostic against both the rebuilt raw PCA object and the vector-normalized standardized PCA object.
- Updated project documentation to distinguish PCA-basis sampling from downstream manual metric exclusions.

## Key Outputs

- `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`
- `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds`
- `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm_variance_explained.csv`
- `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm_variance_explained.csv`
- `Quad_Values/*_PCA_metrics_smooth_masked_5nm_summary.csv`
- `Quad_Values/*_standardized_PCA_metrics_smooth_masked_5nm_summary.csv`
- `Quad_Values/*_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
- `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_*_smooth_masked_5nm.*`
- `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`

## Notes

The PCA basis now samples the available 10 m footprint only. The downstream PCA metric tables still retain the existing manual atmospheric/cloud exclusion policy for metric reporting, so excluded quadrats have missing raw and standardized PCA-dependent metric values while remaining part of the 10 m footprint PCA basis when their source raster is available.

The multiscale biodiversity script was updated to include standardized PCA sensitivity responses, but biodiversity/environment relationship analyses were intentionally not rerun as part of the current PCA-only diagnostic update.
