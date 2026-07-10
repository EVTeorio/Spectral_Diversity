# ExecPlan: 50 m PC1 Mean Reflectance Correlation

## Objective

Quantify the relationship between mean spectral reflectance and PC1 using the existing global PCA basis and a 20,000-pixel sample from current 50 m quadrat spectra.

## Requested Task

- Use the existing PCA object in `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`.
- Sample 250 valid sunlit pixels from each 50 m quadrat, for 20,000 pixels total.
- Compare sampled-pixel mean reflectance to PC1.
- Produce Pearson `r` and `R2` values.

## Files To Review

- `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R`
- `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`
- `Quad_Spectra/50m_smooth_5nm/`

## Proposed Changes

- Add a focused analysis script under `scripts/3_Analysis/`.
- Write output tables under `reports/tables/pc1_mean_reflectance_correlation/`.
- Write figures under `reports/figures/pc1_mean_reflectance_correlation/`.
- Write a short analysis report under `reports/analysis/`.
- Add task and validation notes.

## Validation Plan

- Confirm 80 current 50 m raster files are found.
- Confirm exactly 250 pixels are sampled per 50 m quadrat.
- Confirm the pixel-level sample contains 20,000 rows.
- Confirm output tables, figures, and report are created.
- Confirm the script parses under R 4.2.3.

## Risks

- PC axis sign is arbitrary, so the sign of `r` should be interpreted as orientation of this fitted PCA object, not as an absolute ecological direction.
- Manual atmospheric/cloud-excluded 50 m quadrats are included because the requested design targets all 80 50 m quadrats; the output flags these quadrats for sensitivity review.
