# ExecPlan: 50 m PC1-PC3 Mean Reflectance Correlation

## Objective

Extend the existing 50 m mean-reflectance diagnostic from PC1 to PC1, PC2, and PC3 using the same deterministic 20,000-pixel sample design.

## Requested Task

- Use the existing PCA object in `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`.
- Sample 250 valid sunlit pixels from each 50 m quadrat, for 20,000 pixels total.
- Compare sampled-pixel mean reflectance to PC1, PC2, and PC3.
- Produce Pearson `r`, `R2`, p-values, and confidence intervals.

## Implementation

- Generalize `scripts/3_Analysis/pc1_mean_reflectance_correlation_50m.R` so the same sampled pixels are projected onto PC1-PC3.
- Write combined PC1-PC3 output tables under `reports/tables/pc1_mean_reflectance_correlation/`.
- Write pixel-level and quadrat-mean scatterplots for each PC axis under `reports/figures/pc1_mean_reflectance_correlation/`.
- Write a short analysis report under `reports/analysis/`.

## Validation Plan

- Confirm 80 current 50 m raster files are found.
- Confirm exactly 250 pixels are sampled per 50 m quadrat.
- Confirm the pixel-level sample contains 20,000 rows.
- Confirm the correlation summary contains six rows: pixel-level and quadrat-mean for PC1, PC2, and PC3.
- Confirm the script parses under R 4.2.3.

## Risks

- PC axis signs are arbitrary, so correlation direction should be interpreted relative to the saved PCA object's orientation.
- Pixel-level p-values can become very small because the sample has 20,000 pixels; effect sizes and quadrat-level checks should carry most of the interpretation.
