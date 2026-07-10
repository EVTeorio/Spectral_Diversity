# 50 m PC1-PC3 and Mean Reflectance Correlation

Date: 2026-07-07

## Purpose

Evaluate how strongly PC1, PC2, and PC3 are associated with simple mean spectral reflectance for both the regular PCA and vector-normalized standardized PCA. The key interpretation is based on standardized PCA because it is the brightness-reduced PCA basis.

## Design

- Sampled pixels: 250 valid sunlit pixels per 50 m quadrat.
- Total sample per PCA basis: 20,000 pixels from 80 quadrats.
- Spectra source: `Quad_Spectra/50m_smooth_5nm/`.
- PCA sources: regular and standardized PCA RDS files in `Quad_Values/Spectral_diversitySHPs/`.
- Shadow mask: retain pixels with `563 nm` reflectance greater than 0.0305476.
- Manual PCA-excluded 50 m quadrats were included in this diagnostic sample and are flagged in the outputs.

## Key Results: Standardized PCA

| pca_basis | pca_label | analysis_level | pc_axis | n_rows | n_quadrats | pearson_r | r_squared | p_value | conf_low | conf_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| standardized_PCA | Standardized PCA | pixel_level | PC1 | 20000.000000 | 80.000000 | -0.117768 | 0.013869 | 0.000000 | -0.131412 | -0.104078 |
| standardized_PCA | Standardized PCA | quadrat_mean_level | PC1 | 80.000000 | 80.000000 | -0.348246 | 0.121275 | 0.001548 | -0.527594 | -0.139178 |
| standardized_PCA | Standardized PCA | pixel_level | PC2 | 20000.000000 | 80.000000 | -0.331173 | 0.109676 | 0.000000 | -0.343456 | -0.318777 |
| standardized_PCA | Standardized PCA | quadrat_mean_level | PC2 | 80.000000 | 80.000000 | 0.059775 | 0.003573 | 0.598404 | -0.162071 | 0.275868 |
| standardized_PCA | Standardized PCA | pixel_level | PC3 | 20000.000000 | 80.000000 | 0.110240 | 0.012153 | 0.000000 | 0.096528 | 0.123910 |
| standardized_PCA | Standardized PCA | quadrat_mean_level | PC3 | 80.000000 | 80.000000 | -0.265605 | 0.070546 | 0.017257 | -0.458561 | -0.048732 |

## Supporting Regular PCA Results

| pca_basis | pca_label | analysis_level | pc_axis | n_rows | n_quadrats | pearson_r | r_squared | p_value | conf_low | conf_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| regular_PCA | Regular PCA | pixel_level | PC1 | 20000.000000 | 80.000000 | 0.893074 | 0.797582 | 0.000000 | 0.890234 | 0.895845 |
| regular_PCA | Regular PCA | quadrat_mean_level | PC1 | 80.000000 | 80.000000 | 0.925360 | 0.856291 | 0.000000 | 0.885726 | 0.951600 |
| regular_PCA | Regular PCA | pixel_level | PC2 | 20000.000000 | 80.000000 | -0.426810 | 0.182166 | 0.000000 | -0.438077 | -0.415408 |
| regular_PCA | Regular PCA | quadrat_mean_level | PC2 | 80.000000 | 80.000000 | -0.627785 | 0.394114 | 0.000000 | -0.744772 | -0.473362 |
| regular_PCA | Regular PCA | pixel_level | PC3 | 20000.000000 | 80.000000 | 0.033676 | 0.001134 | 0.000002 | 0.019826 | 0.047513 |
| regular_PCA | Regular PCA | quadrat_mean_level | PC3 | 80.000000 | 80.000000 | 0.565265 | 0.319525 | 0.000000 | 0.394551 | 0.698259 |

## Interpretation Notes

- The pixel-level result is the requested 20,000-pixel diagnostic per PCA basis.
- The quadrat-mean result is a sensitivity check where each quadrat contributes one mean point.
- PC axis signs are arbitrary, so the sign of the correlation reflects the orientation of each saved PCA object.

## Output Tables

- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_pixel_sample.csv`
- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_quadrat_means.csv`
- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_correlation_summary.csv`
- `reports/tables/pc1_mean_reflectance_correlation/50m_regular_PCA_pc_mean_reflectance_correlation_summary.csv`
- `reports/tables/pc1_mean_reflectance_correlation/50m_standardized_PCA_pc_mean_reflectance_correlation_summary.csv`
