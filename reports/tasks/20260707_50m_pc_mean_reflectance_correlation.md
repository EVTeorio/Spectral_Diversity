# Task Report: 50 m PC1-PC3 Mean Reflectance Correlation

Date: 2026-07-07

## Objective

Compare mean spectral reflectance to PC1, PC2, and PC3 for regular and standardized PCA using 250 sampled sunlit pixels from each current 50 m quadrat.

## Key Output

The combined summary table is `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_correlation_summary.csv`; basis-specific summaries are also written for regular and standardized PCA.

## Standardized PCA Results

| pca_basis | pca_label | analysis_level | pc_axis | n_rows | n_quadrats | pearson_r | r_squared | p_value | conf_low | conf_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| standardized_PCA | Standardized PCA | pixel_level | PC1 | 20000.000000 | 80.000000 | -0.117768 | 0.013869 | 0.000000 | -0.131412 | -0.104078 |
| standardized_PCA | Standardized PCA | quadrat_mean_level | PC1 | 80.000000 | 80.000000 | -0.348246 | 0.121275 | 0.001548 | -0.527594 | -0.139178 |
| standardized_PCA | Standardized PCA | pixel_level | PC2 | 20000.000000 | 80.000000 | -0.331173 | 0.109676 | 0.000000 | -0.343456 | -0.318777 |
| standardized_PCA | Standardized PCA | quadrat_mean_level | PC2 | 80.000000 | 80.000000 | 0.059775 | 0.003573 | 0.598404 | -0.162071 | 0.275868 |
| standardized_PCA | Standardized PCA | pixel_level | PC3 | 20000.000000 | 80.000000 | 0.110240 | 0.012153 | 0.000000 | 0.096528 | 0.123910 |
| standardized_PCA | Standardized PCA | quadrat_mean_level | PC3 | 80.000000 | 80.000000 | -0.265605 | 0.070546 | 0.017257 | -0.458561 | -0.048732 |
