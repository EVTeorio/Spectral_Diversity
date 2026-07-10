# Task Report: Standardized PCA Diagnostics and Sample-Size Completion

Date: 2026-07-09

## Objective

Resume the interrupted PCA-only work after OneDrive access was restored. Complete regular and standardized PCA diagnostic outputs without rerunning biodiversity or environmental relationship analyses.

## Completed

- Patched `scripts/3_Analysis/spectral_metric_sample_size_effects.R` so standardized alpha-hull sample-size runs correctly request the alpha-hull calculation after metric IDs receive the `standardized_PCA_` prefix.
- Reran the 50 m mean-reflectance correlation diagnostic for both regular PCA and standardized PCA.
- Reran PCA loading spectral-region diagnostics for both regular PCA and standardized PCA.
- Reran PCA-derived sample-size analysis for:
  - regular PCA mean distance
  - regular PCA spectral Rao's Q
  - regular PCA alpha-hull area
  - standardized PCA mean distance
  - standardized PCA spectral Rao's Q
  - standardized PCA alpha-hull area
- Regenerated the regular PCA and standardized PCA sample-size Markdown reports in report-only mode after all metric folders were current.

## Key Standardized PCA Findings

- Standardized PCA loading analysis identifies the main brightness-reduced PC2 region at 753-778 nm.
- Standardized PC2 has broader nonuniform structure at 753-813 nm.
- Standardized PC1 has its strongest nonuniform loading region at 843-903 nm.
- In the 50 m mean-reflectance diagnostic, standardized PC1 is only weakly associated with pixel-level mean reflectance (r = -0.1178, R2 = 0.0139), supporting the use of standardized PCA as the brightness-reduced interpretation.

## Outputs

- `reports/analysis/20260707_50m_pc_mean_reflectance_correlation.md`
- `reports/analysis/20260707_pca_loading_spectral_regions.md`
- `reports/analysis/20260704_pca_metric_sample_size_effects.md`
- `reports/analysis/20260707_standardized_pca_metric_sample_size_effects.md`
- `reports/tables/pc1_mean_reflectance_correlation/`
- `reports/tables/pca_loading_spectral_regions/`
- `reports/tables/sample_size_effects/standardized_PCA_pca_mean_distance/`
- `reports/tables/sample_size_effects/standardized_PCA_spectral_rao_q/`
- `reports/tables/sample_size_effects/standardized_PCA_alpha_hull_area/`

