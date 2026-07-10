# ExecPlan: PCA Loading Spectral Region Analysis

## Objective

Identify wavelength regions where PC1 and PC2 loadings are not uniform, then isolate PC2 loading regions that are strong without being dominated by broad PC1 brightness structure.

## Requested Structure

- Determine where PC1 loadings depart from uniform brightness.
- Determine where PC2 loadings are nonuniform.
- Compare PC2 loadings against PC1 departures to find wavelength regions where PC2 has strong structure not explained by broad brightness.
- Record the finding in project state as a key result.

## Implementation

- Use `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`.
- Extract PC1 and PC2 rotation/loadings by wavelength.
- Compare PC1 to a flat loading vector of `1 / sqrt(n_bands)`.
- Flag PC1 and PC2 nonuniform bands using the top 10% of absolute z-scored departure/loading values.
- Flag non-brightness PC2 bands where PC2 is strong and PC1 departure from uniform brightness is less than 1 z-score.
- Collapse adjacent 5 nm bands into interpretable spectral regions.

## Outputs

- `scripts/3_Analysis/pca_loading_spectral_region_analysis.R`
- `reports/analysis/20260707_pca_loading_spectral_regions.md`
- `reports/tables/pca_loading_spectral_regions/`
- `reports/figures/pca_loading_spectral_regions/`
- `reports/tasks/20260707_pca_loading_spectral_regions.md`
- `reports/validation/20260707_pca_loading_spectral_regions_validation.md`

## Validation Plan

- Confirm the PCA rotation matrix is read from the saved PCA object.
- Confirm all 121 wavelength bands are represented in the per-band loading table.
- Confirm PC1, PC2, and non-brightness PC2 regions are written.
- Confirm the script parses under R 4.2.3.
