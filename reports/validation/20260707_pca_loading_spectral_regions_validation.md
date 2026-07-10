# Validation: PCA Loading Spectral Region Analysis

Date: 2026-07-07

## Checks

| pca_basis | pca_label | metric | value |
| --- | --- | --- | --- |
| regular_PCA | Regular PCA | n_bands | 121.0000 |
| regular_PCA | Regular PCA | wavelength_min_nm | 398.0000 |
| regular_PCA | Regular PCA | wavelength_max_nm | 998.0000 |
| regular_PCA | Regular PCA | pc1_uniform_loading | 0.0909 |
| regular_PCA | Regular PCA | pc1_nonuniform_abs_z_threshold | 1.4736 |
| regular_PCA | Regular PCA | pc2_nonuniform_abs_z_threshold | 1.1654 |
| regular_PCA | Regular PCA | brightness_flat_abs_pc1_departure_z_threshold | 1.0000 |
| regular_PCA | Regular PCA | pc1_nonuniform_band_count | 13.0000 |
| regular_PCA | Regular PCA | pc2_nonuniform_band_count | 13.0000 |
| regular_PCA | Regular PCA | pc2_strong_not_brightness_band_count | 12.0000 |
| standardized_PCA | Standardized PCA | n_bands | 121.0000 |
| standardized_PCA | Standardized PCA | wavelength_min_nm | 398.0000 |
| standardized_PCA | Standardized PCA | wavelength_max_nm | 998.0000 |
| standardized_PCA | Standardized PCA | pc1_uniform_loading | 0.0909 |
| standardized_PCA | Standardized PCA | pc1_nonuniform_abs_z_threshold | 1.3402 |
| standardized_PCA | Standardized PCA | pc2_nonuniform_abs_z_threshold | 1.6531 |
| standardized_PCA | Standardized PCA | brightness_flat_abs_pc1_departure_z_threshold | 1.0000 |
| standardized_PCA | Standardized PCA | pc1_nonuniform_band_count | 13.0000 |
| standardized_PCA | Standardized PCA | pc2_nonuniform_band_count | 13.0000 |
| standardized_PCA | Standardized PCA | pc2_strong_not_brightness_band_count | 6.0000 |

## Result

Pass. The analysis read both saved PCA rotation matrices, extracted PC1 and PC2 loadings for all wavelengths, wrote combined per-band and region-level tables, and generated basis-specific figures.
