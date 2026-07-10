# Validation Report: 10 m Footprint PCA and Standardized PCA Rebuild

Date: 2026-07-07

## Scope

Validate the replacement PCA workflow requested on 2026-07-07:

- raw PCA sampled from current 10 m rasters only
- standardized PCA using the same sample after vector normalization
- regenerated downstream spectral heterogeneity metrics
- standardized PCA fields added to final spectral heterogeneity and combined analysis outputs

## PCA Basis

Both PCA objects used the same source footprint:

| Field | Value |
|---|---:|
| Source scale | 10 m |
| Source rasters | 1,909 |
| Requested pixels per raster | 450 |
| Sample rows | 854,767 |
| Bands | 121 |
| Shadow mask | 563 nm, threshold 0.0305476, retain greater-than-threshold pixels |

The sample rows are slightly below 1,909 x 450 because some 10 m rasters have fewer than 450 retained illuminated pixels after masking.

## Variance Checks

Raw PCA:

| Axis | Percent variance | Cumulative percent |
|---:|---:|---:|
| PC1 | 66.23% | 66.23% |
| PC2 | 20.72% | 86.95% |
| PC3 | 3.76% | 90.71% |

Standardized PCA:

| Axis | Percent variance | Cumulative percent |
|---:|---:|---:|
| PC1 | 45.46% | 45.46% |
| PC2 | 22.54% | 68.00% |
| PC3 | 6.80% | 74.80% |

## CSV Validation

Combined spectral heterogeneity CSVs:

| File | Rows | Columns | Standardized PCA columns |
|---|---:|---:|---:|
| `Quad_Values/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv` | 1,909 | 47 | 17 |
| `Quad_Values/20m_spectral_heterogeneity_smooth_masked_5nm_summary.csv` | 485 | 47 | 17 |
| `Quad_Values/50m_spectral_heterogeneity_smooth_masked_5nm_summary.csv` | 80 | 47 | 17 |

Root combined analysis tables:

| File | Rows | Columns | Standardized PCA fields |
|---|---:|---:|---|
| `quadrat_analysis_10m.csv` | 2,000 | 32 | `spec_spca_mean`, `spec_spca_med`, `spec_spca_sd`, `spec_spca_rao`, `spec_spca_alpha`, `spec_spca_convex`, `spec_spca_hull3d_v`, `spec_spca_hull3d_a` |
| `quadrat_analysis_20m.csv` | 500 | 32 | same as 10 m |
| `quadrat_analysis_50m.csv` | 80 | 32 | same as 10 m |

Missing values in root combined analysis tables:

| Scale | Missing `spec_sa` | Missing `spec_pca_mean` | Missing `spec_spca_mean` |
|---|---:|---:|---:|
| 10 m | 97 | 256 | 256 |
| 20 m | 15 | 64 | 64 |
| 50 m | 0 | 6 | 6 |

Shapefile DBF header checks:

| Scale | DBF records | DBF fields | Standardized PCA DBF fields |
|---|---:|---:|---|
| 10 m | 2,000 | 23 | `spca_eucl`, `spca_alph`, `spca_h3dv`, `spca_h3da`, `spca_rao`, `spca_ahm`, `spca_h3dm` |
| 20 m | 500 | 22 | `spca_eucl`, `spca_alph`, `spca_h3dv`, `spca_h3da`, `spca_rao`, `spca_ahm`, `spca_h3dm` |
| 50 m | 80 | 22 | `spca_eucl`, `spca_alph`, `spca_h3dv`, `spca_h3da`, `spca_rao`, `spca_ahm`, `spca_h3dm` |

## Downstream Diagnostics

The 50 m PC mean-reflectance diagnostic was rerun against both the rebuilt raw PCA and the vector-normalized standardized PCA. The key brightness-reduced interpretation comes from standardized PCA:

| Axis | Pixel r | Pixel R2 | Quadrat-mean r | Quadrat-mean R2 |
|---|---:|---:|---:|---:|
| Standardized PC1 | -0.1178 | 0.0139 | -0.3482 | 0.1213 |
| Standardized PC2 | -0.3312 | 0.1097 | 0.0598 | 0.0036 |
| Standardized PC3 | 0.1102 | 0.0122 | -0.2656 | 0.0705 |

The PC loading spectral-region diagnostic was rerun against both regular and standardized PCA. The key finding is stated from standardized PCA: standardized PC1 has its strongest nonuniform loading region at 843-903 nm, standardized PC2 has concentrated structure at 753-813 nm, and the brightness-reduced standardized PC2 region to carry forward is 753-778 nm.

The PCA-derived sample-size sensitivity analysis was rerun for regular and standardized PCA mean distance, spectral Rao's Q, and alpha-hull area. Each of the six metric folders has 392 summary rows with finite metric means.

## Result

Pass for the PCA rebuild, standardized PCA creation, combined spectral heterogeneity CSV updates, root combined analysis table updates, PC mean-reflectance diagnostic, PC loading diagnostic, and PCA-derived sample-size sensitivity outputs.
