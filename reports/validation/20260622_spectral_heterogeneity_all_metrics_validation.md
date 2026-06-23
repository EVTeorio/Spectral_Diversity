# Validation Report: Spectral Heterogeneity All-Metrics Workflow

Date: 2026-06-22

Updated: 2026-06-22, exclusion-aware PCA rebuild

## Scope

Validate the all-scale spectral heterogeneity outputs that combine the existing spectral angle entropy values with three additional PCA-dependent measures:

- mean Euclidean distance from the quadrat centroid in global PC1-PC3 space
- alpha-hull area in global PC1-PC2 space
- Rao's Q using equal pixel weights and squared Euclidean distance in global PC1-PC3 space

A supplemental 3-axis hull diagnostic was also generated as convex hull volume and surface area in global PC1-PC3 space. A true 3D alpha shape was not used because `alphashape3d` is not installed; the implemented 3D value is therefore documented as a PCA convex hull volume, not as a 3D alpha hull.

## Supersession Note

The PCA-dependent outputs created earlier on 2026-06-22 before applying the manual atmospheric/cloud exclusion set should be disregarded. The current PCA object, PCA diagnostics, PCA metric CSVs, and combined shapefiles were regenerated after excluding the documented affected quadrats from the global PCA sample and from PCA-dependent metric calculation.

## Inputs

- `Quad_Spectra/10m_smooth_5nm`
- `Quad_Spectra/20m_smooth_5nm`
- `Quad_Spectra/50m_smooth_5nm`
- `Indices_SHPs/*_SA_entropy_smooth_masked_5nm_summary.csv`
- `Quad_Scale_SHPs/PR_10m.shp`
- `Quad_Scale_SHPs/PR_20m.shp`
- `Quad_Scale_SHPs/PR_50m.shp`

The workflow used the validated spectral-angle entropy shadow mask: `563 nm`, threshold `0.0305476`, retaining pixels greater than the threshold.

## Manual Exclusions

The exclusion policy is encoded in `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R` as `manual_atmospheric_cloud_exclusions_20260622`.

Excluded current raster outputs:

| Scale | Excluded rasters | Rule |
|---|---:|---|
| 10 m | 165 | Current 10 m rasters whose parent 20 m quadrat is in the documented atmospheric/cloud list |
| 20 m | 49 | Current 20 m rasters in the documented atmospheric/cloud list |
| 50 m | 6 | `sub50_80`, `sub50_79`, `sub50_71`, `sub50_70`, `sub50_62`, `sub50_53` |

Some documented 20 m IDs do not have current smoothed 5 nm raster outputs, so the executable current-raster exclusion count is smaller than the full documented ID list.

## PCA Basis

The global PCA was rebuilt after removing the manual exclusions. It used a uniform per-eligible-quadrat sample across all scales. The requested target was up to 500 retained sunlit pixels per eligible quadrat raster; the configured maximum PCA sample of 800,000 rows produced an effective uniform sample of 354 pixels per eligible raster.

PCA sample summary:

| Field | Value |
|---|---:|
| Eligible rasters | 2,254 |
| Sample rows | 797,916 |
| Bands | 121 |
| Effective sample per eligible raster | 354 |

Variance explained by the first axes:

| Axis | Percent variance | Cumulative percent |
|---:|---:|---:|
| PC1 | 66.98% | 66.98% |
| PC2 | 19.93% | 86.91% |
| PC3 | 3.80% | 90.71% |

PCA diagnostic outputs:

- `Indices_SHPs/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`
- `Indices_SHPs/Spectral_diversitySHPs/global_pca_smooth_masked_5nm_variance_explained.csv`
- `reports/figures/spectral_heterogeneity/global_pca_smooth_masked_5nm_variance_explained.png`

## Output Validation

| Scale | CSV rows | Shapefile features | Manual exclusions | Missing SA entropy | Missing PCA Euclidean | Missing alpha hull | Missing 3D hull volume | Missing Rao Q |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 m | 1,909 | 2,000 | 165 | 6 | 165 | 168 | 165 | 165 |
| 20 m | 485 | 500 | 49 | 0 | 49 | 49 | 49 | 49 |
| 50 m | 80 | 80 | 6 | 0 | 6 | 6 | 6 | 6 |

Shapefile missing values reflect unmatched quadrat polygons plus metric-level missingness:

| Scale | Missing `sa_ent` | Missing `pca_eucl` | Missing `alpha_h` | Missing `hull3d_v` | Missing `rao_q` | `pca_excl` true |
|---|---:|---:|---:|---:|---:|---:|
| 10 m | 97 | 256 | 259 | 256 | 256 | 165 |
| 20 m | 15 | 64 | 64 | 64 | 64 | 49 |
| 50 m | 0 | 6 | 6 | 6 | 6 | 6 |

## Method Counts

PCA metric method:

| Scale | All pixels | Manual excluded |
|---|---:|---:|
| 10 m | 1,744 | 165 |
| 20 m | 436 | 49 |
| 50 m | 74 | 6 |

Alpha hull:

| Scale | All pixels | Sampled pixels | Fallback sampled pixels | Failed fallback | Manual excluded |
|---|---:|---:|---:|---:|---:|
| 10 m | 1,731 | 0 | 10 | 3 | 165 |
| 20 m | 2 | 434 | 0 | 0 | 49 |
| 50 m | 0 | 74 | 0 | 0 | 6 |

## Output Files

Per-scale combined CSVs were written both beside the earlier SA entropy CSVs and inside the requested shapefile output directory:

- `Indices_SHPs/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/20m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/50m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_10m_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_20m_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_50m_smooth_masked_5nm_summary.csv`

Per-scale combined shapefiles:

- `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_10m_smooth_masked_5nm.shp`
- `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_20m_smooth_masked_5nm.shp`
- `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_50m_smooth_masked_5nm.shp`

PCA-only metric CSVs:

- `Indices_SHPs/10m_PCA_metrics_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/20m_PCA_metrics_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/50m_PCA_metrics_smooth_masked_5nm_summary.csv`

## Result

Pass. The current outputs contain four primary spectral heterogeneity measures: SA entropy, PCA Euclidean heterogeneity, alpha-hull area, and Rao's Q. The PCA-dependent values were regenerated with manual atmospheric/cloud exclusions applied; excluded quadrats retain SA entropy values where available but have missing PCA-dependent metrics.
