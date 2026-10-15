# Data Dictionary

Last updated: 2026-10-14

This is the concise inventory of current and archived project data. It identifies one canonical location for each active product and marks known redundant or historical variants. It does not modify protected manuscript materials in `Documents/`.

## Active Inputs and Planned Scales

| Dataset or location | Status | Purpose | Key or identifier | Notes |
|---|---|---|---|---|
| `HSI_NA_trimmed/` | current | Trimmed hyperspectral source rasters | raster tile name | Raw imagery source for quadrat extraction. |
| `Quad_Scale_SHPs/PR_10m.*`, `PR_20m.*`, `PR_50m.*` | current | Established quadrat boundaries | scale-specific quadrat ID | Current established scales. |
| `Quad_Scale_SHPs/PR_25m.*`, `PR_33m.*` | planned | Expansion quadrat boundaries | scale-specific quadrat ID | Do not remove; 25 m and 33 m are planned scales. |
| `Quad_Spectra/10m/`, `20m/`, `50m/` | current | Confirmed base quadrat spectra | raster filename / quadrat ID | Use as the reproducible base spectral inputs. |
| `Quad_Spectra/10m_smooth_5nm/` | pending restore | Current intended 10 m resampled spectra | raster filename / quadrat ID | Accidentally deleted; user will restore. |
| `Quad_Spectra/20m_smooth_5nm/`, `50m_smooth_5nm/` | current | Current resampled spectra | raster filename / quadrat ID | Inputs for existing 20 m and 50 m spectral outputs. |
| `Quad_Spectra/25m/`, `25m_smooth_5nm/`, `33m/` | planned or partial | 25 m and 33 m scale expansion | raster filename / quadrat ID | Empty or partial state is intentional. |

Intermediate `*_smooth/` folders were deliberately removed to save space. They are reproducible products, not missing canonical source data.

## Current Derived Biodiversity and Environmental Data

| Canonical dataset | Location | Rows currently expected | Primary key | Content |
|---|---|---:|---|---|
| Plant diversity, 10 m | `Quad_Values/Diversity_SHPs/plant_diversity_10m.csv` | 2,000 | `scale`, `quad_id` | Species composition plus richness, Shannon, Simpson, evenness, Faith PD, phylogenetic Rao's Q, and abundance-weighted Faith PD. |
| Plant diversity, 20 m | `Quad_Values/Diversity_SHPs/plant_diversity_20m.csv` | 500 | `scale`, `quad_id` | Same diversity fields at 20 m. |
| Plant diversity, 50 m | `Quad_Values/Diversity_SHPs/plant_diversity_50m.csv` | 80 | `scale`, `quad_id` | Same diversity fields at 50 m. |
| Environmental covariates | `Quad_Values/Enviro_SHPs/enviro_variables_*m.csv` | 2,000 / 500 / 80 | `scale`, `quad_id` | Mean elevation and 5x5, 11x11, and 21x21 topographic roughness summaries. |

The matching shapefiles are spatial delivery copies of these values. Preserve each shapefile component set as a unit (`.shp`, `.shx`, `.dbf`, `.prj`, and `.cpg` where present).

## Current Priority Spectral Data

| Canonical dataset | Location | Primary key | Use |
|---|---|---|---|
| Original PCA basis | `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds` | PCA model object | Preserve as a current comparison basis. |
| Vector-normalized PCA basis | `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds` | PCA model object | Preserve as a current brightness-reduced comparison basis. |
| Original PCA variance description | `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm_variance_explained.csv` | `pc_axis` | Original-PCA provenance and diagnostics. |
| PCA variance description | `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm_variance_explained.csv` | `pc_axis` | PCA provenance and diagnostics. |
| Combined spectral summary, 10 m | `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_10m_smooth_masked_5nm_summary.csv` | `quad_id` | Canonical 10 m spectral source table. |
| Combined spectral summary, 20 m | `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_20m_smooth_masked_5nm_summary.csv` | `quad_id` | Canonical 20 m spectral source table. |
| Combined spectral summary, 50 m | `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_50m_smooth_masked_5nm_summary.csv` | `quad_id` | Canonical 50 m spectral source table. |
| All-pixel alpha-hull values | `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values.csv` | `scale`, `quad_id`, `basis`, `alpha_value` | Appendable 50 m production-style alpha-hull values. |
| Sampled 2D PCA distance and Rao's Q | `reports/tables/sample_2d_pca_mean_rao_values/sample_2d_pca_mean_rao_values.csv` | `scale`, `quad_id`, `basis`, `sample_max_points` | Appendable sampled 50 m comparison values. |

The current priority fields in the canonical combined spectral summaries and joined analysis tables are `pca_euclidean_mean` / `spec_pca_mean`, `alpha_hull_area` / `spec_alpha`, `standardized_PCA_pca_euclidean_mean` / `spec_spca_mean`, and `standardized_PCA_alpha_hull_area` / `spec_spca_alpha`.

## Analysis-Ready Tables

| Dataset | Location | Primary key | Status |
|---|---|---|---|
| 10 m joined table | `quadrat_analysis_10m.csv` | `scale`, `quad_id` | Current established-scale analysis table. |
| 20 m joined table | `quadrat_analysis_20m.csv` | `scale`, `quad_id` | Current established-scale analysis table. |
| 50 m joined table | `quadrat_analysis_50m.csv` | `scale`, `quad_id` | Current established-scale analysis table. |

These tables join current diversity, spectral, environmental, and coordinate fields. They retain broader historical spectral fields for reproducibility; forward analysis should select the current priority fields rather than treat every column as active.

## Archived and Historical Data

| Location or product | Classification | Retention rule |
|---|---|---|
| `Quad_Values/Spectral_diversitySHPs_first/` | removed first-pass spectral archive | Removed during archive cleanup; no longer present in the live directory tree. Current canonical spectral products remain in `Quad_Values/Spectral_diversitySHPs/`. |
| `Quad_Values/*_PCA_metrics_*`, `*_standardized_PCA_metrics_*`, and `*_SA_entropy_*` | component/provenance tables | Preserve as reproducible source components; they are not active primary analysis tables. |
| `Quad_Values/20m_SA_smooth_masked*.csv`, `20m_SA_entrop_boot*.csv`, `20m_spectral_sp.csv`, and `Other_variables/` | earlier 20 m-only products | Preserve as compact historical data; do not use as current integrated analysis inputs. |
| Bootstrap, sample-size, hull benchmark, and relationship tables under `reports/tables/` | diagnostics and sensitivity provenance | Retain unless byte-identical duplicates are verified and unreferenced. |
| Convex-hull, 3D hull, and spectral-angle-entropy values | non-primary metric archive | Keep available for later comparison, but exclude from default forward analysis. |

## Verified Redundant Copies

The following root-level files are byte-identical copies of the canonical combined spectral summaries under `Quad_Values/Spectral_diversitySHPs/`:

- `Quad_Values/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
- `Quad_Values/20m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`
- `Quad_Values/50m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`

The copies are retained pending an explicit data-deletion decision. New scripts and documentation should use the `Spectral_diversitySHPs/` copies.

## Data Handling Notes

- Treat an unnamed first CSV column as an index column unless script review proves otherwise.
- Join quadrat-level products with `scale` plus `quad_id`.
- Do not delete individual shapefile sidecars; each is part of a single spatial dataset.
- Do not classify incomplete 25 m or 33 m resources as errors while those scales are under construction.
