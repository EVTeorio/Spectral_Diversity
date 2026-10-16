# Data Dictionary

Last updated: 2026-10-15

This is the concise inventory of current and archived project data. It identifies one canonical location for each active product and marks known redundant or historical variants. It does not modify protected manuscript materials in `Documents/`.

## Active Inputs and Planned Scales

| Dataset or location | Status | Purpose | Key or identifier | Notes |
|---|---|---|---|---|
| `HSI_NA_trimmed/` | current | Trimmed hyperspectral source rasters | raster tile name | Raw imagery source for quadrat extraction. |
| `Quad_Scale_SHPs/PR_10m.*`, `PR_20m.*`, `PR_50m.*` | current | Established quadrat boundaries | scale-specific quadrat ID | Current established scales. |
| `Quad_Scale_SHPs/PR_25m.*` | active expansion | 320 quadrat boundaries | `Name` | Prepared rasters exist for 312 cells. |
| `Quad_Scale_SHPs/PR_33m.*`, `PR_100m.*` | deferred | Future quadrat boundaries | `Name` | Retain; no processing at these scales in the current task. |
| `Quad_Spectra/10m/`, `20m/`, `50m/` | current | Confirmed base quadrat spectra | raster filename / quadrat ID | Use as the reproducible base spectral inputs. |
| `Quad_Spectra/10m_smooth_5nm/` | restored and verified | Current 10 m resampled spectra | raster filename / quadrat ID | 1,909 rasters; all 1,744 eligible reference quadrats reproduce saved vector-normalized mean, Rao, and pixel counts at CSV precision. |
| `Quad_Spectra/20m_smooth_5nm/`, `50m_smooth_5nm/` | current | Current resampled spectra | raster filename / quadrat ID | Inputs for existing 20 m and 50 m spectral outputs. |
| `Quad_Spectra/25m/`, `25m_smooth_5nm/` | prepared | Active 25 m scale expansion | raster filename / quadrat ID | 312 rasters; eight grid cells have no prepared raster. Exclusion decisions belong to the user. |
| `Quad_Spectra/33m/` | deferred | Future 33 m inputs | raster filename / quadrat ID | Empty state is intentional. |

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
| Distance/Rao summary, 25 m | `Quad_Values/Spectral_diversitySHPs/spectral_mean_rao_25m_smooth_masked_5nm_summary.csv` | `scale`, `quad_id` | Canonical distance/Rao-only table in both PCA bases: 320 grid rows, 312 finite values per metric/basis, eight `missing_raster` rows. No exclusions applied or hull values calculated. |
| All-pixel alpha-hull values | `reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values.csv` | `scale`, `quad_id`, `basis`, `alpha_value` | Appendable 50 m production-style alpha-hull values. |
| Sampled 2D PCA distance and Rao's Q | `reports/tables/sample_2d_pca_mean_rao_values/sample_2d_pca_mean_rao_values.csv` | `scale`, `quad_id`, `basis`, `sample_max_points` | Appendable sampled 50 m comparison values. |
| Restored 10 m verification | `reports/tables/20261015_pca_mean_rao_completion/restored_10m_vector_normalized_comparison.csv` | `scale`, `quad_id`, `basis` | Independent rerun compared with saved values, counts, differences, and tolerance checks. Existing reference products are preserved. |
| Alpha = 0.5 small-scale run | `reports/tables/20261015_alpha_half_small_scales/` | `scale`, `quad_id`, `basis`, `alpha_value` | Resumable per-quadrat checkpoints for the authorized 10 m, 20 m, and 25 m run. 50 m is outside scope. |
| Alpha = 0.5, 10 m | `Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_10m_all_pixels.csv` | `scale`, `quad_id`, `basis`, `alpha_value` | 4,000 full-grid/basis rows. Finite areas: original 1,743; vector-normalized 1,744. Original-PCA 1000_c remains NA; 1402_b, 1019_a, and 1816_c have user-reviewed manual orientation repairs with explicit provenance. |
| Alpha = 0.5, 20 m | `Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_20m_all_pixels.csv` | same | 1,000 full-grid/basis rows. Finite areas: original 436; vector-normalized 436. Original-PCA 413 includes a user-authorized manual arc-orientation repair, explicitly identified in `source` and `warning`. |
| Alpha = 0.5, 25 m | `Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_25m_all_pixels.csv` | same | 640 full-grid/basis rows; 312 finite areas per basis and eight missing-raster rows per basis. No exclusions applied. |

The current small-scale hull workflow writes canonical long-format tables named
`Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_*m_all_pixels.csv` when each
scale finishes. These retain both PCA bases, all boundary IDs, explicit statuses,
pixel/unique-point counts, runtimes, and `source` provenance distinguishing reused
finite historical values from newly computed areas. They do not replace the
older alpha = 1 fields. Status summaries and unresolved cases are recorded in the
dated run folder. No 25 m exclusion flags are added.

Thirteen initially failed areas were recovered through exact integrals of the
same returned circular-arc boundary. `source` distinguishes center/arc and
endpoint/chord recovery; the latter verifies exact oriented vertex-ID balance
and endpoint-circle consistency. Only original-PCA 10 m 1000_c remains an explicit geometry-failure NA after the reviewed manual repairs.
First attempts, numerical warnings, validation tables, coverage, and a recovery
audit are retained in the dated run directory. Native triangulation jitter
warnings appear on three successful rows; no manual jitter was introduced.

20 m original-PCA quadrat 413 was subsequently repaired by user authorization:
the traversal of arc 423 in an isolated three-arc loop was reversed, preserving
the input points, alpha, and circular-arc coordinates. Its area is
460.7109822782754 PCA-score squared units. Source suffix
`manual_arc_orientation_repair` distinguishes it from unmodified package output.
The original NA, corrected signed traversal, affected-table backups, audit, and
validation are under `reports/tables/20261015_alpha_half_small_scales/manual_orientation_20m_413/`.
Earlier validation tables describe their pre-manual-repair snapshots; the manual
repair's `validation.csv` records the subsequent focused checks.

The user subsequently authorized the reviewed arc-orientation repairs of 10 m
1402_b, 1019_a, and 1816_c. They use the same manual source suffix and retain
all original point/curve geometry. Audit, signed traversals, before-table
backups, and validation are under
`reports/tables/20261015_alpha_half_small_scales/manual_orientation_three_10m/`.
The user deferred 1000_c for closer inspection; its NA is unchanged.

The current priority fields in the canonical combined spectral summaries and joined analysis tables are `pca_euclidean_mean` / `spec_pca_mean`, `alpha_hull_area` / `spec_alpha`, `standardized_PCA_pca_euclidean_mean` / `spec_spca_mean`, and `standardized_PCA_alpha_hull_area` / `spec_spca_alpha`.

The October 15 completion task also includes `rao_q_pca` / `spec_rao_q` and
`standardized_PCA_rao_q_pca` / `spec_spca_rao`. These use all retained pixels in
PC1-PC3 with equal weights and squared Euclidean distance, computed as
`2 * mean(squared_radius)`. They are distinct from the sampled PC1-PC2 products.

The 25 m component tables are `Quad_Values/25m_PCA_metrics_smooth_masked_5nm_summary.csv`
and `Quad_Values/25m_standardized_PCA_metrics_smooth_masked_5nm_summary.csv`. They
have 320 rows each and intentionally contain only distance/Rao and processing
metadata. Missing-raster rows have NA metric values; no manual-exclusion fields
or flags are added. There is not yet a joined `quadrat_analysis_25m.csv`.

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
