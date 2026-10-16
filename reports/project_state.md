# Project State

Last updated: 2026-10-15

## Current Objective

Maintain a reproducible spectral-biodiversity workflow at 10 m, 20 m, 25 m, and 50 m. Mean distance and Rao's Q are complete for available target inputs. The authorized alpha = 0.5 run at 10 m, 20 m, and 25 m is finished; only original-PCA 10 m 1000_c remains unresolved after three reviewed manual repairs. User-authorized manual arc-orientation repair resolved 20 m quadrat 413. Processing stopped before 50 m. 33 m and 100 m are deferred.

## Active Analysis

- Current spectral products: original-PCA and vector-normalized PCA all-pixel alpha-hull area in PC1-PC2 at alpha = 0.5 for 10 m, 20 m, and 25 m, plus PC1-PC3 mean distance and Rao's Q. New alpha-only tables are canonical for alpha = 0.5; older combined-table alpha fields retain alpha = 1.
- Biodiversity metrics: richness, Shannon, Simpson, evenness, Faith PD, phylogenetic Rao's Q, and abundance-weighted Faith PD.
- Canonical spectral summaries: `Quad_Values/Spectral_diversitySHPs/`.
- Canonical biodiversity summaries: `Quad_Values/Diversity_SHPs/`.
- Current joined inputs: `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`.

## Recent Direction Changes

- The workflow now prioritizes both original-PCA and vector-normalized PCA distance and alpha-hull area over spectral angle entropy, convex hulls, and 3D hull products.
- Historical metric products are retained as compact provenance, not deleted by default.
- Intermediate `*_smooth/` spectra were intentionally removed to conserve space because they can be regenerated.
- `10m_smooth_5nm` has been restored from a data repository. All 1,744 non-excluded reference quadrats reproduce the saved vector-normalized mean distance, Rao's Q, and pixel counts exactly at CSV precision.
- 25 m is the active scale addition, with 312 prepared rasters for 320 grid cells. The user will handle exclusions manually; no new exclusion flags are added.
- 33 m and 100 m are deferred; retain their resources.
- Repository consolidation must not alter `Documents/`, manuscript-state assets, or PDFs.

## Completed Work

- Established current diversity, environmental, and spectral output locations.
- Generated and documented 10 m, 20 m, and 50 m biodiversity and environmental summaries.
- Produced global original and vector-normalized PCA bases and spectral-summary tables.
- Created joined analysis tables for the established scales.
- Completed historical sensitivity, bootstrap, correlation, regression, and alpha-hull exploratory work retained under `reports/`.
- Refreshed governance, research objectives, directory mapping, and the data dictionary to distinguish current products from archived provenance.
- Verified the restored 10 m input against every eligible saved vector-normalized mean-distance and Rao value; preserved existing summaries and both PCA bases.
- Completed original and vector-normalized 25 m mean distance and Rao's Q for all 312 prepared rasters. The canonical `spectral_mean_rao_25m_smooth_masked_5nm_summary.csv` has 320 grid rows, including eight explicit missing-raster rows. No new exclusion flags were added.
- Finished the alpha = 0.5 production passes at 10 m, 20 m, and 25 m (approximately 31, 20, and 49 minutes). Canonical alpha-only tables cover their full boundary grids. Vector-normalized finite counts are 1,744 / 436 / 312; original-PCA finite counts are 1,743 / 436 / 312 after the user-authorized manual repair of 20 m 413 (area 460.7109822782754). Its original NA and signed traversal audit are preserved in `reports/tables/20261015_alpha_half_small_scales/manual_orientation_20m_413/`.
- Recovered 13 of 18 initial geometry/area failures using validated exact integrals of returned hull arcs, preserving successful package values and first-attempt provenance. Both integral forms were tested on 100 synthetic cases; actual-data checks and bounds also passed.

## Pending Work

- Inspect original-PCA 10 m `1000_c` later as requested by the user; its NA is unchanged. User-reviewed manual repairs resolved `1402_b`, `1019_a`, and `1816_c`. Their original failures, corrected traversals, and validation are in `reports/tables/20261015_alpha_half_small_scales/manual_orientation_three_10m/`. These and the earlier 20 m `413` repair are specific reviewed exceptions. No new exclusions added. See the October 15 small-scale hull report and `unresolved_geometry.csv`.
- 50 m original-PCA alpha = 0.5 values remain future work; no 50 m processing was performed or initiated in this task.
- 25 m diversity, environmental outputs, and the joined analysis table remain future work.
- Eight 25 m grid cells lack prepared rasters; existing 10 m and 20 m spectral footprints also omit 91 and 15 full-grid cells, respectively. Do not invent values for absent rasters.
- Decide whether the three root-level duplicated combined spectral-summary CSVs should be removed after downstream reference checks.
- Decide which compact archive products are worth retaining once 25 m and 33 m outputs exist.
- Apply spatial/environmental sensitivity work before final inferential claims where needed.

## Known Issues and Decisions Needed

- `Quad_Values/` contains historical 20 m-only outputs that are not current inputs but may be useful provenance. The earlier `Spectral_diversitySHPs_first/` archive has been removed; current canonical spectral products remain in `Spectral_diversitySHPs/`.
- The root-level and `Spectral_diversitySHPs/` combined spectral-summary CSVs are exact duplicates; the latter is canonical.
- Preserve both PCA `.rds` basis objects and their variance-explained tables. Future modeling should explicitly select `spec_pca_mean`, `spec_alpha`, `spec_spca_mean`, and `spec_spca_alpha`.
- Current alpha = 0.5 values live in `alpha_hull_0_5_*m_all_pixels.csv`; existing combined-table alpha fields retain alpha = 1. Three successful rows retain library automatic-jitter warnings for manual review: original 10 m `413_d`, and both bases for 25 m `sub25_69`.
- No byte-identical figures were found in `reports/figures/`; do not remove figures based solely on topic or filename similarity.
