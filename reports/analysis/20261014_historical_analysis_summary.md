# Historical Analysis Summary

Last updated: 2026-10-14

This summary replaces detailed pre-October analysis notes that are no longer the active direction. Underlying data tables, figures, scripts, original-PCA metrics, vector-normalized PCA metrics, and both PCA basis objects are retained unless separately noted.

## Workflow Foundation: June-July 2026

- Confirmed 10 m, 20 m, and 50 m quadrat spectra; created smoothed and 5 nm resampled derivatives.
- Calculated spectral-angle entropy with a 563 nm shadow mask. Most large quadrats use a 70-iteration bootstrap mean; this remains archive/QC information, not a current primary metric.
- Generated plant-diversity, phylogenetic-diversity, environmental, and joined quadrat tables for the three established scales.
- Rebuilt both PCA bases from the 10 m footprint after masking: original PCA and vector-normalized PCA. Both bases, their `.rds` files, variance tables, PCA distance, and alpha-hull metrics are retained.
- Established direct scale-specific spectral-biodiversity correlations. Relationships were generally stronger for phylogenetic than species-diversity measures and strongest at 50 m.

## Diagnostics and Direction Changes: July-August 2026

- Bootstrap, edge, brightness, sample-size, loading-region, and metric-relationship analyses documented sensitivity and interpretation limits.
- Original and vector-normalized PCA were retained as complementary bases. Convex hull, 3D hull, and spectral-angle entropy were moved out of the primary forward metric set.
- The research direction converged on PCA distance and alpha-hull area with current species and phylogenetic biodiversity metrics, while retaining the wider metric suite as provenance.
- Spatial autocorrelation and environmental context remain sensitivity considerations rather than grounds for causal interpretation.

## Current October Sequence

The following detailed analysis reports remain because they directly support the current PCA and alpha-hull direction:

- `20260819_final_research_direction_analysis.md`
- `20261009_pca_phylo_shannon_multiple_regression.md`
- `20261009_pca_phylogenetic_only_multiple_regression.md`
- `20261009_pca_sa_phylogenetic_only_multiple_regression.md`
- `20261010_alpha_hull_all_pixels_benchmark.md`
- `20261011_sample_2d_pca_mean_rao_50m.md`
- `20261012_alpha_hull_all_pixels_values_multiscale_alpha0_5_1_1_5.md`

The multiscale alpha-hull report supersedes the individual October 11 notes for separate scale/alpha runs and their processing-speed summaries.
