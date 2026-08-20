# Task Report: Spectral Heterogeneity Metric Relationship Analysis

Date: 2026-08-17

## Objective

Create all-pair spectral heterogeneity scatterplot sets and matching elevation, species-composition, individual-count, species-count, whole-brightness, and region-specific pixel-brightness highlight figures.

## Outputs Created

- `reports/analysis/20260817_spectral_heterogeneity_relationship_analysis.md`
- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_pairwise_correlations.csv`
- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_elevation_adjusted_models.csv`
- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_summary.csv`
- `reports/tables/spectral_heterogeneity_relationships/spectral_composition_scatterplot_silhouette.csv`
- `reports/tables/spectral_heterogeneity_relationships/quadrat_pixel_brightness_summary.csv`
- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_regional_illumination_correlations.csv`
- `reports/figures/spectral_heterogeneity_relationships/`

## Result Size

- Spectral pairwise correlation rows: 234
- Elevation-adjusted model rows: 234
- Spectral composition silhouette rows: 1170

## Notes

- Spectral metric pairings use complete cases for each pair, preserving upstream missingness.
- Pixel brightness is calculated from each quadrat raster as mean retained sunlit reflectance at 563 nm, the same band used for shadow masking.
- Region-specific brightness uses retained sunlit pixel means for blue 450-495 nm, green 500-570 nm, red 620-750 nm, and near-infrared 750-998 nm.
- Composition types and quadrat species/individual counts are reused from the presence-based biodiversity composition analysis.
