# Combined Quadrat Analysis Tables Variable Guide

Last updated: 2026-06-24

This guide documents the shortened columns in the root-level combined CSVs:

- `quadrat_analysis_10m.csv`
- `quadrat_analysis_20m.csv`
- `quadrat_analysis_50m.csv`

Center coordinates are polygon centroids from the plant-diversity shapefiles. CRS: `NAD83 / UTM zone 16N`.

Pixel-count, pair-count, bootstrap replicate, method, exclusion, geometry metadata, and per-species composition fields are intentionally excluded from the combined CSVs.

## Column Definitions

| Column | Group | Source Name | Definition |
|---|---|---|---|
| `quad_id` | Identifier and coordinates | `quad_id` | Scale-aware quadrat identifier used for joins across spectral, diversity, and environmental products. |
| `scale` | Identifier and coordinates | `scale` | Quadrat grain: 10m, 20m, or 50m. |
| `center_x` | Identifier and coordinates | `centroid X` | X coordinate of the quadrat polygon centroid in the plant-diversity shapefile CRS. |
| `center_y` | Identifier and coordinates | `centroid Y` | Y coordinate of the quadrat polygon centroid in the plant-diversity shapefile CRS. |
| `spec_sa` | Spectral heterogeneity | `sa_entropy` | Primary spectral angle entropy value from sunlit, shadow-masked smoothed 5 nm spectra; exact all-pixel entropy where feasible, otherwise the bootstrap mean. |
| `spec_pca_mean` | Spectral heterogeneity | `pca_euclidean_mean` | Mean Euclidean distance of retained pixels from the quadrat centroid in global PC1-PC3 spectral space. |
| `spec_pca_med` | Spectral heterogeneity | `pca_euclidean_median` | Median Euclidean distance of retained pixels from the quadrat centroid in global PC1-PC3 spectral space. |
| `spec_pca_sd` | Spectral heterogeneity | `pca_euclidean_sd` | Standard deviation of Euclidean distances of retained pixels from the quadrat centroid in global PC1-PC3 spectral space. |
| `spec_rao_q` | Spectral heterogeneity | `rao_q_pca` | Rao's Q spectral heterogeneity metric using equal pixel weights and squared Euclidean distance in global PC1-PC3 space. |
| `spec_alpha` | Spectral heterogeneity | `alpha_hull_area` | Alpha-hull area in global PC1-PC2 spectral space. |
| `spec_convex` | Spectral heterogeneity | `pca_convex_hull_area` | Convex-hull area in global PC1-PC2 spectral space, used as a supplemental hull summary. |
| `spec_hull3d_v` | Spectral heterogeneity | `pca_hull_volume_3d` | Convex-hull volume in global PC1-PC3 spectral space, used as a supplemental three-axis hull summary. |
| `spec_hull3d_a` | Spectral heterogeneity | `pca_hull_area_3d` | Convex-hull surface area in global PC1-PC3 spectral space, used as a supplemental three-axis hull summary. |
| `sp_rich` | Species and phylogenetic diversity | `richness` | Number of species with positive crown-overlap values in the quadrat. |
| `sp_shannon` | Species and phylogenetic diversity | `shannon` | Shannon species diversity calculated from positive species crown-overlap proportions. |
| `sp_simpson` | Species and phylogenetic diversity | `simpson` | Simpson diversity calculated from positive species crown-overlap proportions as 1 minus the sum of squared proportional abundances. |
| `sp_even` | Species and phylogenetic diversity | `evenness` | Species evenness calculated as Shannon diversity divided by log richness when richness is greater than 1. |
| `phy_faith` | Species and phylogenetic diversity | `faith_pd` | Faith's phylogenetic diversity calculated as summed branch length for species present in the quadrat. |
| `phy_rao` | Species and phylogenetic diversity | `rao_pd` | Rao phylogenetic diversity calculated from species crown-overlap weights and phylogenetic distances. |
| `phy_afaith` | Species and phylogenetic diversity | `afaith_pd` | Abundance-weighted Faith's phylogenetic diversity calculated by weighting phylogenetic branch lengths by species crown-overlap values. |
| `env_elev` | Environmental/topographic | `elev_mean` | Mean DTM elevation within the quadrat, in meters. |
| `env_tri5` | Environmental/topographic | `tri5_mean` | Mean Riley topographic roughness index from a 5x5 DTM moving window. |
| `env_tri11` | Environmental/topographic | `tri11_mean` | Mean Riley topographic roughness index from an 11x11 DTM moving window. |
| `env_tri21` | Environmental/topographic | `tri21_mean` | Mean Riley topographic roughness index from a 21x21 DTM moving window. |
