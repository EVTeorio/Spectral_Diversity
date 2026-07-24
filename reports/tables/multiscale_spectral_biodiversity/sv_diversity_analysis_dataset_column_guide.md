# SV Diversity Analysis Dataset Column Guide

Last updated: 2026-07-24

This guide explains the columns in `sv_diversity_analysis_dataset.csv` in plain language. The dataset combines quadrat identity, plant diversity, phylogenetic diversity, spectral heterogeneity, elevation/topographic variables, and analysis flags used in the multiscale spectral-diversity workflow.

## Quadrat Data

Columns 1-4 identify each quadrat and locate it in space.

| Column number | Column | Plain-language meaning |
|---:|---|---|
| 1 | `quad_id` | The quadrat identifier used to match records across spectral, plant-diversity, phylogenetic-diversity, and environmental tables. |
| 2 | `scale` | The quadrat size or spatial grain: `10m`, `20m`, or `50m`. |
| 3 | `center_x` | The x coordinate of the quadrat center, calculated from the quadrat polygon centroid. |
| 4 | `center_y` | The y coordinate of the quadrat center, calculated from the quadrat polygon centroid. |

## Species And Phylogenetic Diversity

Columns 5-11 describe plant species diversity and phylogenetic diversity within each quadrat.

| Column number | Column | Plain-language meaning |
|---:|---|---|
| 5 | `sp_rich` | Species richness: the number of species present in the quadrat. |
| 6 | `sp_shannon` | Shannon species diversity, which increases when there are more species and when abundance is more evenly spread among species. |
| 7 | `sp_simpson` | Simpson species diversity, calculated as 1 minus the summed squared species proportions. Higher values mean the quadrat is less dominated by one or a few species. |
| 8 | `sp_even` | Species evenness: how evenly abundance is distributed among the species present. |
| 9 | `phy_faith` | Faith's phylogenetic diversity: the total phylogenetic branch length represented by the species present in the quadrat. |
| 10 | `phy_rao` | Phylogenetic Rao's Q: a phylogenetic diversity value that accounts for both species abundance and phylogenetic distances among species. |
| 11 | `phy_afaith` | Abundance-weighted Faith's phylogenetic diversity: Faith's PD adjusted by the relative abundance or crown-overlap contribution of species. |

## Spectral Heterogeneity Measures

Columns 12-28 describe spectral variation within each quadrat. These are the main hyperspectral variation metrics.

| Column number | Column | Plain-language meaning |
|---:|---|---|
| 12 | `spec_sa` | Spectral angle entropy. This is the primary spectral-angle heterogeneity value from sunlit, shadow-masked, smoothed 5 nm spectra. For most quadrats it is the mean of bootstrap iterations; for small enough quadrats it may be exact all-pixel entropy. |
| 13 | `spec_pca_mean` | Mean Euclidean distance of retained pixels from their quadrat centroid in raw global PC1-PC3 spectral space. Larger values mean pixels are more spread out spectrally. |
| 14 | `spec_pca_med` | Median Euclidean distance of retained pixels from their quadrat centroid in raw global PC1-PC3 spectral space. |
| 15 | `spec_pca_sd` | Standard deviation of the raw PCA Euclidean distances. This describes how variable the pixel-to-centroid distances are within the quadrat. |
| 16 | `spec_rao_q` | Spectral Rao's Q calculated in raw global PC1-PC3 space using equal pixel weights and squared Euclidean distance. |
| 17 | `spec_alpha` | Alpha-hull area of the quadrat's pixels in raw global PC1-PC2 spectral space. This describes the two-dimensional area occupied by the pixel cloud. |
| 18 | `spec_convex` | Convex-hull area of the quadrat's pixels in raw global PC1-PC2 spectral space. This is a supplemental hull summary. |
| 19 | `spec_hull3d_v` | Convex-hull volume of the quadrat's pixels in raw global PC1-PC3 spectral space. This is a supplemental three-axis hull summary. |
| 20 | `spec_hull3d_a` | Convex-hull surface area of the quadrat's pixels in raw global PC1-PC3 spectral space. This is a supplemental three-axis hull summary. |
| 21 | `spec_spca_mean` | Mean Euclidean distance of retained vector-normalized spectra from their quadrat centroid in standardized PCA PC1-PC3 spectral space. This is the current primary PCA-based spectral variation measure. |
| 22 | `spec_spca_med` | Median Euclidean distance of retained vector-normalized spectra from their quadrat centroid in standardized PCA PC1-PC3 spectral space. |
| 23 | `spec_spca_sd` | Standard deviation of standardized PCA Euclidean distances. This describes how variable the standardized pixel-to-centroid distances are within the quadrat. |
| 24 | `spec_spca_rao` | Spectral Rao's Q calculated in standardized PCA PC1-PC3 space using equal pixel weights and squared Euclidean distance. |
| 25 | `spec_spca_alpha` | Alpha-hull area of vector-normalized spectra in standardized PCA PC1-PC2 spectral space. |
| 26 | `spec_spca_convex` | Convex-hull area of vector-normalized spectra in standardized PCA PC1-PC2 spectral space. This is a supplemental hull summary. |
| 27 | `spec_spca_hull3d_v` | Convex-hull volume of vector-normalized spectra in standardized PCA PC1-PC3 spectral space. This is a supplemental three-axis hull summary. |
| 28 | `spec_spca_hull3d_a` | Convex-hull surface area of vector-normalized spectra in standardized PCA PC1-PC3 spectral space. This is a supplemental three-axis hull summary. |

## Elevation And Topographic Measures

Columns 29-32 describe elevation and terrain roughness within each quadrat.

| Column number | Column | Plain-language meaning |
|---:|---|---|
| 29 | `env_elev` | Mean elevation within the quadrat, derived from the leaf-off digital terrain model. |
| 30 | `env_tri5` | Mean Riley topographic roughness index within the quadrat, calculated from a 5x5 moving window. |
| 31 | `env_tri11` | Mean Riley topographic roughness index within the quadrat, calculated from an 11x11 moving window. |
| 32 | `env_tri21` | Mean Riley topographic roughness index within the quadrat, calculated from a 21x21 moving window. |

## Other Related Info

Columns 33 through the end contain processing metadata and analysis flags.

| Column number | Column | Plain-language meaning |
|---:|---|---|
| 33 | `sa_n_pixels` | Number of retained sunlit pixels used or available for the spectral angle entropy calculation. |
| 34 | `sa_method` | Method used for the spectral angle entropy value, such as bootstrap mean, exact all-pixel entropy, or insufficient pixels. |
| 35 | `sa_all_pixels_sampled` | Indicates whether spectral angle entropy used all retained pixels. `TRUE` means all retained pixels were sampled; `FALSE` means a capped or sampled approach was used. |
| 36 | `spca_n_pixels` | Number of retained sunlit pixels used or available for the standardized PCA spectral heterogeneity calculation. |
| 37 | `spca_metric_method` | Method used for the standardized PCA metric, all retained pixels unless the quadrat was excluded or had too few valid pixels. |
| 38 | `spca_euclidean_all_pixels_sampled` | Indicates whether standardized PCA Euclidean distance used all retained pixels. |
| 39 | `parent_20m` | The parent 20 m quadrat ID. This is most useful for 10 m quadrats, which nest inside 20 m quadrats. |
| 40 | `edge_flag` | Indicates whether the quadrat is flagged as an edge quadrat for the analysis workflow. |
| 41 | `primary_analysis` | Indicates whether the row is included in the primary analysis dataset after applying the current analysis flags and exclusions. |
