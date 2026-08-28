# MDPI Remote Sensing Figure and Table Package

Tables use manuscript-facing terminology from SVH_Paper.docx. Internal variable names are included only where they help connect the table to analysis outputs.

## Table 1. Metric definitions and calculation roles

| Metric family | Manuscript label | Variable | Definition / role |
| --- | --- | --- | --- |
| Spectral | standardized PCA alpha hull | spec_spca_alpha | Area occupied by quadrat pixels in standardized PCA spectral space; emphasizes the spectral envelope. |
| Spectral | standardized PCA mean distance | spec_spca_mean | Mean Euclidean distance of quadrat pixels from the quadrat spectral centroid in standardized PCA space; emphasizes typical dispersion. |
| Spectral | standardized PCA Rao's Q | spec_spca_rao | Rao-style spectral dissimilarity in standardized PCA space; emphasizes pairwise squared spectral distances. |
| Biodiversity | abundance-weighted Faith's PD | phy_afaith | Phylogenetic diversity weighted by canopy crown-overlap abundance; primary phylogenetic response emphasized in the abstract. |
| Biodiversity | phylogenetic Rao's Q | phy_rao | Abundance-weighted phylogenetic dissimilarity using pairwise cophenetic distances among species. |
| Biodiversity | Shannon diversity | sp_shannon | Conventional species-diversity contrast based on proportional crown-overlap abundance. |


## Table 2. Quadrat coverage and complete-case counts

| Scale | Spectral-biodiversity complete cases (n) | Primary spectral metrics | Primary biodiversity contrast |
| --- | --- | --- | --- |
| 10 m | 1744 | standardized PCA alpha hull; standardized PCA mean distance; standardized PCA Rao's Q | abundance-weighted Faith's PD; phylogenetic Rao's Q; Shannon diversity |
| 20 m | 436 | standardized PCA alpha hull; standardized PCA mean distance; standardized PCA Rao's Q | abundance-weighted Faith's PD; phylogenetic Rao's Q; Shannon diversity |
| 50 m | 74 | standardized PCA alpha hull; standardized PCA mean distance; standardized PCA Rao's Q | abundance-weighted Faith's PD; phylogenetic Rao's Q; Shannon diversity |


## Table 3. Top standardized PCA spectral-biodiversity correlations by scale

| Scale | Biodiversity metric | Spectral metric | n | r | R2 | p |
| --- | --- | --- | --- | --- | --- | --- |
| 10 m | abundance-weighted Faith's PD | standardized PCA alpha hull | 1744 | 0.188 | 0.035 | <0.001 |
| 10 m | phylogenetic Rao's Q | standardized PCA alpha hull | 1744 | 0.177 | 0.031 | <0.001 |
| 10 m | Faith's PD | standardized PCA alpha hull | 1744 | 0.153 | 0.023 | <0.001 |
| 10 m | phylogenetic Rao's Q | standardized PCA mean distance | 1744 | 0.140 | 0.020 | <0.001 |
| 10 m | abundance-weighted Faith's PD | standardized PCA mean distance | 1744 | 0.122 | 0.015 | <0.001 |
| 10 m | Faith's PD | standardized PCA mean distance | 1744 | 0.119 | 0.014 | <0.001 |
| 20 m | abundance-weighted Faith's PD | standardized PCA alpha hull | 436 | 0.384 | 0.147 | <0.001 |
| 20 m | phylogenetic Rao's Q | standardized PCA alpha hull | 436 | 0.336 | 0.113 | <0.001 |
| 20 m | abundance-weighted Faith's PD | standardized PCA mean distance | 436 | 0.301 | 0.091 | <0.001 |
| 20 m | Faith's PD | standardized PCA alpha hull | 436 | 0.285 | 0.081 | <0.001 |
| 20 m | phylogenetic Rao's Q | standardized PCA mean distance | 436 | 0.274 | 0.075 | <0.001 |
| 20 m | Faith's PD | standardized PCA mean distance | 436 | 0.240 | 0.058 | <0.001 |
| 50 m | abundance-weighted Faith's PD | standardized PCA alpha hull | 74 | 0.497 | 0.247 | <0.001 |
| 50 m | phylogenetic Rao's Q | standardized PCA alpha hull | 74 | 0.487 | 0.238 | <0.001 |
| 50 m | phylogenetic Rao's Q | standardized PCA mean distance | 74 | 0.444 | 0.197 | <0.001 |
| 50 m | abundance-weighted Faith's PD | standardized PCA mean distance | 74 | 0.419 | 0.175 | <0.001 |
| 50 m | phylogenetic Rao's Q | standardized PCA Rao's Q | 74 | 0.352 | 0.124 | 0.003 |
| 50 m | Faith's PD | standardized PCA alpha hull | 74 | 0.312 | 0.097 | 0.013 |


## Table 4. Biodiversity metric concordance by scale

| Scale | Metric pair | n | r | R2 | p |
| --- | --- | --- | --- | --- | --- |
| 10 m | Faith's PD vs. species richness | 2000 | 0.710 | 0.504 | <0.001 |
| 10 m | phylogenetic Rao's Q vs. abundance-weighted Faith's PD | 2000 | 0.750 | 0.562 | <0.001 |
| 10 m | abundance-weighted Faith's PD vs. Shannon diversity | 2000 | 0.426 | 0.182 | <0.001 |
| 10 m | species richness vs. Shannon diversity | 2000 | 0.824 | 0.679 | <0.001 |
| 10 m | Shannon diversity vs. Simpson diversity | 2000 | 0.969 | 0.939 | <0.001 |
| 20 m | Faith's PD vs. species richness | 500 | 0.727 | 0.528 | <0.001 |
| 20 m | phylogenetic Rao's Q vs. abundance-weighted Faith's PD | 500 | 0.742 | 0.551 | <0.001 |
| 20 m | abundance-weighted Faith's PD vs. Shannon diversity | 500 | 0.152 | 0.023 | <0.001 |
| 20 m | species richness vs. Shannon diversity | 500 | 0.856 | 0.732 | <0.001 |
| 20 m | Shannon diversity vs. Simpson diversity | 500 | 0.949 | 0.901 | <0.001 |
| 50 m | Faith's PD vs. species richness | 80 | 0.761 | 0.580 | <0.001 |
| 50 m | phylogenetic Rao's Q vs. abundance-weighted Faith's PD | 80 | 0.688 | 0.474 | <0.001 |
| 50 m | abundance-weighted Faith's PD vs. Shannon diversity | 80 | -0.120 | 0.014 | 0.293 |
| 50 m | species richness vs. Shannon diversity | 80 | 0.847 | 0.718 | <0.001 |
| 50 m | Shannon diversity vs. Simpson diversity | 80 | 0.939 | 0.883 | <0.001 |


## Table 5. Final main-versus-supplement figure inventory

| Figure | File | Recommended placement | Purpose |
| --- | --- | --- | --- |
| 1 | 01_study_site_quadrat_scales_workflow.png | Main | Introduces the site, shared 10, 20, and 50 m quadrat grains, and the analysis workflow. |
| 2 | 02_standardized_pca_spectral_metric_relationships_by_scale.png | Main | Shows how standardized PCA spectral heterogeneity metrics relate to one another across scale. |
| 3 | 03_biodiversity_metric_relationships_by_scale.png | Main | Shows concordance and divergence among species and phylogenetic biodiversity metrics. |
| 4 | 04_standardized_pca_spectral_biodiversity_correlations_by_scale.png | Main | Summarizes all primary standardized PCA spectral-biodiversity correlations by scale. |
| 5 | 05_abundance_weighted_faith_pd_scale_response.png | Main | Highlights the scale-dependent relationship with abundance-weighted Faith's PD. |
| S1 | 06_raw_vs_vector_normalized_pca_brightness_relationship.png | Supplement | Documents brightness dominance in raw PCA and reduction after vector normalization. |
| S2 | 07_spatial_moran_i_diagnostics.png | Supplement | Summarizes Moran's I diagnostics for priority variables and residuals. |
| S3 | 08_transformation_sensitivity_priority_relationships.png | Supplement | Shows that transformations mainly affected weaker priority relationships. |
